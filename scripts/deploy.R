# Generic Rcpp BiocStyle deployment script
#
# Usage:
#   #bash
#   $ R --slave < deploy.R
#  or
#   #R shell in package directory
#    > source("deploy.R")
#     > deployR()
#
# mammana, helmuth 2015-01-13

#action
deployR <- function(pckg=".") {
  chgDir <- F
  if (pckg == ".") {#assuming current directory
    pckg <- tail(strsplit(getwd(), split="/")[[1]], 1)
    chgDir <- T
    setwd("..")
  }

  if (!file.exists(file.path(pckg, "DESCRIPTION"))) {
    if (chgDir) {
      setwd(pckg)
    }
    stop(paste(pckg, "is not a valid R package path"))
    stop("Package directory does not exist.")
  }


  message("Deploying ", pckg, "...")

  vsion <- read.dcf(file.path(pckg, "DESCRIPTION"))[,"Version"]
  message("Version: ", vsion)

  library(Rcpp)
  compileAttributes(pckg)

  library(roxygen2)
  roxygenize(pckg, clean=T)

  setwd(pckg)
  require(Kmisc)
  Kmisc:::registerFunctions(pckg, prefix="")

  library(devtools)
  build_vignettes()
  test()

  setwd("../")


  #  try installing the package in a stub library
  tmp.lib <- paste0(pckg, "_BUILDtmp")
  dir.create(tmp.lib) 
  system(paste0("R3 CMD build ", pckg))
  install.packages(paste0(pckg, "_", vsion, ".tar.gz"), lib=tmp.lib)
  library(pckg, character.only=T)
  unlink(tmp.lib, recursive=T)


  if (chgDir) {
    setwd(pckg)
  }
}

#call
if (!interactive()) {
  args=(commandArgs(TRUE))
  if ( length(args) > 0 ) {
    for(i in 1:length(args)){
      eval(parse(text=args[[i]])) #reads origin
    }
    message(pckg)
    if (!exists("pckg"))
      stop("Provide a valid packagename with '--args pckg=\"packagename\"")
  } else {
    pckg="."
  }
  deployR(pckg)
}
