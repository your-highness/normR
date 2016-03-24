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

registerFunctions <- function (pkg_name, prefix="") {
  if (missing(pkg_name))
    pkg_name <- basename(getwd())
  files <- list.files("src", pattern = "[cC]$", full.names = TRUE)
  files <- files[files != paste0("src/", pkg_name, "_init.c")]
  c_files <- lapply(files, readLines)
  get_c_prototypes <- function(x) {
    export_lines <- grep("// \\[\\[export\\]\\]", x)
    sapply(export_lines, function(line) {
           return(gsub("(.*)// \\[\\[export\\]\\]\n(.*?) ?\\{(.*)",
                       "\\2;", paste(x[line:length(x)], collapse = "\n")))
  })}
  c_prototypes <- sapply(c_files, get_c_prototypes)
  if (length(c_prototypes) != 0)
    c_prototypes <- c_prototypes[sapply(c_prototypes, function(x) {
                                        !identical(x, list()) })]
  rcpp_exports <- readLines("src/RcppExports.cpp")
  fn_lines <- grep("^RcppExport", rcpp_exports, value = TRUE)
  cpp_prototypes <- sapply(fn_lines, USE.NAMES = FALSE, function(x) {
                           gsub("RcppExport (.*) \\{", "\\1;", x)  })
  all_prototypes <- unlist(c(c_prototypes, cpp_prototypes))
  all_names <- sapply(all_prototypes, function(x) {
                      gsub("SEXP (.*)\\(.*", "\\1", x)})
  all_nargs <- sapply(all_prototypes, function(x) {
                      defn <- gsub("SEXP (.*)\\((.*)\\).*", "\\2", x)
                      m <- gregexpr("SEXP +", defn)
                      if (identical(as.integer(m[[1]]), -1L)) {
                        return(0)
                      }
                      else {
                        return(length(m[[1]]))
                      }})
  Cnames <- paste0(prefix, all_names)
  cmd_lines <- paste0("{\"", Cnames, "\", (DL_FUNC) &", all_names,
                      ", ", all_nargs, "},")
  R_CallMethodsDef <- c("R_CallMethodDef callMethods[]  = {",
                        paste0("  ", cmd_lines), "  {NULL, NULL, 0}", "};")
  R_RegisterRoutines <- c(paste0("void R_init_", pkg_name,
                                 "(DllInfo *info) {"), "  R_registerRoutines(info, NULL, callMethods, NULL, NULL);",
                          "  R_useDynamicSymbols(info, FALSE);", "}")
  init.c <- c("#include <R.h>", "#include <Rinternals.h>",
              "", "#include <R_ext/Rdynload.h>", "", all_prototypes,
              "", R_CallMethodsDef, "", R_RegisterRoutines, "")
  cat(init.c, file = paste0("src/", pkg_name, "_init.c"), sep = "\n")
  return(invisible(NULL))
}

#action
deployR <- function(pckg=".", registerFun=T, check=T) {
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
  if (registerFun) {
    registerFunctions(pckg, prefix="")
  }

  library(devtools)
  build_vignettes()
  test()

  setwd("../")


  #  try installing the package in a stub library
  tmp.lib <- paste0(pckg, "_BUILDtmp")
  dir.create(tmp.lib) 
  if (check) {
    system(paste0("~/R-devel/bin/R CMD check ", pckg))
    system(paste0("~/R-devel/bin/R CMD BiocCheck ", pckg))
  }
  system(paste0("~/R-devel/bin/R CMD build ", pckg))
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
