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

  message(" #### Deploying ", pckg, "...")
  vsion <- read.dcf(file.path(pckg, "DESCRIPTION"))[,"Version"]
  message(" #### Version: ", vsion)

  #message(" #### devtools: Entering dev_mode()")
  #library(devtools, quietly=T)
  #dev_mode()
  #load_all()

  message(" #### Rcpp: compileAttributes(pckg) to generate object files")
  library(Rcpp, quietly=T)
  compileAttributes(pckg)

  message(" #### roxygen2: roxygenize(pckg, clean=T) to generate man/*.Rd")
  library(roxygen2, quietly=T)
  roxygenize(pckg, clean=T)

  message(paste0(" ### Kmisc: registerFunctions(pckg, prefix=\"\") to generate",
                 " src/", pckg, "_init.c"))
  setwd(pckg)
  if (registerFun) {
    unlink(paste0("src/", pckg, "_init.c"))
    registerFunctions(pckg, prefix="")
  }

  #message(" #### devtools: build_vignettes() to compile Rmarkdown")
  #build_vignettes()
  #message(" #### devtools: test() to run testthat routines")
  #test()
  setwd("../")

  #  try installing the package in a stub library
  message(" #### Trying to install the package")
  system(paste0("R CMD build ", pckg))
  install.packages(paste0(pckg, "_", vsion, ".tar.gz"), lib=tmp.lib)

  #extensive checking
  tmp.lib <- paste0(pckg, "_BUILDtmp")
  dir.create(tmp.lib)
  if (check) {
    message(" #### Testing requested: running 'R CMD check'")
    system(paste0("R CMD check ", pckg))
    message(" #### Testing requested: running 'R CMD BiocCheck'")
    system(paste0("R CMD BiocCheck ", pckg))
  }

  #message(" #### Trying to load the package")
  #library(pckg, character.only=T, quietly=T)
  #unlink(tmp.lib, recursive=T)

  if (chgDir) {
    setwd(pckg)
  }
  message(" #### SUCCESS! Refer to ../normr.Rcheck/00check.log to fix warnings")
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
