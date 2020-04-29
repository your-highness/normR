#
# Automatic code generation (for package maintainers)
#
# You need to call this script if you change:
# 1. Rcpp functions (add, remove or change C function with \\[[Rcpp::export]])
# 2. Roxygen documentation (including the roxygen import statements)
# 
# R --slave < normr/scripts/generateCode.R
#
# Followed by
# R CMD build --keep-empty-dirs --no-resave-data normr 
# R CMD check --library=/home/helmuth/miniconda3/envs/R-4.0.0/lib/R/library \
#   --no-vignettes --timings \
#   normr_${VERSION}.tar.gz

load_or_get <- function(pkg){
  if (!require(pkg, character.only=T)) install.packages(pkg, dependencies=TRUE)
  library(pkg, character.only=T)
}

#wrap C++ functions into R functions
pckg <- "normr"
load_or_get("Rcpp")
compileAttributes(pckg)

#create man pages and the NAMESPACE file
load_or_get("devtools")
devtools::install_deps(pckg, dependencies=TRUE)
load_or_get("roxygen2")
roxygenize(pckg)
