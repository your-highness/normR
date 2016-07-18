#~/R-devel/bin/R

#CRAN packages
install.packages(c("testthat", "Kmisc", "devtools", "Rcpp", "roxygen2", "knitr", "rmarkdown"))

#github packages
devtools::install_github("Kmisc", "kevinushey")

#addtional bioconductor packages
source("https://bioconductor.org/biocLite.R")
#bamsignals
biocLite(c("BiocStyle", "bamsignals", "Rsamtools", "BiocGenerics", "zlibbioc", "Rhtslib"))
#normR
biocLite(c("GenomeInfoDb", "GenomicRanges", "IRanges", "BiocParallel"))


#deployment in current directory (building, registering, testing, packaging)
source("~/software/R/deploy.R")
deployR()
