# Copyright (C) 2016 Johannes Helmuth
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#' Container for a NormR fit
#'
#' This S4 class wraps a NormR fit containing counts, fit configuration and 
#' results of the fit. 
#'
#' This class provides printing, summary and accessor methods for retrieving
#' data from the object. Internally this class uses a map structure which saves
#' memory, such that accessing data through \code{slot()} gives undesired
#' output. Please use accessor methods. A NormRFit object is read-only, i.e. it 
#' can not be modified.
#'
#' @param obj A NormRFit object.
#' @param i \code{integer()}-index for subsetting.
#' @param B The index of the background component.
#'
#' @slot type A \code{character()} representing the type of fit ("enrichR", 
#' "diffR", "regimeR")
#' @slot n Number of regions analyzed in total.
#' @slot ranges A \code{GenomicRanges}-object specifying the analysed regions.
#' @slot k An \code{integer()} giving the number of binomial mixture 
#' components to be fit to the data.
#' @slot B The index of the background component.
#' @slot map A \code{integer}-vector to map back \code{counts},
#' \code{lnposteriors}, \code{lnenrichment}, \code{lnpvals} and \code{lnqvals}. 
#' See also \code{\link{map2uniquePairs()}} for map generation.
#' @slot counts A \code{list()} of length 2 containing unique count pairs. Use 
#' accessor \code{\link{getCounts}} to retrieve original count matrix. 
#' @slot names The corresponding names for counts.
#' @slot thetaStar A \code{numeric()} representing a naive background 
#' estimation, i.e. \code{sum(getCounts(obj)[2,])/sum(getCounts(obj))}
#' @slot theta A \code{k}-length \code{numeric()}-vector representing 
#' NormR-fitted parametrization of binomial mixture components.
#' @slot mixtures A \code{k}-length \code{numeric()}-vector representing 
#' NormR-fitted mixture proportions of binomial mixture
#' components.
#' @slot lnL A \code{numeric()}-vector representing the log-likelihood-trace of
#' the NormR model fit.
#' @slot eps A \code{numeric{}} representing preset threshold for NormR fit 
#' convergence.
#' @slot lnposteriors A \code{k}-dimensional \code{matrix()} containing ln
#' posterior probabilities for every unique \code{n} 2-tupel in counts given 
#' each \code{theta} and \code{mixtures}. Use accessor
#' \code{\link{getPosterior()}} to get the complete posterior matrix.
#' @slot lnenrichment A \code{numeric()}-vector of normalized ln enrichment over 
#' component \code{B} for each unique count-tupel. Use accessor
#' \code{\link{getEnrichment}} to retrieve enrichment for original count sequence.
#' @slot lnpvals A \code{numeric()}-vector containing ln P-values every unique 
#' 2-tupel in counts given \code{theta} of \code{B}.
#' @slot filteredT A \code{integer()}-vector giving indices of P-values to be
#' considered for FDR correction. These indices have been obtainted by filtering
#' P-values using the T method.
#' @slot lnqvals A \code{numeric()}-vector containing ln q-values. These are 
#' P-values corrected for multiple testing using Storey's method.
#' @aliases NormRFit
#' @seealso \code{\link{normr-methods}} for the functions that produce 
#' this object
#' @return return values are described in the Methods section.
#' @example inst/examples/class_example.R
#' @export
setClass("NormRFit", 
    representation = representation(type = "character",
                                    n = "integer",
                                    ranges = "GenomicRanges",
                                    k = "integer",
                                    B = "integer",
                                    map = "integer",
                                    counts = "list",
                                    names = "character",
                                    thetastar = "numeric",
                                    theta = "numeric",
                                    mixtures = "numeric",
                                    lnL = "numeric",
                                    eps = "numeric",
                                    lnposteriors = "matrix",
                                    lnenrichment = "numeric",
                                    lnpvals = "numeric",
                                    filteredT = "numeric",
                                    lnqvals = "numeric" )
)

setValidity("NormRFit",
    function(object) {
      if (!(object@type %in% c("enrichR", "diffR", "regimeR"))) {
          return("invalid type slot")
      }
      if (object@n != length(object@ranges)) return("invalid n and ranges")
      if (length(object@k) == 0 | object@k <= 0) return("invalid k slot")
      if (object@B <= 0 || object@B > k) return("invalid B slot")
      if (length(object@map) != length(object@counts[[1]])) {
        return("incompatible map and count slot")
      }
      if (length(object@counts) != 2 ||
          length(object@counts[[1]]) != length(object@counts[[2]])) {
        return("invalid counts slot")
      }
      if (length(object@thetastar) != 1 || object@thetastar < 0 || 
          object@thetastar > 1) {
        return("invalid thetastar slot")
      }
      if (length(object@theta) != object@k) return("invald theta slot")
      if (length(object@mixtures) != object@k) return("invald mixtures slot")
      n <- length(object@counts[[1]])
      if (NCOL(object@lnposteriors) != object$k || 
          NROW(object@lnposteriors) != length(object@counts[[1]])) {
        return("invald lnposterios slot")
      }
      if (length(object@lnenrichment) != length(object@counts[[1]])) {
        return("invaled lnenrichment slot")
      }
      if (length(object@lnpvals) != length(object@counts[[1]])) { 
          return("invaled lnpvals slot")
      }
      if (max(object@filteredT) > length(object@counts[[1]])) {
        return("invaled filteredT slot")
      }
      if (length(object@lnqvals) != length(object@filteredT)) {
        return("invaled lnqvals slot")
      }
      TRUE
    }
)

setMethod("print", "NormRFit",
  function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("NormRFit-class xect\n\n",
        "Type:\t\t'", x@type, "'\n",
        "Number of Regions:\t", x@n, "\n",
        "Theta* (naive bg):\t", x@thetastar, "\n\n")
    if (length(x@theta)) {
      cat("Results of fit\n", 
          "Mixture Proporitons:\n")
      print.default(format(x@mixtures,digits=digits), print.gap=2L, quote=F)
      cat("Theta:\n")
      print.default(format(x@theta,digits=digits), print.gap=2L, quote=F)
    } else {
      cat("No results of fit.")
    }
    cat("\n")
    invisible(x)
  }
) 
setMethod("show", "NormRFit", function(object) print(object))

setMethod("plot", "NormRFit",
     function(x, ...) {
       stop("not implemented yet")
     }
)

#' @describeIn NormRFit Number of regions analyzed.
#' @aliases length
#' @export
setMethod("length", "NormRFit", function(x) x@n)

#' @describeIn NormRFit Prints a summary of the NormRFit object.
#' @param object A NormRFit objectect
#' @param ... Not used
#' @return NULL
#' @export
setMethod("summary", "NormRFit", 
  function(object, print=T, ...) {
    ans <- paste0("NormRFit-class objectect\n\n",
                  "Type:\t\t'", object@type, "'\n",
                  "Number of Regions:\t", object@n, "\n",
                  "Number of components:\t", object@k, "\n",
                  "Theta* (naive bg):\t", object@thetastar, "\n",
                  "Backgroundcomponent B:\t", object@B, "\n\n")
    if (length(object@theta)) {
      ans <- paste(ans, 
                   "Results of fit\n", 
                   "Mixture Proporitons:\n")
      print.default(format(object@mixtures,digits=digits), print.gap=2L, quote=F)
      cat("Theta:\n")
      print.default(format(object@theta,digits=digits), print.gap=2L, quote=F)
      cat("Bayesian Information Criterion:\n")
      print.default(format(
                           (-2 * object@lnL[length(object@lnL)] + length(object@theta) * log(object@n)),
                           quote=F))
      cat("Significantly different from background B:\n")
      cts <- c("***"=length(which(object@lnqvals <= log(0))),
               "*****" =length(which(object@lnqvals <= log(0.001))),
               "*"  =length(which(object@lnqvals <= log(0.01))),
               "."  =length(which(object@lnqvals <= log(0.05))),
               " "  =length(which(object@lnqvals <= log(0.1))))
      print.default(format(cts, digits=digits), print.gap=2L, quote=F)
      cat("---\nSignif. codes:  ", sleg, sep = "", fill = w + 
          4 + max(nchar(sleg, "bytes") - nchar(sleg)))
    } else {
      cat("No results of fit.")
    }
    cat("\n")
  }
  )

#'@export
setGeneric("getCounts", function(obj) standardGeneric("getCounts"))
#' @describeIn NormRFit Retrieve used counts.
#' @aliases getCounts
#' @export
setMethod("getCounts", "NormRFit", function(obj) { 
    list("control"=obj@counts[[1]][obj@map], 
         "treatment"=obj@counts[[2]][obj@map])
})

#'@export
setGeneric("getEnrichment", function(obj) standardGeneric("getEnrichment"))
#' @describeIn NormRFit Retrieve NormR calculated enrichment. Ln enrichment for 
#' a difference call.
#' @aliases getEnrichment
#' @export
setMethod("getEnrichment", "NormRFit", function(obj) {
    if (obj@tupe == "diffR") obj@lnenrichment[obj@map]
    else exp(obj@lnenrichment)[obj@map]
})

#'@export
setGeneric("getPosterior", function(obj) standardGeneric("getPosterior"))
#' @describeIn NormRFit Retrieve computed posteriors.
#' @aliases getPosterior
#' @export
setMethod("getPosterior", "NormRFit", function(obj) {
    exp(obj@lnposteriors)[obj@map,]
})

#'@export
setGeneric("getPvalues", function(obj, ...) standardGeneric("getPvalues"))
#' @describeIn NormRFit Retrieve computed P-values
#' @aliases getPosterior
#' @export
setMethod("getPvalues", "NormRFit", function(obj, filtered=F) {
  if (filtered) {
    idx = obj@map[which(obj@map %in% obj@filteredT)]
    exp(obj@lnpvals)[idx]
  } else {
    exp(obj@lnpvals)[obj@map]
  }
})

#'@export
setGeneric("getPvalues", function(obj) standardGeneric("getPvalues"))
#' @describeIn NormRFit Retrieve computed P-values
#' @aliases getPosterior
#' @export
setMethod("getPvalues", "NormRFit", function(obj) {
    exp(obj@lnpvals)[obj@map]
})

#'@export
setGeneric("getQvalues", function(obj) standardGeneric("getQvalues"))
#' @describeIn NormRFit Retrieve computed posteriors.
#' @aliases getPosterior
#' @export
setMethod("getQvalues", "NormRFit", function(obj) {
    exp(obj@lnqvals)
})

#' @export
setGeneric("writeEnrichment", function(obj,...) 
           standardGeneric("writeEnrichment"))
#' describeIn NormRFit Writes the calculated enrichment of a NormRFit to disk.
#' @param obj A NormRFit object
#' 
#' @aliases writeEnrichment
#' @export
setMethod("writeEnrichment", "NormRFit",
    function(obj, gr=NULL, filename=NULL) {
      if (is.null(gr) & is.null(obj@ranges)) stop("no chromosomes stored in obj")
      if (is.null(gr)) gr <- obj@ranges
      if (length(ranges) != obj@n) stop("gr not corresponding to obj")
      enr <- obj@enrichment
      if (is.null(enr)) enr <- getEnrichment(obj)
      require(rtracklayer)
      gr$score <- enr
      file <- filename
      if (grep('.bw$|.bigWig$', file, perl=T)) file.type="BigWig"
      else if (grep('.bg|bedGraph', file, perl=T)) file.type="bedGraph"
      else stop("filetype could not be determined by filename suffix")
      export(gr, filename, con=file.type)
    }
)

