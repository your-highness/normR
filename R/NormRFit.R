# Copyright (C) 2017 Johannes Helmuth & Ho-Ryun Chung
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
#' Container for a fit done with normR
#'
#' This S4 class wraps a \code{\link{normR}} fit containing counts, fit
#' configuration and results of the fit. Herein, functions for printing,
#' summarization and accessing are provided. The
#' functions \code{\link{enrichR}}, \code{\link{diffR}} and
#' \code{\link{regimeR}} generate a container of this class to save results of
#' a normR binomial mixture fitting. Please refer to their documentation for
#' conventional usage of the normR package.
#'
#' When working with instances of this S4 class, it is recommended to only use
#' functions to access contents of this object. Internally, the class holds a
#' map structure of unique elements to reduce memory requirements. #'
#'
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#'
#' @slot type A \code{character} representing the type of fit. One of
#' \code{c("enrichR","diffR", "regimeR")}.
#' @slot n An \code{integer} specifying the number of regions.
#' @slot ranges A \code{GenomicRanges} specifying the genomic coordinates of
#' the regions.
#' @slot k An \code{integer} giving the number of binomial mixture components.
#' @slot B An \code{integer} specifying the index of the background component.
#' @slot map A vector of \code{integer} holding a map to map back
#' \code{counts}, \code{lnposteriors}, \code{lnenrichment}, \code{lnpvals},
#' \code{lnqvals} and \code{classes}. See low level function
#' \code{normr:::map2uniquePairs} for how the map is generated.
#' @slot counts A \code{list} of length two containing a vector of
#' \code{integer} holding unique counts for control and treatment each. Use
#' \code{\link{getCounts}} to retrieve original count matrix.
#' @slot amount A vector of \code{integer} specifying the number of occurences
#' of each unique control / treatment count pair.
#' @slot names A \code{character} of length two specifying the names for
#' control and treatment.
#' @slot thetastar A \code{numeric} giving the calculated naive background
#' estimation, \emph{i.e.} \code{sum(getCounts(obj)[2,])/sum(getCounts(obj))}
#' @slot theta A \code{numeric} of length \code{k} giving the normR fitted
#' parametrization of \code{k} binomial mixture components.
#' @slot mixtures A \code{numeric} of length \code{k} giving the normR fitted
#' mixture proportions of \code{k} binomial mixture components. Should add up
#' to one.
#' @slot lnL A vector of \code{numeric} holding the log-likelihood-trace of
#' a normR model fit.
#' @slot eps A \code{numeric} used as threshold for normR fit EM convergence.
#' @slot lnposteriors A \code{matrix} with \code{length(amount)} rows and
#' \code{k} columns. It contains the ln posterior probabilities for each unique
#' control / treatment count pair. Use \code{\link{getPosteriors}} to get the
#' posterior matrix for the original data.
#' @slot lnenrichment A \code{numeric} of length \code{length(amount)} holding
#' calculared normalized enrichment for each unique control / treatment count
#' pair. The enrichment is calculated with respect to the fitted component
#' \code{B}. Use \code{\link{getEnrichment}} to retrieve enrichment for the
#' original data.
#' @slot lnpvals A \code{numeric} of length \code{length(amount)} holding ln
#' P-values for each unique control / treatment count pair. Given
#' \code{theta} of \code{B} the signifcane of enrichment is assigned. Use
#' \code{\link{getPvalues}} to retrieve Pvalues for original data.
#' @slot thresholdT An \code{integer} giving the threshold used to filter
#' P-values for FDR correction. The T-Filter threshold is a calculated
#' population size for which the null hypothesis (\code{theta} of \code{B}) can
#' be rejected. \code{eps} specifies the significance level.
#' @slot filteredT A vector of \code{integer} giving indices of P-values
#' passing \code{thresholdT}. Only these P-values will be considered for FDR
#' correction.
#' @slot lnqvals A \code{numeric} of length \code{length(filteredT)} holding
#' ln q-values (FDR correction). P-values are corrected for multiple testing
#' using Storey's method.
#' @slot classes A \code{integer} of length \code{length(amount)} specifying
#' the class assignments for each unique control / treatment count pair. These
#' class assignments are based on the normR model fit. For \code{type ==
#' "enrichR"}, this vector contains either \code{NA} (not enriched) or \code{1}
#' (enriched). For \code{type == "diffR"}, this vector contains \code{NA}
#' (unchanged), \code{1} (differential in ChIP-seq 1) and \code{2}
#' (differential in ChIP-seq 2). For \code{type == "regimeR"}, this vector
#' contains \code{NA} (not enriched) and an arbitary number of enrichment class
#' \code{>= 1}.
#'
#' @aliases NormRFit normRFit normrfit NormRfit
#'
#' @seealso \link{normr} for function creating this container
#'
#' @example inst/examples/NormRFit_example.R
#'
#' @import GenomicRanges
#' @import IRanges
#' @import grDevices
#' @import methods
#' @import parallel
#' @importFrom stats na.omit
#' @importFrom utils capture.output
#' @importFrom utils write.table
#'
#' @include methods.R
#' @export
setClass("NormRFit",
  representation = representation(type = "character",
                                  n = "integer",
                                  ranges = "GenomicRanges",
                                  k = "integer",
                                  B = "integer",
                                  map = "integer",
                                  counts = "list",
                                  amount = "integer",
                                  names = "character",
                                  thetastar = "numeric",
                                  theta = "numeric",
                                  mixtures = "numeric",
                                  lnL = "numeric",
                                  eps = "numeric",
                                  lnposteriors = "matrix",
                                  lnenrichment = "numeric",
                                  lnpvals = "numeric",
                                  thresholdT = "integer",
                                  filteredT = "integer",
                                  lnqvals = "numeric",
                                  classes = "integer")
)
setValidity("NormRFit",
  function(object) {
    if (!(object@type %in% c("enrichR", "diffR", "regimeR"))) {
      return("invalid type slot")
    }
    if (!is.null(object@ranges) &&
        object@n != length(object@ranges)) return("invalid n and ranges")
    if (length(object@k) == 0 | object@k <= 0) return("invalid k slot")
    if (object@B <= 0 || object@B > object@k) return("invalid B slot")
    if (max(object@map) != length(object@counts[[1]])) {
      return("incompatible map and count slot")
    }
    if (length(object@counts) != 2 ||
        length(object@counts[[1]]) != length(object@counts[[2]])) {
      return("invalid counts slot")
    }
    if (length(object@counts[[1]]) != length(object@amount) ||
        sum(object@amount) != object@n) {
      return("invalid amount slot")
    }
    if (length(object@thetastar) != 1 || object@thetastar < 0 ||
        object@thetastar > 1) {
      return("invalid thetastar slot")
    }
    if (length(object@theta) != object@k) return("invalid theta slot")
    if (length(object@mixtures) != object@k) return("invalid mixtures slot")
    n <- length(object@counts[[1]])
    if (NCOL(object@lnposteriors) != object@k ||
        NROW(object@lnposteriors) != length(object@counts[[1]])) {
      return("invalid lnposterios slot")
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
    if (object@thresholdT < 0) {
      return("invaled thresholdT slot")
    }
    if (length(object@lnqvals) != length(object@counts[[1]])) {
      return("invaled lnqvals slot")
    }
    if (length(object@classes) != length(object@counts[[1]])) {
      return("invaled classes slot")
    }
    TRUE
  }
)

#' @export
setGeneric("exportR", function(x, filename, ...)
           standardGeneric("exportR"))
#' @describeIn NormRFit Export results of a normR fit to common file formats.
#'
#' @param filename A \code{character} specifying the file to write to.
#' @param fdr \code{NA} or a \code{numeric} between \code{0} and \code{1}
#' specifying a FDR-level. Only regions with a q-value smaller than \code{fdr}
#' will be returned. If set to \code{NA}, all regions analyzed will be
#' returned (\code{\link{getRanges}}), classes are assigned by Maximum A
#' Posteriori (\code{\link{exportR}}).
#' @param color Specified color(s) when printing a bed file. If \code{x@type ==
#' "enrichR"}, \code{color} should of length 1 and color shading will be done
#' on this color. If \code{x@type == "diffR"}, \code{color} should be of length
#' 2 giving, firstly, the color for \code{control} and, secondly, the color for
#' \code{treatment}. If \code{x@type == "regimeR"}, \code{color} should be of
#' length \code{x@k-1}, specifying a color for each enrichment component. Per
#' default an appropriate color palette is used.
#' @param type A \code{character} specifying the filetype for exporting
#' results. If \code{NA}, format is guessed from \code{filename}'s extension.
#'
#' @aliases exportR
#'
#' @export
setMethod("exportR", signature=c("NormRFit", "character"),
  function(x, filename, fdr=.01, color=NA,
           type=c(NA, "bed", "bedGraph", "bigWig")) {
    typ <- match.arg(type)
    filename <- path.expand(filename)
    if (is.null(getRanges(x)))
      stop("no ranges set in x. Please set via 'x@ranges <- ranges'")
    if (!is.na(fdr))
      if (fdr < 0 | fdr > 1) stop("invalid fdr specified (0<=fdr<=1)")
    if (is.na(typ)) {
      ext <- tools::file_ext(filename)
      if (ext == "bed") typ <- "bed"
      else if (ext %in% c("bw", "bigWig")) typ <- "bigWig"
      else if (ext %in% c("bg", "bedGraph")) typ <- "bedGraph"
      else stop("type could not be inferred from filename extension. Specify!")
    }

    if (typ == "bed") {#qualitative output
      #clzzez contains integers (assigned regimes) and <NA> for background
      clzzez <- getClasses(x, fdr)
      nna <- which(!is.na(clzzez))
      clzzez <- clzzez[nna]#reduce classes to integer only
      if (is.na(fdr)) {
        score <- as.integer(1e3-getPosteriors(x)[nna,x@B]*1e3)
      } else {
        score <- as.integer(((fdr-getQvalues(x)[nna])/fdr)*1e3)
      }
      score[score < 0] <- 0
      score[score > 1e3] <- 1e3

      #retrieve coordinates
      gr <- getRanges(x)[nna]
      gr$component <- NULL
      gr$score <- score

      #name and color dependent on type of NormRFit
      getColRamp <- function(col, factor=5, steps=6) {
        colTab <- t(col2rgb(col))*factor
        colTab[colTab > 255] <- 255
        newCol <- rgb(colTab,maxColorValue=255)
        pal <- colorRampPalette(c(newCol, col))(steps)
        return(apply(col2rgb(pal), 2, paste, collapse=","))
      }
      if (x@type == "enrichR") {
        if (is.na(color[1])) color <- "gray50"
        if (length(color) != 1) {
          stop("invalid color argument for type 'enrichR' (length!=1)")
        }
        gr$name <- paste0("enrichR_score:", gr$score)
        gr$col <- getColRamp(color)[as.integer(gr$score/250)+1]

      } else if (x@type == "diffR") {
        if (is.na(color[1])) color <- c("darkblue", "darkred")
        if (length(color) != 2) {
          stop("invalid color argument for type 'diffR' (length!=2)")
        }
        clzzName <- c("ctrl", "treat")[clzzez]
        gr$name <- paste0("diffR_",clzzName,"_score:", gr$score)
        col <- rep(NA, length(gr))
        #Cond1
        col[which(clzzez==1)] <-
          getColRamp(color[1])[as.integer(gr$score[which(clzzez==1)]/250)+1]
        #Cond1
        col[which(clzzez==2)] <-
          getColRamp(color[2])[as.integer(gr$score[which(clzzez==2)]/250)+1]
        gr$col <- col

      } else if (x@type == "regimeR") {
        if (is.na(color[1])) color <- rev(terrain.colors(x@k+1))[-1]
        if (length(color) != x@k) {
          stop(paste0("invalid color argument for type 'regimeR' (length!=",
            x@k, ")"))
        }
        gr$name <- paste0("regimeR_",clzzez, "_score:", gr$score)
        #Loop through regimes and colors for creating color vector
        col <- rep(NA, length(gr))
        for (i in 1:x@k) {
          col[which(clzzez==i)] <-
            getColRamp(color[i])[as.integer(gr$score[which(clzzez==i)]/250)+1]
        }
        gr$col <- col
      }

      #prepare for output
      gr <- sort(gr)
      out <- as.data.frame(gr)[,c(1,2,3,7,6,5,2,3,8)]
      out[,c(2,7)] <- out[,c(2,7)]-1#0-based bed coordinates [start,end)
      out[,6] <- "."

      #writing
      cat(paste0('track name=', basename(filename), ' description="',
        filename, '" visibility=dense itemRgb="On"\n'), file=filename)
      utils::write.table(file=filename, x=out, sep="\t", col.names=FALSE,
        row.names=FALSE, quote=FALSE, append=TRUE)

    } else { #quantitative output
      gr <- getRanges(x)
      #e <- getEnrichment(x)
      #gr$score <- as.numeric(format(e,1,1))
      #idx <- which(getCounts(x)$control > 0 & getCounts(x)$treatment > 0)
      #gr <- gr[idx]
      gr$score <- getEnrichment(x)
      gr$score[which(abs(gr$score) < x@eps)] <- 0

      if (typ == "bedGraph") { #bedGraph - configure the trackline
        cat(paste0("track type=bedGraph name='", basename(filename),
                   "' description='", filename, "' visibility=full ",
                   "color=128,128,128 altColor=128,128,128 autoScale=off ",
                   "alwaysZero=on graphType=bar viewLimits=0:1.5 ",
                   "windowingFunction=mean\n"),
            file=filename)
        rtracklayer::export(object=gr, con=filename, format=typ, append=TRUE)

      } else { #bigWig - no trackline
        rtracklayer::export(object=gr, con=filename, format=typ)
      }
    }
})

#' @describeIn NormRFit Plot a NormRFit.
#'
#' @param y not used.
#' @aliases plot.NormRFit
#'
#' @export
setMethod("plot", c("NormRFit", "missing"),
  function(x, y, ...) {
    stop("not implemented yet")
   }
)


#' @export
setGeneric("getCounts", function(x) standardGeneric("getCounts"))
#' @describeIn NormRFit Get count data for control and treatment.
#'
#' @param x A \code{NormRFit} object.
#'
#' @return getCounts: A \code{list} of length 2 with \code{integer} for control
#' and treatment each.
#'
#' @aliases getCounts
#'
#' @export
setMethod("getCounts", "NormRFit", function(x) {
  return(list("control"=x@counts[[1]][x@map],"treatment"=x@counts[[2]][x@map]))
})

#' @export
setGeneric("getRanges", function(x, ...) standardGeneric("getRanges"))
#' @describeIn NormRFit Get the genomic coordinates of regions analyzed with
#' information about component assignment.
#'
#' @param k \code{NULL} or a \code{integer} specifying a model component for
#' which regions have to be returned. If set to \code{NULL}, regions are not
#' filtered on component assignments. If \code{fdr} is set and \code{k == x@B},
#' the function stops.
#'
#' @return getRanges: A \code{GenomicRanges} object.
#'
#' @aliases getRanges
#'
#' @export
setMethod("getRanges", "NormRFit", function(x, fdr=NA, k=NULL ) {
   if (!is.na(fdr)) if(fdr < 0 | fdr > 1) stop("invalid fdr specified")
   if (!is.null(k)) {
     if(k < 0 | k > x@k) stop("invalid k specified")
     if (k == 0 & !is.na(fdr)) stop("impossible to filter k with fdr")
   }

   gr <- x@ranges
   if (is.na(fdr)) {
     gr$component <- getClasses(x) #MAP class assignments
     if (!is.null(k)) {
       gr <- gr[which(gr$component == k)]
     }
     return(gr)
   } else if (!is.na(fdr)) {
     gr$component <- getClasses(x, fdr) #FDR-based class assignments
     gr <- gr[!is.na(gr$component)]
     if (!is.null(k)) {
       gr <- gr[which(gr$component == k)]
     }
     return(gr)
   }
})

#' @export
setGeneric("getPosteriors", function(x) standardGeneric("getPosteriors"))
#' @describeIn NormRFit Get computed posteriors for each mixture component.
#'
#' @return getPosteriors: A \code{matrix} of posteriors for \code{x@k} mixture
#' components
#'
#' @aliases getPosteriors
#'
#' @export
setMethod("getPosteriors", "NormRFit", function(x) {
  return(exp(x@lnposteriors)[x@map,])
})

#' @export
setGeneric("getEnrichment", function(x, ...) standardGeneric("getEnrichment"))
#' @describeIn NormRFit Get normalized enrichment.
#'
#' @param B An \code{integer} specifying the index of a mixture component. The
#' enrichment is calculated relative to this component used as a background
#' component. If \code{<NA>} (default), the background is determined by normR.
#' @param F An \code{integer} specifying the index of a mixture component. The
#' enrichment is calculated for this component over background \code{B}. If
#' \code{<NA>} (default), the component with theta closest to B is used
#' (\code{enrichR}, \code{regimeR}). For \code{diffR}, \code{F} is not
#' effective.
#' @param standardized A \code{logical} indicating if the enrichment should be
#' standardized betwen 0 and 1. A non-standardized enrichment is particular
#' useful when comparing intensities for ChIP-seq against the same antigen in
#' different conditions (default = TRUE).
#' @param procs An \code{integer} specifying the number of threads to use.
#'
#' @return getEnrichment: A \code{numeric} of length \code{length(x@n)} giving
#' the normR computed enrichment.
#'
#' @aliases getEnrichment
#'
#' @export
setMethod("getEnrichment", "NormRFit", function(x, B=NA, F=NA,
                                                standardized=TRUE,  procs=1L) {
  if (is.na(B) & is.na(F) & standardized) {
    return(x@lnenrichment[x@map])
  } else { #recomputation needed for specified B or non-standardized enrichment
    if (is.na(B)) B <- x@B
    if (B < 1 | B > x@k) stop("invalid B specified")
    if (is.na(F)) F <- x@k
    if (F == B | F < 1 | F > x@k) stop("invalid F specified")
    e <- normr:::computeEnrichmentWithMap(x@lnposteriors,
       list("values"=do.call(rbind, x@counts), "amount"=x@amount), x@theta,
      (F-1), (B-1), (x@type == "diffR"), standardized, procs)
    return(e[x@map])
  }
})

#' @export
setGeneric("getPvalues", function(x, ...) standardGeneric("getPvalues"))
#' @describeIn NormRFit Get normR-computed P-values.
#'
#' @param filtered A \code{logical} specifying if T-filtered or all
#' (default) P-values should be returned.
#'
#' @return getPvalues: A \code{numeric} of length \code{length(x@n)} giving
#' the normR computed Pvalues.
#'
#' @aliases getPvalues
#'
#' @export
setMethod("getPvalues", "NormRFit", function(x, filtered = FALSE) {
  if (filtered) {
    idx = x@map[which(x@map %in% x@filteredT)]
    exp(x@lnpvals)[idx]
  } else {
    exp(x@lnpvals)[x@map]
  }
})

#' @export
setGeneric("getQvalues", function(x) standardGeneric("getQvalues"))
#' @describeIn NormRFit Get FDR-corrected q-values.
#'
#' @return getQvalues: A \code{numeric} of length \code{length(x@filteredT)}
#' giving the FDR-corrected q-values using Storey's method.
#'
#' @aliases getQvalues
#'
#' @export
setMethod("getQvalues", "NormRFit", function(x) {
  exp(x@lnqvals)[x@map]
})

#' @export
setGeneric("getClasses", function(x, ...) standardGeneric("getClasses"))
#' @describeIn NormRFit Get component assignments for each region analyzed.
#'
#' @return getClasses: A \code{integer} specifying assignments of regions to
#' the mixture model. If \code{x@type == "enrichR"}, it contains \code{1} for
#' enriched regions and \code{NA} for non-enriched regions. If \code{x@type ==
#' "diffR"}, it contains \code{1} for control-enriched regions, \code{2} for
#' treatment-enriched regions and \code{NA} for non-enriched regions. If
#' \code{x@type == "regimeR"}, it contains \code{>= 1} for regime-enriched
#' regions and \code{NA} for non-enriched regions.
#'
#' @aliases getClasses
#'
#' @export
setMethod("getClasses", c("NormRFit"), function(x, fdr=NA) {
  if (is.na(fdr)) {
    x@classes[x@map]
  } else {
    #Filter on qvalues
    signf <- which(x@lnqvals <= log(fdr))
    clzzezSignf <- x@classes[signf]
    #assign classes based on Maximum A Posteriori for <NA>
    na <- which(is.na(clzzezSignf))
    if (length(na) > 0) {
      if (x@k == 2) {
        clzzezSignf[na] <- 1
      } else if (length(na) == 1) {
        clzzezSignf[na] <- which.max(x@lnposteriors[signf,][na,-x@B])
      } else {
        clzzezSignf[na] <-
          apply(x@lnposteriors[signf,][na,-x@B], 1, which.max)
      }
    }
    #create output vector with reverse mapping
    out <- rep(NA_integer_, length(x@classes))
    out[signf] <- clzzezSignf
    out[x@map]
  }
})

#' @describeIn NormRFit Returns the number of regions analyzed.
#'
#' @export
setMethod("length", "NormRFit", function(x) x@n)

#' @describeIn NormRFit Prints a small summary on a NormRFit.
#'
#' @export
setMethod("print", "NormRFit",
  function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("NormRFit-class object\n\n",
        "Type:                   ", x@type, "\n",
        "Number of Regions:      ", x@n, "\n",
        "Theta* (naive bg):      ", format(x@thetastar,digits=digits), "\n",
        "Background component B: ", x@B, "\n\n")
    if (length(x@theta)) {
      cat("+++ Results of fit +++ \nMixture Proportions:\n")
      mixtures <- paste0(format(x@mixtures*100, digits=digits), "%")
      names(mixtures)[x@B] <- "Background"
      names(mixtures)[-x@B] <- paste("Class", 1:(x@k-1))
      mixtures.string <- utils::capture.output(print.default(
        mixtures,print.gap=2L,quote=FALSE))
      cat(paste(mixtures.string, collapse="\n"), "\n")

      cat("Theta:\n")
      theta <- format(x@theta, digits=digits)
      names(theta)[x@B] <- "Background"
      names(theta)[-x@B] <- paste("Class", 1:(x@k-1))
      theta.string <- utils::capture.output(print.default(
        theta,print.gap=2L,quote=FALSE))
      cat(paste(theta.string, collapse="\n"), "\n")
    } else {
      cat("No results of fit.")
    }
    cat("\n")
    invisible(x)
  }
)

#' @describeIn NormRFit Shows a small summary on a NormRFit.
#'
#' @param object A \code{NormRCountConfig} object.
#'
#' @export
setMethod("show", "NormRFit", function(object) print(object))

#' @describeIn NormRFit Prints a concise summary of a NormRFit.
#'
#' @param print \code{logical()} indicating if summary should be print to
#' screen
#' @param digits Number of digits to show in number formatting.
#' @param ... optional arguments to be passed directly to the inherited
#' function without alteration and with the original names preserved.
#'
#' @export
setMethod("summary", "NormRFit",
   function(object, print=TRUE, digits=3, ...) {
     ans <- paste0("NormRFit-class object\n\n",
                  "Type:                  '", object@type, "'\n",
                  "Number of Regions:     ", object@n, "\n",
                  "Number of Components:  ", object@k, "\n",
                  "Theta* (naive bg):     ",
                  format(object@thetastar, digits=digits), "\n",
                  "Background component B: ", object@B, "\n\n")
    if (length(object@theta)) {
      ans <- paste0(ans, "+++ Results of fit +++ \n",
                    "Mixture Proportions:\n")
      mixtures <- paste0(format(object@mixtures*100, digits=digits), "%")
      names(mixtures)[object@B] <- "Background"
      names(mixtures)[-object@B] <- paste("Class", 1:(object@k-1))
      mixtures.string <- utils::capture.output(print.default(
        mixtures,print.gap=4L,quote=FALSE))
      ans <- paste0(ans, paste(mixtures.string, collapse="\n"), "\n")

      ans <- paste0(ans, "Theta:\n")
      theta <- format(object@theta, digits=digits)
      names(theta)[object@B] <- "Background"
      names(theta)[-object@B] <- paste("Class", 1:(object@k-1))
      theta.string <- utils::capture.output(print.default(
        theta,print.gap=4L,quote=FALSE))
      ans <- paste0(ans, paste(theta.string, collapse="\n"), "\n")

      ans <- paste0(ans, "\nBayesian Information Criterion:\t", format(
        (-2*object@lnL[length(object@lnL)]+length(object@theta)*log(object@n)),
        digits=digits), "\n\n")

      qvals <- getQvalues(object)
      nfiltered <- sum(is.na(qvals))
      ans <- paste0(ans, "+++ Results of binomial test +++ \n",
        "T-Filter threshold: ", object@thresholdT, "\n",
        "Number of Regions filtered out: ", nfiltered, "\n")

      ans <- paste0(ans,
        "Significantly different from background B based on q-values:\n",
        "TOTAL:\n")
      cts <- c(sum(qvals <= 0, na.rm=TRUE),
               sum(qvals <= 0.001, na.rm=TRUE),
               sum(qvals <= 0.01, na.rm=TRUE),
               sum(qvals <= 0.05, na.rm=TRUE),
               sum(qvals <= 0.1, na.rm=TRUE),
               sum(qvals > 0.1, na.rm=TRUE))
      cts <- matrix(c(format(cts - c(0,cts[1:4],0), digits=digits),
                    format(cts/sum(cts)*100,digits= digits)),
                    ncol=6, byrow=TRUE)
      colnames(cts) <- c("***", "**" , "*"  , "."  , " "  , "n.s." )
      rownames(cts) <- c("Bins", "%")
      cts.string <-
        utils::capture.output(
          print.default(cts,print.gap=4L,quote=FALSE, right=TRUE))
      ans <- paste0(ans, paste(cts.string, collapse="\n"), "\n")

      if (object@type %in% c("diffR", "regimeR")) { #print regime statistics
        clcts <- matrix(0, ncol=6, nrow=object@k-1,
          dimnames=list(1:(object@k-1), names(cts)))
        i = 1
        for (p in c(0,.001,.01,.05,.1)) {
          cl <- getClasses(object, p)
          tab <- table(na.omit(cl))
          for (n in names(tab)) {
            clcts[n,i] <- tab[n]
          }
          i = i + 1
        }
        clcts <- clcts - cbind(0, clcts[,1:4], 0)
        clcts[,6] <- sum(as.numeric(cts[1,])) - rowSums(clcts[,1:5])

        for (i in 1:(object@k-1)) {
          ans <- paste0(ans, "Class ", i, ":\n")
          out <- matrix(c(format(clcts[i,], digits=digits),
                          format(clcts[i,]/sum(clcts[i,])*100,digits=digits)),
                          ncol=6, byrow=TRUE)
          colnames(out) <- c("***", "**" , "*"  , "."  , " "  , "n.s." )
          rownames(out) <- c("Bins", "%")
          cts.string <- utils::capture.output(
            print.default(out, print.gap=4L,quote=FALSE, right=TRUE))
          ans <- paste0(ans, paste(cts.string, collapse="\n"), "\n")
        }
      }
      ans <- paste0(ans, "---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*'",
                    " 0.05 '.' 0.1 '  ' 1 'n.s.'\n\n")
    } else {#no fit found
      ans <- paste0(ans, "No results of fit.\n\n")
    }
    if (print) cat(ans)
    invisible(ans)
  }
)
