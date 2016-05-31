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
#' See also \code{map2uniquePairs()} for map generation.
#' @slot counts A \code{list()} of length 2 containing unique count pairs. Use
#' accessor \code{\link{getCounts}} to retrieve original count matrix.
#' @slot amount A \code{integer()} specifying the number of occurences of each
#' count pair
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
#' \code{getPosteriors()} to get the complete posterior matrix.
#' @slot lnenrichment A \code{numeric()}-vector of normalized ln enrichment over
#' component \code{B} for each unique count-tupel. Use accessor
#' \code{\link{getEnrichment}} to retrieve enrichment for original count sequence.
#' @slot lnpvals A \code{numeric()}-vector containing ln P-values every unique
#' 2-tupel in counts given \code{theta} of \code{B}.
#' @slot filteredT A \code{integer()}-vector giving indices of P-values to be
#' considered for FDR correction. These indices have been obtainted by filtering
#' P-values using the T method.
#' @slot thresholdT The threshold for filtering regions based on the sum of
#' counts in treatment and control.
#' @slot lnqvals A \code{numeric()}-vector containing ln q-values. These are
#' P-values corrected for multiple testing using Storey's method.
#' @slot classes A \code{integer()}-vector giving class assignments based on
#' model fit. For \code{obj@type == "enrichR"}, this vector contains either
#' \code{NA} (not enriched) or \code{1} (enriched). For \code{obj@type ==
#' "diffR"}, this vector contains \code{NA} (unchanged), \code{1} (differential
#' in ChIP-seq 1) and \code{2} (differential in ChIP-seq 2). For
#' \code{obj@type == "regimeR"}, this vector contains \code{NA} (not enriched)
#' and an arbitary number of enrichment class \code{>= 1}.
#' @aliases NormRFit
#' @seealso \code{\link{normr-methods}} for the functions that produce
#' this object
#' @return return values are described in the methods section.
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
                                  filteredT = "integer",
                                  thresholdT = "integer",
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

###
#DIAGNOSTIC FUNCTIONS FOR S4 DEFINED
###
#' @export
setMethod("print", "NormRFit",
  function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("NormRFit-class object\n\n",
        "Type:                  ", x@type, "\n",
        "Number of Regions:     ", x@n, "\n",
        "Theta* (naive bg):     ", format(x@thetastar,digits=digits), "\n",
        "Backgroundcomponent B: ", x@B, "\n\n")
    if (length(x@theta)) {
      cat("+++ Results of fit +++ \nMixture Proporitons:\n")
      cat(paste(format(x@mixtures,digits=digits), collapse="  "))
      cat("\nTheta:\n")
      cat(paste(format(x@theta,digits=digits), collapse="  "))
    } else {
      cat("No results of fit.")
    }
    cat("\n")
    invisible(x)
  }
)

#' @export
setMethod("show", "NormRFit", function(object) print(object))

#' @describeIn NormRFit Number of regions analyzed.
#' @aliases length
#' @export
setMethod("length", "NormRFit", function(x) x@n)

#' @describeIn NormRFit Prints a summary of the NormRFit object.
#' @param object A NormRFit object
#' @param print \code{logical()} indicating if summary should be print to screen
#' @param digits Number of digits to show in number formatting.
#' @param ... Not used
#' @return NULL
#' @export
setMethod("summary", "NormRFit",
   function(object, print=T, digits=3, ...) {
     ans <- paste0("NormRFit-class object\n\n",
                  "Type:                  '", object@type, "'\n",
                  "Number of Regions:     ", object@n, "\n",
                  "Number of Components:  ", object@k, "\n",
                  "Theta* (naive bg):     ",
                  format(object@thetastar, digits=digits), "\n",
                  "Backgroundcomponent B: ", object@B, "\n\n")
    if (length(object@theta)) {
       ans <- paste0(ans, "+++ Results of fit +++ \nMixture Proporitons:\n")
      ans <- paste0(ans,
        paste(format(object@mixtures,digits=digits), collapse="  "))
      ans <- paste0(ans, "\nTheta:\n")
      ans <- paste0(ans,
        paste(format(object@theta,digits=digits), collapse="  "))
      ans <- paste0(ans, "\nBayesian Information Criterion:\t", format(
        (-2*object@lnL[length(object@lnL)]+length(object@theta)*log(object@n)),
        digits=digits), "\n\n")
      
      qvals <- getQvalues(object)
      ans <- paste0(ans, "+++ Results of binomial test +++ \n", 
        "T-Filter threshold: ", object@thresholdT, "\n", 
        "Number of Regions filtered out: ", sum(is.na(qvals)), "\n")
      ans <- paste0(ans, 
        "Significantly different from background B based on q-values:\n")
      cts <- c(sum(qvals <= 0, na.rm=T),
               sum(qvals <= 0.001, na.rm=T),
               sum(qvals <= 0.01, na.rm=T),
               sum(qvals <= 0.05, na.rm=T),
               sum(qvals <= 0.1, na.rm=T),
               sum(qvals > 0.1, na.rm=T))
      names(cts) <- c("***", "**" , "*"  , "."  , " "  , "n.s." )
      cts <- cts - c(0,cts[1:4],0)
      cts.string <-
         capture.output(print.default(format(cts, digits=digits), print.gap=2L,
                                     quote=F))
      ans <- paste0(ans, paste(cts.string, collapse="\n"), "\n")
      ans <- paste0(ans, "---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*'",
                          " 0.05 '.' 0.1 '  ' 1 'n.s.'\n\n")
#TODO
#      if (object@type != "diffR") { #print regime statistics
#        for (p in c(0,.001,.01,.05,.1)) {
#          cl = getClasses(object, p)
#
#        }
#      }
#      if (object@type != "regimeR") { #print regime statistics
#        for (p in c(0,.001,.01,.05,.1)) {
#          cl = getClasses(object, p)
#
#        }
#      }
    } else {
      ans <- paste0(ans, "No results of fit.\n\n")
    }
    if (print) cat(ans)
    invisible(ans)
  }
)

###
#ACCESSOR METHODS
###
#' @export
setGeneric("getCounts", function(obj) standardGeneric("getCounts"))
#' @describeIn NormRFit Retrieve read count data.
#' @aliases getCounts
#' @export
setMethod("getCounts", "NormRFit", function(obj) {
  return(list("control"=obj@counts[[1]][obj@map],"treatment"=obj@counts[[2]][obj@map]))
})

#' @export
setGeneric("getPosteriors", function(obj) standardGeneric("getPosteriors"))
#' @describeIn NormRFit Retrieve computed posteriors.
#' @aliases getPosteriors
#' @export
setMethod("getPosteriors", "NormRFit", function(obj) {
  return(exp(obj@lnposteriors)[obj@map,])
})

#' @export
setGeneric("getEnrichment", function(obj, ...) standardGeneric("getEnrichment"))
#' @describeIn NormRFit Retrieve NormR calculated enrichment.
#' @param B Index of component to compute the enrichment to. If \code{<NA>}, the
#' background is determined by normR.(DEFAULT NA)
#' @param procs Number of threads.
#' @aliases getEnrichment
#' @export
setMethod("getEnrichment", "NormRFit", function(obj, B=NA, procs=1) {
  if (is.na(B)) {
    return(obj@lnenrichment[obj@map])
  } else {
    if (B < 1 | B > obj@k) stop("invalid B specified")
    e <- normr::computeEnrichmentWithMap(obj@lnposteriors, 
      list("values"=do.call(rbind, obj@counts), "amount"=obj@amount), obj@theta,
      (obj@k-1), (B-1), (obj@type == "diffR"), 1L)
    return(e[obj@map])
  }
})

#' @export
setGeneric("getPvalues", function(obj, ...) standardGeneric("getPvalues"))
#' @describeIn NormRFit Retrieve computed P-values.
#' @aliases getPvalues
#' @export
setMethod("getPvalues", "NormRFit", function(obj, filtered=F) {
  if (filtered) {
    idx = obj@map[which(obj@map %in% obj@filteredT)]
    exp(obj@lnpvals)[idx]
  } else {
    exp(obj@lnpvals)[obj@map]
  }
})

#' @export
setGeneric("getQvalues", function(obj) standardGeneric("getQvalues"))
#' @describeIn NormRFit Retrieve computed Q-values. See Bioconductor package
#' \link{qvalue}.
#' @aliases getQvalues
#' @export
setMethod("getQvalues", "NormRFit", function(obj) {
  exp(obj@lnqvals)[obj@map]
})

#' @export
setGeneric("getClasses", function(obj, ...) standardGeneric("getClasses"))
#' @describeIn NormRFit Retrieve classes for each bin, i.e. enriched vs
#' non-enriched (enrichR), differential enrichment (diffR) and enrichment
#' regimes (regimeR).
#' @aliases getClasses
#' @export
setMethod("getClasses", c("NormRFit"), function(obj, fdr=NA) {
  if (is.na(fdr)) {
    obj@classes[obj@map]
  } else {
    #Filter on qvalues
    signf <- which(obj@lnqvals <= log(fdr))
    clzzezSignf <- obj@classes[signf]
    #assign classes based on Maximum A Posteriori for <NA>
    na <- which(is.na(clzzezSignf))
    if (length(na) > 0) {
      if (obj@k == 2) {
        clzzezSignf[na] <- 1
      } else if (length(na) == 1) {
        clzzezSignf[na] <- which.max(obj@lnposteriors[signf,][na,-obj@B])
      } else {
        clzzezSignf[na] <-
          apply(obj@lnposteriors[signf,][na,-obj@B], 1, which.max)
      }
    }
    #create output vector with reverse mapping
    out <- rep(NA_integer_, length(obj@classes))
    out[signf] <- clzzezSignf
    out[obj@map]
  }
})

#' @export
setGeneric("getRanges", function(obj, ...) standardGeneric("getRanges"))
#' @describeIn NormRFit Retrieve read count data.
#' @aliases getRanges
#' @export
setMethod("getRanges", "NormRFit", function(obj, fdr=NA, k=NULL ) {
   if (!is.na(fdr)) if(fdr < 0 | fdr > 1) stop("invalid fdr specified")
   if (!is.null(k)) {
     if(k < 0) stop("invalid k specified")
     if (k == obj@B & !is.na(fdr)) stop("impossible to filter k with fdr")
   }

   gr <- obj@ranges
   if (is.na(fdr)) {
     gr$component <- getClasses(obj) #MAP class assignments
     if (!is.null(k)) {
       gr <- gr[gr$component == k]
     }
     return(gr)
   } else if (!is.na(fdr)) {
     gr$component <- getClasses(obj, fdr) #FDR-based class assignments
     gr <- gr[!is.na(gr$component)]
     if (!is.null(k)) {
       gr <- gr[which(gr$component == k)]
     }
     return(gr)
   }
})

###
#EXPORTING
###
##' @export
#setGeneric("recomputeR", function(obj, filename, type, ...)
#           standardGeneric("recomputeR"))

###
#PLOTTING
###
setMethod("plot", "NormRFit",
  function(x, ...) {
    stop("not implemented yet")
  }
)

###
#EXPORTING
###
#' @export
setGeneric("exportR", function(obj, filename, ...)
           standardGeneric("exportR"))
#' @describeIn NormRFit Export results of a \link{NormRFit} to either bed,
#' bedGraph or bigWig.
#' @param filename File to write results to.
#' @param fdr The q-value threshold for calling a region significant. If
#' \code{NA}, no filtering is done on q-values but regions are reported by
#' Maximum A Posteriori.
#' @param color Specified color(s) for printing a ped file. If obj is of type
#' "enrichR", this argument should of length 1 #' and bed color shading will be
#' done on this color. If \code{obj} is of type
#' "diffR", this argument should be of length2
#' @param type Type of file used for exporting results. If <NA>, format is
#' guessed from \code{filename}'s extension.
#' @aliases exportR
#' @export
setMethod("exportR", signature=c("NormRFit", "character"),
  function(obj, filename, fdr=.01, color=NA, 
           type=c(NA, "bed", "bedGraph", "bigWig")) {
    typ <- match.arg(type)
    filename <- path.expand(filename)
    if (is.null(getRanges(obj)))
      stop("no ranges set in obj. Please set via 'obj@ranges <- ranges'")
    if (!is.na(fdr) || fdr < 0 | fdr > 1) stop("invalid fdr specified (0<=fdr<=1)")
    if (is.na(typ)) {
      ext <- tools::file_ext(filename)
      if (ext == "bed") typ <- "bed"
      else if (ext %in% c("bw", "bigWig")) typ <- "bigWig"
      else if (ext %in% c("bg", "bedGraph")) typ <- "bedGraph"
      else stop("type could not be inferred from filename extension. Specify!")
    }

    if (typ == "bed") {#qualitative output
      #clzzez contains integers (assigned regimes) and <NA> for background
      clzzez <- getClasses(obj, fdr)
      nna <- which(!is.na(clzzez))
      clzzez <- clzzez[nna]#reduce classes to integer only
      if (is.na(fdr)) {
        score <- as.integer(1e3-getPosteriors(obj)[nna,obj@B]*1e3)
      } else {
        score <- as.integer(((fdr-getQvalues(obj)[nna])/fdr)*1e3)
      }
      score[score < 0] <- 0
      score[score > 1e3] <- 1e3

      #retrieve coordinates
      gr <- getRanges(obj)[nna]
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
      if (obj@type == "enrichR") {
        if (is.na(color)) color <- "gray50"
        if (length(color) != 1) {
          stop("invalid color argument for type 'enrichR' (length!=1)")
        }
        gr$name <- paste0("enrichR_score:", gr$score)
        gr$col <- getColRamp(color)[as.integer(gr$score/250)+1]

      } else if (obj@type == "diffR") {
        if (is.na(color)) color <- c("darkred", "darkblue")
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

      } else if (obj@type == "regimeR") {
        if (is.na(color)) color <- rev(terrain.colors(obj@k+1))[-1]
        if (length(color) != obj@k) {
          stop(paste0("invalid color argument for type 'regimeR' (length!=",
            obj@k, ")"))
        }
        gr$name <- paste0("regimeR_",clzzez, "_score:", gr$score)
        #Loop through regimes and colors for creating color vector
        col <- rep(NA, length(gr))
        for (i in 1:obj@k) {
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
      write.table(file=filename, x=out, sep="\t", col.names=F, row.names=F,
        quote=F, append=T)

    } else { #quantitative output
      gr <- getRanges(obj)
      #e <- getEnrichment(obj)
      #gr$score <- as.numeric(format(e,1,1))
      #idx <- which(getCounts(obj)$control > 0 & getCounts(obj)$treatment > 0)
      #gr <- gr[idx]
      gr$score <- getEnrichment(obj)
      gr$score[which(abs(gr$score) < obj@eps)] <- 0

      if (typ == "bedGraph") { #bedGraph - configure the trackline
        cat(paste0("track type=bedGraph name='", basename(filename),
                   "' description='", filename, "' visibility=full ",
                   "color=128,128,128 altColor=128,128,128 autoScale=off alwaysZero=on ",
                   "graphType=bar viewLimits=0:1.5 windowingFunction=mean\n"),
            file=filename)
        rtracklayer::export(object=gr, con=filename, format=typ, append=T)

      } else { #bigWig - no trackline
        rtracklayer::export(object=gr, con=filename, format=typ)
      }
    }
})
