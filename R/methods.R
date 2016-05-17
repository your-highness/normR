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
#' Normalization and difference calling in ChIP-seq data
#'
#' Robust normalization and difference calling procedures for ChIP-seq and
#' alike data. Count vectors can be obtained from experiment bam files in two
#' experiments. Counts are jointly modeled as a multinomial sampling trial.
#' A binomial mixture model with a given number of components is fit and used
#' for identifying enriched or depleted regions in two given data tracks.
#' Log-space multinomial model is fit by Expectation maximization in C++.
#'
#' @name normr
#' @aliases normR
#' @aliases diffR
#' @aliases enrichR
#' @aliases PeakCalling
#' @aliases DifferentialPeakCalling
#' @aliases EnrichmentCalling
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import IRanges
#' @import Rcpp
#' @import qvalue
#' @import parallel
#' @import bamsignals
#' @docType package
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#' @seealso \code{\link{normr-methods}} for available functions
#' @include NormRFit.R
#' @include BamCountConfig.R
#' @useDynLib normr, .registration=TRUE
NULL

#' Enrichment, Difference and Regime Calling in ChIP-seq data.
#'
#' A correct background estimation is crucial for calling enrichment in ChIP-seq
#' data.  These functions implement the mixture modeling for two given ChIP-seq
#' tracks by joint multinomial modeling. A predefined number of model components
#' is fit simultaneously to binned read count data via an efficient Expectation
#' Maximization implementation in C++. Convergence is achieved by a threshold on
#' the minimum change in model loglikelihood. After the model is fit, every bin
#' is tested for significance against the fitted background component. Moreover,
#' a standardized enrichment for each bin is calculated based on the fitted
#' background component.
#'
#' @param treatment A \link{integer} vector of treatment counts or a
#' \link{character} pointing to the treatment bam file. In the latter case an
#' \code{treatment}.bai index file should exist in the same folder. Should be
#' consistent with control.
#' @param control A \link{integer} vector of control counts or a
#' \link{character} pointing to the control bam file. In the latter case an
#' \code{control}.bai index file should exist in the same folder. Should be
#' consistent with treatment.
#' @param genome Either \code{NULL}, an USCS genome identifier for
#' \code{GenomeInfoDb::fetchExtendedChromInfoFromUCSC} (only assembled
#' molecules, circular omitted) or a \link{data.frame} consisting of two
#' columns(1st column=chromosome names, 2nd column=lengths). Chromosome names
#' given should be present in the bam file header. If \code{NULL}, chromosome
#' names will be read from treatment bamheader. Please be aware that bamheader
#' might contain irregular contigs and chrM which influence the diffR fit.
#' (DEFAUlT=NULL).
#' @param models Number of model components for \code{regimeR()}. Should be an
#' \code{integer() >= 3}. (DEFAULT=3)
#' @param eps Threshold for EM convergence.
#' @param procs Number of threads to use
#' @param countConfig A BamCountConfig object specifying bam counting parameters
#' for read count retrieval. See \link{BamCountConfig}.
#' @param verbose A logical value indicating whether verbose output is desired
#' @return a \link{NormRFit} object
#'
#' @seealso BamCountConfig-class
#' @seealso NormRFit-class
#' @name normr-methods
#' @example inst/examples/methods_example.R
NULL

#HELPER FUNCTIONS
handleCharCharChar <- function(treatment, control, genome, verbose) {
  treatment <- path.expand(treatment); control <- path.expand(control)
  if (!file.exists(treatment)) stop("treatment is not a file")
  if (!file.exists(control)) stop("control is not a file")
  if (treatment == control) stop("treatment and control are identical")

  if (verbose) {
    message(paste0("Getting genome coordinates for ", genome, " ..."))
  }

  #UseGenomeInfoDb to fetch information on regular, non-circular chroms
  genome <- fetchExtendedChromInfoFromUCSC(genome)
  idx <- which(!genome$circular & genome$SequenceRole=="assembled-molecule")
  genome <- genome[idx,1:2]

  return(genome)
}
handleCharCharDf <- function(treatment, control, genome, countConfig, procs,
                             verbose) {
  treatment <- path.expand(treatment); control <- path.expand(control)
  if(!file.exists(paste(treatment, ".bai", sep=""))) {
    stop("No index file for", treatment, ".\n")
  }
  if(!file.exists(paste(control, ".bai", sep=""))) {
    stop("No index file for", control, ".\n")
  }
  if (NCOL(genome) != 2) stop("invalid genome data.frame")

  if (verbose) {
    message(paste0("Counting on ", control, " & ", treatment,
                   " for specified genome coordinates..."))
  }

  gr <- GenomicRanges::GRanges(genome[,1], IRanges(1,genome[,2]))
  counts <- parallel::mcmapply(
    bamsignals::bamProfile, bampath=c(treatment, control),
    MoreArgs=list(gr=gr, binsize=countConfig@binsize,
                  mapq=countConfig@mapq,
                  shift=countConfig@shift,
                  paired.end=getFilter(countConfig),
                  tlenFilter=countConfig@tlenFilter,
                  filteredFlag=countConfig@filteredFlag,
                  verbose=F),
    mc.cores=procs, SIMPLIFY=F
  )
  counts[[1]] <- unlist(as.list(counts[[1]]))
  counts[[2]] <- unlist(as.list(counts[[2]]))

  #Give bins across the supplied genome
  gr <- unlist(tile(gr, width=countConfig@binsize))
  seqinfo(gr) <- Seqinfo(as.character(genome[,1]), genome[,2])

  return(list(counts=counts, gr=gr))
}
handleCharCharGR <- function(treatment, control, gr, countConfig, procs,
                             verbose) {
  treatment <- path.expand(treatment); control <- path.expand(control)
  if(!file.exists(paste(treatment, ".bai", sep=""))) {
    stop("No index file for", treatment, ".\n")
  }
  if(!file.exists(paste(control, ".bai", sep=""))) {
    stop("No index file for", control, ".\n")
  }
  if (NCOL(genome) != 2) stop("invalid genome data.frame")

  if (verbose) {
    message(paste0("Counting on ", control, " & ", treatment,
                   " for specified GenomicRanges..."))
  }

  counts <- parallel::mcmapply(
    bamsignals::bamProfile, bampath=c(treatment, control),
    MoreArgs=list(gr=gr, binsize=countConfig@binsize,
                  mapq=countConfig@mapq,
                  shift=countConfig@shift,
                  paired.end=getFilter(countConfig),
                  tlenFilter=countConfig@tlenFilter,
                  filteredFlag=countConfig@filteredFlag,
                  verbose=F),
    mc.cores=procs, SIMPLIFY=F
  )
  counts[[1]] <- unlist(as.list(counts[[1]]))
  counts[[2]] <- unlist(as.list(counts[[2]]))

  return(list(counts=counts, gr=gr))
}

#' \code{enrichR}: Enrichment calling between \code{treatment} (ChIP-seq) and
#' \code{control} (Input) for either two given read count vectors and an
#' optional \link{GenomicRanges} object that defines the original regions
#' ('integer,integer,GenomicRangesOrNULL'), two given bam filepaths and a
#' genome \code{data.frame} defining chromosomes which will be binned according
#' to 'countConfig' ('character,character,data.frame') or two given bam
#' filepaths and a genome \code{character} specifying a UCSC genome identifier
#' (e.g. "hg19") that is used for retrieval of chromosome annotation which will
#' be binned according to 'countConfig' ('character,character,character').
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
#' @export
setGeneric("enrichR", function(treatment, control, genome, ...)
  standardGeneric("enrichR"))
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("integer", "integer", "GenomicRanges"),
  function(treatment, control, genome, eps=1e-5, iterations=10,
           procs=1L, verbose=TRUE) {
    if (length(treatment) != length(control)) {
      stop("incompatible treatment and control count vectors")
    }
    fit <- normr:::normr_core(control, treatment, 2L, eps, iterations, 0, F,
                              verbose, procs)

    #Storey's q-value on T filtered P-values
    if (verbose) message("... computing Q-values.")
    idx <- which(fit$map$map %in% fit$filteredT)
    lnqvals <- as.numeric(rep(NA,length(treatment)))
    lnqvals[idx] <-
      log(qvalue(exp(normr:::mapToOriginal(fit$lnpvals, fit$map)[idx]))$qvalues)
    lnqvals <- normr:::mapToUniqueWithMap(lnqvals, fit$map)

    #Create classes vector w/ Maximum A Posteriori (<NA> isset for background)
    #Note that Maximum A Posteriori and P-/Q-values are sometimes not in
    #correspondence!
    classes <- apply(fit$lnpost,1,which.max)
    classes[classes == 1] <- NA_integer_
    classes <- classes-1

    #NormRFit-class object
    o <- new("NormRFit", type="enrichR", n=length(treatment), ranges=genome,
             k=2L, B=1L, map=fit$map$map,
             counts=list(as.integer(fit$map$values[1,]),
                         as.integer(fit$map$values[2,])),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             lnqvals=lnqvals, classes=as.integer(classes))

    #Print logging information
    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=T)
    }
    return(o)
})
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("character", "character", "GenomicRanges"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharGR(treatment, control, genome, countConfig, procs,
                                 verbose)
    return(enrichR(countsGr$counts[[1]], countsGr$counts[[2]], genome, eps,
      iterations, procs, verbose))
})
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("character", "character", "data.frame"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig, procs,
                                 verbose)
    return(enrichR(countsGr$counts[[1]], countsGr$counts[[2]], countsGr$gr, eps,
      iterations, procs, verbose))
})
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("character", "character", "character"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome, verbose)
    return(enrichR(treatment, control, genome, countConfig, eps,
      iterations, procs, verbose))
})

#' \code{diffR}: Difference calling between \code{treatment} (ChIP-seq 1) and
#' \code{control} (ChIP-seq 2) for either two given read count vectors and an
#' optional \link{GenomicRanges} object that defines the original regions
#' ('integer,integer,GenomicRangesOrNULL'), two given bam filepaths and a
#' genome \code{data.frame} defining chromosomes which will be binned according
#' to 'countConfig' ('character,character,data.frame') or two given bam
#' filepaths and a genome \code{character} specifying a UCSC genome identifier
#' (e.g. "hg19") that is used for retrieval of chromosome annotation which will
#' be binned according to 'countConfig' ('character,character,character').
#' @aliases diffR
#' @aliases differenceCall
#' @rdname normr-methods
#' @export
setGeneric("diffR", function(treatment, control, genome, ...)
  standardGeneric("diffR"))
#' @rdname normr-methods
#' @export
setMethod("diffR", signature("integer", "integer", "GenomicRanges"),
  function(treatment, control, genome, eps=1e-5, iterations=10,
           procs=1L, verbose=TRUE) {
    if (length(treatment) != length(control)) {
      stop("incompatible treatment and control count vectors")
    }
    fit <- normr:::normr_core(control, treatment, 3L, eps, iterations, 1, T,
                              verbose, procs)

    #Pvalues from two sided test are marginally above 1
    fit$lnpvals[which(fit$lnpvals > 0)] <- 0

    #Storey's q-value on T filtered P-values
    if (verbose) message("... computing Q-values.")
    idx <- which(fit$map$map %in% fit$filteredT)
    lnqvals <- as.numeric(rep(NA,length(treatment)))
    lnqvals[idx] <-
      log(qvalue(exp(normr:::mapToOriginal(fit$lnpvals, fit$map)[idx]))$qvalues)
    lnqvals <- normr:::mapToUniqueWithMap(lnqvals, fit$map)

    #Create classes vector w/ Maximum A Posteriori (<NA> isset for background)
    #Note that Maximum A Posteriori and P-/Q-values are sometimes not in
    #correspondence!
    classes <- apply(fit$lnpost[,c(1,3,2)],1,which.max)
    classes[classes == 3] <- NA_integer_

    o <- new("NormRFit", type="diffR", n=length(treatment), ranges=genome,
             k=3L, B=2L, map=fit$map$map,
             counts=list(as.integer(fit$map$values[1,]),
                         as.integer(fit$map$values[2,])),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             lnqvals=lnqvals, classes=as.integer(classes))

    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=T)
    }
    return(o)
})
#' @rdname normr-methods
#' @export
setMethod("diffR", signature("character", "character", "GenomicRanges"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharGR(treatment, control, genome, countConfig, procs,
                                 verbose)
    return(diffR(countsGr$counts[[1]], countsGr$counts[[2]], genome, eps,
      iterations, procs, verbose))
})
#' @rdname normr-methods
#' @export
setMethod("diffR", signature("character", "character", "data.frame"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig, procs,
                                 verbose)
    return(diffR(countsGr$counts[[1]], countsGr$counts[[2]], countsGr$gr, eps,
      iterations, procs, verbose))
})
#' @rdname normr-methods
#' @export
setMethod("diffR", signature("character", "character", "character"),
  function(treatment, control, genome="", countConfig=countConfigSingleEnd(),
           eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome, verbose)
    return(diffR(treatment, control, genome, countConfig, eps, iterations,
      procs, verbose))
})

#' \code{regimeR}: Enrichment regime calling between \code{treatment}
#' (ChIP-seq 1) and \code{control} (ChIP-seq 2) with with \code{models} number
#' of models for either two given read count vectors and an optional
#' \link{GenomicRanges} object that defines the original regions
#' ('integer,integer,GenomicRangesOrNULL'), two given bam filepaths and a
#' genome \code{data.frame} defining chromosomes which will be binned according
#' to 'countConfig' ('character,character,data.frame') or two given bam
#' filepaths and a genome \code{character} specifying a UCSC genome identifier
#' (e.g. "hg19") that is used for retrieval of chromosome annotation which will
#' be binned according to 'countConfig' ('character,character,character').
#' @aliases regimeR
#' @aliases regimeCall
#' @rdname normr-methods
#' @export
setGeneric("regimeR", function(treatment, control, genome, models, ...)
  standardGeneric("regimeR"))
#' @rdname normr-methods
#' @export
setMethod("regimeR",
          signature("integer", "integer", "GenomicRanges", "integer"),
  function(treatment, control, genome, models=3L, eps=1e-5,
            iterations=10, procs=1L, verbose=TRUE) {
    if (models <= 2) stop("invalid number of models specified")
    if (length(treatment) != length(control)) {
      stop("incompatible treatment and control count vectors")
    }
    models = as.integer(models)
    fit <- normr:::normr_core(control, treatment, models, eps, iterations, 0, F,
                              verbose, procs)

    #Storey's q-value on T filtered P-values
    if (verbose) message("... computing Q-values.")
    idx <- which(fit$map$map %in% fit$filteredT)
    lnqvals <- as.numeric(rep(NA,length(treatment)))
    lnqvals[idx] <-
      log(qvalue(exp(normr:::mapToOriginal(fit$lnpvals, fit$map)[idx]))$qvalues)
    lnqvals <- normr:::mapToUniqueWithMap(lnqvals, fit$map)

    #Create classes vector w/ Maximum A Posteriori (<NA> isset for background)
    #Note that Maximum A Posteriori and P-/Q-values are sometimes not in
    #correspondence!
    classes <- apply(fit$lnpost,1,which.max)
    classes[classes == 1] <- NA_integer_
    classes <- classes-1

    #NormRFit-class object
    o <- new("NormRFit", type="regimeR", n=length(treatment), ranges=genome,
             k=models, B=1L, map=fit$map$map,
             counts=list(as.integer(fit$map$values[1,]),
                         as.integer(fit$map$values[2,])),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             lnqvals=lnqvals, classes=as.integer(classes))

    #Print logging information
    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=T)
    }
    return(o)
})
#' @rdname normr-methods
#' @export
setMethod("regimeR",
          signature("character", "character", "GenomicRanges", "integer"),
  function(treatment, control, genome, models=3L,
           countConfig=countConfigSingleEnd(), eps=1e-5,
           iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharGR(treatment, control, genome, countConfig, procs,
                                 verbose)
    return(regimeR(countsGr$counts[[1]], countsGr$counts[[2]], genome, models,
      eps, iterations, procs, verbose))
})
#' @rdname normr-methods
#' @export
setMethod("regimeR", signature("character", "character", "data.frame", "integer"),
  function(treatment, control, genome, models=3L,
           countConfig=countConfigSingleEnd(), eps=1e-5, iterations=10,
           procs=1L, verbose=TRUE) {
    if (models <= 2) stop("invalid number of models specified")
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig, procs,
                                 verbose)
    return(regimeR(countsGr$counts[[1]], countsGr$counts[[2]], countsGr$gr,
      models, eps, iterations, procs, verbose))
})
#' @rdname normr-methods
#' @export
setMethod("regimeR", signature("character", "character", "character", "integer"),
  function(treatment, control, genome="", models=3L,
           countConfig=countConfigSingleEnd(), eps=1e-5,
           iterations=10, procs=1L, verbose=TRUE) {
    if (models <= 2) stop("invalid number of models specified")
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome, verbose)
    return(regimeR(treatment, control, genome, models, countConfig, eps,
      iterations, procs, verbose))
})
