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
#' @import Rsamtools
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
#' @param fdr The FDR threshold for class definition. See \link{NormRFit}.
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

#new class uniting GenomicRanges and NULL
setClassUnion("GenomicRangesOrNULL", c("GenomicRanges", "NULL"))

#HELPER FUNCTIONS
handleCharCharChar <- function(treatment, control, genome) {
  treatment <- path.expand(treatment); control <- path.expand(control)
  if (!file.exists(treatment)) stop("treatment is not a file")
  if (!file.exists(control)) stop("control is not a file")
  if (treatment == control) stop("treatment and control are identical")

  #UseGenomeInfoDb to fetch information on regular, non-circular chroms
  genome <- fetchExtendedChromInfoFromUCSC(genome)
  idx <- which(!genome$circular & genome$SequenceRole=="assembled-molecule")
  genome <- genome[idx,1:2]

  return(genome)
}
handleCharCharDf <- function(treatment, control, genome, countConfig, procs) {
  treatment <- path.expand(treatment); control <- path.expand(control)
  if(!file.exists(paste(treatment, ".bai", sep=""))) {
    stop("No index file for", treatment, ".\n")
  }
  if(!file.exists(paste(control, ".bai", sep=""))) {
    stop("No index file for", control, ".\n")
  }
  if (NCOL(genome) != 2) stop("invalid genome data.frame")

  gr <- GRanges(genome[,1], IRanges(1,genome[,2]))
  counts <- mcmapply(bamsignals::bamProfile, bampath=c(treatment, control),
                     MoreArgs=list(gr=gr, binsize=countConfig@binsize,
                                   mapq=countConfig@mapq,
                                   shift=countConfig@shift,
                                   paired.end=getFilter(countConfig),
                                   #Tlen.filter not yet in Bioconductor
                                   #tlen.filter=countConfig@tlenFilter,
                                   #filteredFlag not yet in Bioconductor
                                   #filteredFlag=countConfig@requiredFlag,
                                   verbose=F),
                     mc.cores=procs, SIMPLIFY=F)

  #Give bins across the supplied genome
  gr <- unlist(tile(gr, width=countConfig@binsize))
  seqlengths(gr) <- genome[,2]

  return(list(counts=counts, gr=gr))
}

#' @export
setGeneric("enrichR", function(treatment, control, genome, ...)
  standardGeneric("enrichR"))
#' \code{enrichR}: Enrichment calling between \code{treatment} (ChIP-seq) and
#' \code{control} (Input) for two given read count vectors and an optional
#' \link{GRanges} object that defines the original regions.
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("integer", "integer", "GenomicRangesOrNULL"),
  function(treatment, control, genome=GRanges(), fdr=5e-2, eps=1e-5, iterations=10,
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

    #Get classes vector
    classes <- as.integer(rep(NA, length(lnqvals)))
    classes[which(lnqvals < log(fdr))] <- 1L

    #NormRFit-class object
    o <- new("NormRFit", type="enrichR", n=length(treatment), ranges=genome,
             k=2L, B=1L, map=fit$map$map,
             counts=list(as.integer(fit$map$values[1,]),
                         as.integer(fit$map$values[2,])),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             lnqvals=lnqvals, classes=classes)

    #Print logging information
    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=T)
    }
    return(o)
})
#' \code{enrichR}: Enrichment calling between \code{treatment} (ChIP-seq) and
#' \code{control} (Input) for two given bam filepaths. The given \code{genome}
#' defines chromosomes which will be binned according to \code{countConfig}.
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("character", "character", "data.frame"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           fdr=5e-2, eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig, procs)
    return(enrichR(countsGr$counts[[1]][1], countsGr$counts[[2]][1],
      countsGr$gr, fdr, eps, iterations, procs, verbose))
})
#' \code{enrichR}: Enrichment calling between \code{treatment} (ChIP-seq) and
#' \code{control} (Input) for two given bam filepaths. The given \code{genome}
#' defines a UCSC genome identifier (e.g. "hg19") that is used for retrieval
#' of chromosome annotation. Chromosomes will be binned according to
#' \code{countConfig}.
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
setMethod("enrichR", signature("character", "character", "character"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           fdr=5e-2, eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome)
    return(enrichR(treatment, control, genome, countConfig, fdr, eps,
      iterations, procs, verbose))
})

#' @export
setGeneric("diffR", function(treatment, control, genome, ...)
  standardGeneric("diffR"))
#' \code{diffR}: Difference calling between \code{treatment} (ChIP-seq 1) and
#' \code{control} (ChIP-seq 2) for two given read count vectors and an optional
#' \link{GRanges} object that defines the original regions.
#' @aliases diffR
#' @aliases differenceCall
#' @rdname normr-methods
#' @export
setMethod("diffR", signature("integer", "integer", "GenomicRangesOrNULL"),
  function(treatment, control, genome=NULL, fdr=5e-2, eps=1e-5, iterations=10,
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

    #Get classes vector
    classes <- as.integer(rep(NA, length(lnqvals)))
    classes[which(lnqvals < log(fdr))] <-
      apply(fit$lnpost[which(lnqvals < log(fdr)),c(1,3)],1,which.max)

    o <- new("NormRFit", type="enrichR", n=length(treatment), ranges=genome,
             k=3L, B=2L, map=fit$map$map,
             counts=list(as.integer(fit$map$values[1,]),
                         as.integer(fit$map$values[2,])),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             lnqvals=lnqvals, classes=classes)

    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=T)
    }
    return(o)
})
#' \code{diffR}: Difference calling between \code{treatment} (ChIP-seq 1) and
#' \code{control} (ChIP-seq 2) for two given bam filepaths. The given
#' \code{genome} defines chromosomes which will be binned according to
#' \code{countConfig}.
#' @aliases diffR
#' @aliases differenceCall
#' @rdname normr-methods
#' @export
setMethod("diffR", signature("character", "character", "data.frame"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           fdr=5e-2, eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig, procs)
    return(diffR(countsGr$counts[[1]][1], countsGr$counts[[2]][1],
      countsGr$gr, fdr, eps, iterations, procs, verbose))
})
#' \code{diffR}: Difference calling between \code{treatment} (ChIP-seq 1) and
#' \code{control} (ChIP-seq 2) for two given bam filepaths. The given
#' \code{genome} defines a UCSC genome identifier (e.g. "hg19") that is used for
#' retrieval of chromosome annotation.  Chromosomes will be binned according to
#' \code{countConfig}.
#' @aliases diffR
#' @aliases differenceCall
#' @rdname normr-methods
setMethod("diffR", signature("character", "character", "character"),
  function(treatment, control, genome="", countConfig=countConfigSingleEnd(),
           fdr=5e-2, eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome)
    return(diffR(treatment, control, genome, countConfig, fdr, eps, iterations,
      procs, verbose))
})

#' @export
setGeneric("regimeR", function(treatment, control, genome, models, ...)
  standardGeneric("regimeR"))
#' \code{regimeR}: Enrichment regime calling for \code{treatment} (ChIP-seq) over
#' \code{control} (Input) for two given read count vectors and an optional
#' \link{GRanges} object that defines the original regions.
#' @aliases regimeR
#' @aliases regimeCall
#' @rdname normr-methods
#' @export
setMethod("regimeR",
          signature("integer", "integer", "GenomicRangesOrNULL", "integer"),
  function(treatment, control, genome=NULL, models=3L, fdr=5e-2, eps=1e-5,
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

    #Get classes vector
    classes <- as.integer(rep(NA, length(lnqvals)))
    classes[which(lnqvals < log(fdr))] <-
      apply(fit$lnpost[which(lnqvals < log(fdr)),2:models],1,which.max)

    #NormRFit-class object
    o <- new("NormRFit", type="regimeR", n=length(treatment), ranges=genome,
             k=models, B=1L, map=fit$map$map,
             counts=list(as.integer(fit$map$values[1,]),
                         as.integer(fit$map$values[2,])),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             lnqvals=lnqvals, classes=classes)

    #Print logging information
    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=T)
    }
    return(o)
})
#' \code{regimeR}: Enrichment regime calling for \code{treatment} (ChIP-seq) over
#' \code{control} (Input) for two given bam filepaths. The given \code{genome}
#' defines chromosomes which will be binned according to \code{countConfig}.
#' @aliases regimeR
#' @aliases regimeCall
#' @rdname normr-methods
#' @export
setMethod("regimeR", signature("character", "character", "data.frame", "integer"),
  function(treatment, control, genome, models=3L,
           countConfig=countConfigSingleEnd(), fdr=5e-2, eps=1e-5,
           iterations=10, procs=1L, verbose=TRUE) {
    if (models <= 2) stop("invalid number of models specified")
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig, procs)
    return(regimeR(countsGr$counts[[1]][1], countsGr$counts[[2]][1],
      countsGr$gr, models, fdr, eps, iterations, procs, verbose))
})
#' \code{regimeR}: Enrichment regime calling for \code{treatment} (ChIP-seq) over
#' \code{control} (Input) for two given bam filepaths. The given \code{genome}
#' defines a UCSC genome identifier (e.g. "hg19") that is used for retrieval
#' of chromosome annotation. Chromosomes will be binned according to
#' \code{countConfig}.
#' @aliases regimeR
#' @aliases regimeCall
#' @rdname normr-methods
setMethod("regimeR", signature("character", "character", "character", "integer"),
  function(treatment, control, genome="", models=3L,
           countConfig=countConfigSingleEnd(), fdr=5e-2, eps=1e-5,
           iterations=10, procs=1L, verbose=TRUE) {
    if (models <= 2) stop("invalid number of models specified")
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome)
    return(regimeR(treatment, control, genome, models, countConfig, fdr, eps,
      iterations, procs, verbose))
})

