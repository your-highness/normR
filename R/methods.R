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
#' @import Rsamtools
#' @import parallel
#' @import bamsignals
#' @docType package
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#' @seealso \code{\link{normr-methods}} for available functions
#' @useDynLib normr, .registration=TRUE
#' @include NormRFit.R 
#' @include BamCountConfig.R
NULL

#' Normalize a NGS experiment with given background data.
#'
#' This function implements the normalization of a treatment track with its
#' given control track by joint multinomial modeling. A given number of model
#' components are fit simultaenously via an efficient Expectation
#' Maximization implementation in C++. Convergence is achieved by a
#' threshold on the minimum change in model ln likelihood. The background
#' component is assumed to be the component with lowest mean and is used
#' to compute treatment fold change and P-Values for statistical significance
#' of enrichment.
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
#' \link{GenomeInfoDb::fetchExtendedChromInfoFromUCSC} (only assembled
#' molecules, circular omitted) or a \link{data.frame} consisting of two
#' columns(1st column=chromosome names, 2nd column=lengths). Chromosome names
#' given should be present in the bam file header. If \code{NULL}, chromosome
#' names will be read from treatment bamheader. Please be aware that bamheader
#' might contain irregular contigs and chrM which influence the diffR fit.
#' (DEFAUlT=NULL).
#' @param bin.size Width of genomic bins in bp.
#' @param models Number of model components.
#' @param eps Threshold for EM convergence.
#' @param procs Number of threads to use
#' @param countConfig A BamCountConfig object specifying bam counting parameters
#' for read count retrieval.
#' @param filename A filename to write data to.
#' @param threshold A threshold on q-values. Set to 1 to write all regions.
#' @param verbose A logical value indicating whether verbose output is desired
#'
#' @return a \link{NormRFit} object 
#'
#' @example inst/examples/methods_example.R
#' @see BamCountConfig-class
#' @see NormRFit-class

#' @export
setGeneric("enrichR", function(treatment, control, genome, ...)
  standardGeneric("enrichR"))
setClassUnion("GenomicRangesOrNULL", c("GenomicRanges", "NULL"))
#' \code{enrichR}: Do enrichment calling between treatment (ChIP-seq) and
#' control (Input) for two given bam files and a genome data.frame
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("integer", "integer", "GenomicRangesOrNULL"),
  function(treatment, control, genome=NULL, eps=1e-5, iterations=10, procs=1L,
           verbose=TRUE) {
    if (length(treatment) != length(control)) {
      stop("incompatible treatment and control count vectors")
    }

    # C++ does computation & construct NormRFit-class object herein
    fit <- normr:::normr_core(control, treatment, 2L, eps, iterations, 1, F, 
                              verbose, procs)

    #Storey's q-value on T filtered P-values
    if (verbose) message("... computing Q-values.")
    idx <- which(fit$map$map %in% fit$filtered)
    qvals <- qvalue(exp(fit$pvals[idx]), eps)
    lnqvals <- rep(NA,length(treatment))
    lnqvals[idx] <- log(qvals$qvalues)
    #FIXME
    lnqvals <- normr:::mapToUniqueWithMap(lnqvals, fit$map)

    #NormRFit-class object
    o <- new("NormRFit", type="enrichR", n=length(treatment), ranges=genome,
             k=2, B=1, map=fit$map$map,
             counts=list(fit$map$values[1,],fit$map$values[2,]),
             names=c(names(treatment), names(control)), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrrichment, 
             lnpvals=fit$lnpvals, filteredT=fit$filtered, lnqvals=qvals)

    #Print logging information
    if (verbose) message("+++ OVERALL RESULT ++++\n\n", summary(obj, print=F))
    return(o)
})
#' \code{enrichR}: Do enrichment calling between treatment (ChIP-seq) and
#' control (Input) for two given bam files and a genome data.frame
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("character", "character", "data.frame"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
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
                                     mapqual=countConfig@mapqual,
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

    return(enrichR(counts[[1]], counts[[2]], gr, countConfig, eps, iterations, 
      procs, verbose))
})
#' \code{enrichR}: Do enrichment calling between treatment (ChIP-seq) and
#' control (Input) for two given bam files and a genome data.frame
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
setMethod("enrichR", signature("character", "character", "character"),
  function(treatment, control, genome="", countConfig=countConfigSingleEnd(),
           eps=1e-5, iterations=10, procs=1L, verbose=TRUE) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    if (!file.exists(treatment)) strop("treatment is not a file")
    if (!file.exists(control)) strop("control is not a file")
    if (treatment == control) stop("treatment and control are identical")

    #UseGenomeInfoDb to fetch information on regular, non-circular chroms
    genome <- fetchExtendedChromInfoFromUCSC(genome)
    idx <- which(!genome$circular & genome$SequenceRole=="assembled-molecule")
    genome <- genome[idx,1:2]

    return(enrichR(treatment, control, genome,countConfig, eps, iterations, 
      procs, verbose))
})

#' @export
setGeneric("recomputeP", function(fit, B) standardGeneric("recomputeP"))
#' \code{recomputeP}: Recompute P-values relative to provided background
#' component B.
#' @aliases computeP
#' @rdname normr-methods
setMethod("recomputeP", signature("NormRFit", "integer"), function(fit, B) {
      stop("not implemented yet")
})

#' @export
setGeneric("exportR", function(fit, filename, format, ...) 
  standardGeneric("exportR"))
#' \code{exportR}: Export results of a normR fit. The provided filetype
#' specifies if significant regions (\code{"bed"}) or calculated enrichment
#' (\code{c("bigWig", "bedGraph")}) is exported
#' @aliases exportR
#' @aliases exportNormRFit
#' @aliases exporter
#' @rdname normr-methods
setMethod("exportR", signature("NormRFit", "character", "character"),
  function(fit, filename, format=c("bed", "bigWig"), threshold=0.05) {
    if (!(format %in% c("bed", "bigWig"))) {
      stop("invalid format parameter")
    }
    if (format == "bed") writeBed(fit,filename,threshold)
    if (format == "bigWig") writeEnrichment(fit,filename)
})
