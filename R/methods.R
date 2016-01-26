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
#' @include NormRFit.R, BamCountConfig.R
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
#' @param mapqual discard reads with mapping quality strictly lower than this
#' parameter. The value 0 ensures that no read will be discarded, the value 254
#' that only reads with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs.
#' This can be handy in the analysis of chip-seq data.
#' @param paired.end a character string indicating how to handle paired-end
#' reads. If \code{paired.end!="ignore"} then only first reads in proper mapped
#' pairs will be consider (i.e. in the flag of the read, the bits in the mask
#' 66 must be all ones). If \code{paired.end=="midpoint"} then the midpoint of a
#' fragment is considered, where \code{mid = fragment_start + int(abs(tlen)/2)},
#' and where tlen is the template length stored in the bam file. For even tlen,
#' the given midpoint will be moved of 0.5 basepairs in the 3' direction.
#' If \code{paired.end=="extend"} then the whole fragment is treated
#' as a single read.
#' @param verbose A logical value indicating whether verbose output is desired
#'
#' @return a \link{list} with the following elements:
#'  \item{posterior}{a matrix containing posteriors for model components.}
#'  \item{fit}{
#'    Result of the multinomial fit. \code{qstar} is naive estimate of
#'    background intensity. \code{theta} gives binomial mixture model
#'    parameters. \code{prior} gives binomial mixture model priors.\code{lnL}
#'    gives lLn Likelihood trace.
#'  }
#'
#' @example inst/examples/methods_example.R

#' @export
setGeneric("enrichR",
           function(treatment, control, genome, ...)
           standardGeneric("enrichR"))
#' \code{enrichR}: Do enrichment calling between treatment (ChIP-seq) and
#' control (Input) for two given bam files and a genome data.frame
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("integer", "integer", "GenomicRanges"),
    function(treatment, control, genome=NULL, eps=1e-5, procs=1, verbose=TRUE) {
      if (length(treatment) != length(control)) {
        stop("invalid treatment and control")
      }

      # C++ does computation & construct NormRFit-class object herein
      fit <- diffr_core(counts[[2]], counts[[1]], 2, eps, verbose, procs)
      obj <- new("NormRFit", type="enrichR", k=2, B=1, eps=eps, ranges=genome,
         names=c(names(treatment), names(control)), thetastar=fit$qstar,
         counts=list("Input"=control, "Treatment"=treatment),
         n=length(treatment), theta=fit$theta, mixtures=fit$prior, lnL=fit$lnL,
         posteriors=fit$post, enrichment=enr, p.vals=fit$p.vals,
         filteredT=fit$filteredidx, q.vals=fit$q.vals)

      #Print logging information
      if (verbose) message(summary(obj, print=F))

      return(obj)
})
#' \code{enrichR}: Do enrichment calling between treatment (ChIP-seq) and
#' control (Input) for two given bam files and a genome data.frame
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname normr-methods
#' @export
setMethod("enrichR", signature("character", "character", "data.frame"),
          function(treatment, control, genome, countConfig=countConfigSingleEnd(),
                   eps=1e-5, procs=1, verbose=TRUE) {
            if(!file.exists(paste(treatment, ".bai", sep=""))) {
              stop("No index file for", treatment, ".\n")
            }
            if(!file.exists(paste(control, ".bai", sep=""))) {
              stop("No index file for", control, ".\n")
            }
            if (NCOL(genome) != 2) stop("invalid genome data.frame")

            require(GenomicRanges)
            gr <- GRanges(genome[,1], IRanges(1,genome[,2]))
            require(parallel)
            counts <- mcmapply(bamsignals::bamProfile, bampath=c(treatment, control),
                               MoreArgs=list(gr=gr, binsize=countConfig@binsize,
                                             mapqual=countConfig@mapqual,
                                             shift=countConfig@shift,
                                             paired.end=getFilter(countConfig),
                                             #Tlen.filter not yet in Bioconductor
                                             #tlen.filter=countConfig@tlenFilter,
                                             verbose=F),
                               mc.cores=procs, SIMPLIFY=F)

            #Give bins across the supplied genome
            gr <- unlist(tile(gr, countConfig@binsize))

            return(enrichR(counts[[1]], counts[[2]], gr, countConfig, eps, procs,
                           verbose))
})
#' \code{enrichR}: Do enrichment calling between treatment (ChIP-seq) and
#' control (Input) for two given bam files and a genome data.frame
#' @aliases enrichR
#' @aliases enrichmentCall
#' @rdname
setMethod("enrichR", signature("character", "character", "character"),
    function(treatment, control, genome, countConfig=countConfigSingleEnd(),
             eps=1e-5, procs=1, verbose=TRUE) {
      #UseGenomeInfoDb to fetch information on regular, non-circular chroms
      genome <- fetchExtendedChromInfoFromUCSC(genome)
      idx <- which(!genome$circular & genome$SequenceRole=="assembled-molecule")
      genome <- genome[idx,1:2]

      return(enrichR(treatment, control, genome, countConfig, eps, procs,
                     verbose))
  }
)




