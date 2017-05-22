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

#' Enrichment, Difference and Regime Calling in ChIP-seq data.
#'
#' A correct background estimation is crucial for calling enrichment and
#' differences in ChIP-seq data. \code{\link{normR}} provides robust
#' normalization and difference calling in ChIP-seq and alike data. In brief, a
#' binomial mixture model with a given number of components is fit to read
#' count data for a treatment and control experiment. Therein, computational
#' performance is improved by fitting a log-space model via Expectation
#' Maximization in C++. Convergence is achieved by a threshold on
#' the minimum change in model loglikelihood. After the model fit has
#' converged, a robust background estimate is obtained. This estimate accounts
#' for the effect of enrichment in certain regions and, therefore, represents
#' an appropriate null hypothesis. This robust background is used to identify
#' significantly enriched or depleted regions with respect to control.
#' Moreover, a standardized enrichment for each bin is calculated based on the
#' fitted background component. For convenience, read count vectors can be
#' obtained directly from bam files when a compliant chromosome annotation is
#' given.  Please refer to the individual documentations of functions for
#' enrichment calling (\code{\link{enrichR}}), difference calling
#' (\code{\link{diffR}}) and enrichment regime calling (\code{\link{regimeR}}).
#'
#' Available functions are
#'
#' \code{\link{enrichR}}: Enrichment calling between \code{treatment}
#' (\emph{e.g.} ChIP-seq) and \code{control} (\emph{e.g.} Input).
#'
#' \code{\link{diffR}}: Difference calling between \code{treatment}
#' (\emph{e.g.} ChIP-seq condition 1) and \code{control} (\emph{e.g.} ChIP-seq
#' condition 2).
#'
#' \code{\link{regimeR}}: Enrichment regime calling between \code{treatment}
#' (\emph{e.g.} ChIP-seq) and \code{control} (\emph{e.g.} Input) with a
#' given number of model components. For example, 3 regimes recover background,
#' broad and peak enrichment.
#'
#' The computational performance is improved by fitting a log-space model in
#' C++. Parallization is achieved in C++ via OpenMP (\url{http://openmp.org}).
#'
#' @seealso \code{\link{NormRFit-class}} for functions on accessing and
#' exporting the normR fit. \code{\link{NormRCountConfig-class}} for
#' configuration of the read counting procedure (binsize, mapping quality,...).
#'
#' @name normR
#' @aliases normr PeakCalling DifferentialPeakCalling
#'
#' @example inst/examples/methods_example.R
#'
#' @import parallel
#' @import grDevices
#' @import methods
#' @import bamsignals
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @import Rcpp
#' @importFrom qvalue qvalue
#' @importFrom rtracklayer export
#'
#' @include NormRFit.R
#' @include NormRCountConfig.R
#'
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#'
#' @docType package
#' @useDynLib normr, .registration=TRUE
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

  #as.character(genome[,1]) to prevent GRanges to set factor levels in seqinfo
  gr <- GenomicRanges::GRanges(as.character(genome[,1]),
                               IRanges(1,as.integer(genome[,2])))
  counts <- parallel::mcmapply(
    bamsignals::bamProfile, bampath=c(treatment, control),
    MoreArgs=list(gr=gr, binsize=countConfig@binsize,
                  mapq=countConfig@mapq,
                  shift=countConfig@shift,
                  paired.end=getFilter(countConfig),
                  tlenFilter=countConfig@tlenFilter,
                  filteredFlag=countConfig@filteredFlag,
                  verbose=FALSE),
    mc.cores=procs, SIMPLIFY=FALSE
  )
  counts[[1]] <- unlist(as.list(counts[[1]]))
  counts[[2]] <- unlist(as.list(counts[[2]]))

  #Give bins across the supplied genome
  seqlengths <- genome[,2]
  names(seqlengths) <- genome[,1]
  gr <- GenomicRanges::tileGenome(seqlengths, tilewidth=countConfig@binsize,
                                  cut.last.tile.in.chrom=T)

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

  if (!all(width(gr) == width(gr)[1])) {
    warning("Supplied GenomicRanges have varying widths! Still doing fit...")
  }


  if (verbose) {
    message(paste0("Counting on ", control, " & ", treatment,
                   " for specified GenomicRanges..."))
  }

  counts <- parallel::mcmapply(
    bamsignals::bamCount, bampath=c(treatment, control),
    MoreArgs=list(gr=gr,
                  mapq=countConfig@mapq,
                  shift=countConfig@shift,
                  paired.end=getFilter(countConfig),
                  tlenFilter=countConfig@tlenFilter,
                  filteredFlag=countConfig@filteredFlag,
                  verbose=FALSE),
    mc.cores=procs, SIMPLIFY=T
  )

  return(list(counts=list(counts[,1], counts[,2]), gr=gr))
}

#' Enrichment Calling on ChIP-seq data in normR with enrichR
#'
#' Enrichment calling between \code{treatment} (ChIP-seq) and \code{control}
#' (Input) in normR is done by fitting background and enrichment
#' simultaenously.  Therefore, a mixture of two binomials is fit to the data
#' with Expectation Maximization (EM). After convergence of the EM, the fitted
#' background component is used to calculate significance for treatment and
#' control count pair. Based on this statistic, user can extract significantly
#' enriched regions with a desired significance level. These regions can be
#' further analyzed within R or exported (see \code{\link{NormRFit-class}}).
#' Furthermore, enrichR calculates a standardized enrichment given the fitted
#' background component. See also Details
#'
#' Supplied count vectors for treatment and control should be of same length
#' and of type \code{integer}.
#'
#' For convenience, read count vectors can be obtained directly from bam files.
#' In this case, please specify a bam file for treatment and control each and a
#' \code{genome}. Bam files should be indexed using samtools (\emph{i.e.}
#' samtools index file file.bai). Furthermore, bam files should contain a valid
#' header with given chromosome names. If \code{genome == NULL}(default),
#' chromosome names will be read from treatment bamheader. Please be aware that
#' bamheader might contain irregular contigs and chrM which influence the fit.
#' Also be sure that treatment and control contain the same chromosomes.
#' Otherwise an error will be thrown. If \code{genome} is a \code{character},
#' \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} is used to
#' resolve this to a valid UCSC genome identifier (see
#' \url{https://genome.ucsc.edu/cgi-bin/hgGateway} for available genomes). In
#' this case, only assembled molecules will be considered (no circular). Please
#' check if your bam files obey this annotation. If \code{genome} is a
#' \code{data.frame}, it represents the chromosome specification. The first
#' column will be used as chromosome ID and the second column will be used as
#' the chromosome lengths. If \code{genome} is a \code{GenomicRanges}, it
#' should contain the equally sized genomic loci to count in, e.g. promoters.
#' The binsize in the supplied NormRCountConfig is ignore in this case.
#'
#' \code{bamCountConfig} is an instance of class \code{\link{NormRCountConfig}}
#' specifying settings for read counting on bam files. You can specify the
#' binsize, minimum mapping quality, shifting of read ends etc.. Please refer
#' to \code{\link{NormRFit-class}} for details.
#'
#' @param treatment An \code{integer} vector of treatment counts or a
#' \code{character} pointing to the treatment bam file. In the latter case an
#' "\code{treatment}.bai" index file should exist in the same folder.
#' @param control An \code{integer} vector of control counts or a
#' \link{character} pointing to the control bam file. In the latter case an
#' "\code{control}.bai" index file should exist in the same folder.
#' @param genome Either \code{NULL} (default), a \code{character} specifying a
#' USCS genome identifier, a \link{data.frame} consisting of two columns or a
#' \link{GenomicRanges} specifying the genomic regions (see Details).
#' @param countConfig A \code{\link{NormRCountConfig}} object specifying bam
#' counting parameters for read count retrieval. See Details.
#' @param procs An \code{integer} giving the number of parallel threads to
#' use.
#' @param verbose A \code{logical} indicating whether verbose output is
#' desired.
#' @param eps A \code{numeric} specifying the T Filter threshold and the
#' threshold for EM convergence, \emph{i.e.} the minimal difference in
#' log-likelihood in two consecutive steps.
#' @param iterations An \code{integer} specifying how many times the EM is
#' initialized with random model parameters.
#' @param minP An \code{integer} controlling the threshold for the T
#' method when filtering low power regions, i.e. regions with low counts.
#' @param ... Optional arguments for the respective implementations of
#' \code{\link{enrichR}}.
#'
#'
#' @return A \code{\link{NormRFit}} container holding results of the fit
#' with type \code{enrichR}.
#'
#' @seealso \code{\link{NormRFit-class}} for functions on accessing and
#' exporting the enrichR fit. \code{\link{NormRCountConfig-class}} for
#' configuration of the read counting procedure (binsize, mapping quality,...).
#'
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#'
#' @example inst/examples/methods_example.R
#'
#' @rdname normR-enrichR
#' @aliases enrichR-generic enrichr enrichmentCall EnrichmentCalling
#' @export
setGeneric("enrichR", function(treatment, control, genome, ...)
  standardGeneric("enrichR"))
#' @rdname normR-enrichR
#'
#' @export
setMethod("enrichR", signature("integer", "integer", "GenomicRanges"),
  function(treatment, control, genome, procs=1L, verbose=TRUE,
           eps=1e-5, iterations=10, minP=5e-2) {
    if (length(treatment) != length(control)) {
      stop("incompatible treatment and control count vectors")
    }
    fit <- normr:::normr_core(control, treatment, 2L, eps, iterations,
                              minP, 0, FALSE, verbose, procs)

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
             amount=as.integer(fit$map$amount),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             thresholdT=fit$Tthreshold, lnqvals=lnqvals,
             classes=as.integer(classes))

    #Print logging information
    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=TRUE)
    }
    return(o)
})
#' @rdname normR-enrichR
#'
#' @export
setMethod("enrichR", signature("character", "character", "GenomicRanges"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
            procs=1L, verbose=TRUE, eps=1e-5, iterations=10,
            minP=5e-2) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharGR(treatment, control, genome, countConfig,
                                 procs, verbose)
    return(enrichR(countsGr$counts[[1]], countsGr$counts[[2]], genome, procs,
      verbose, eps, iterations, minP))
})
#' @rdname normR-enrichR
#'
#' @export
setMethod("enrichR", signature("character", "character", "data.frame"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           procs=1L, verbose=TRUE, eps=1e-5, iterations=10,
           minP=5e-2) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig,
                                 procs, verbose)
    return(enrichR(countsGr$counts[[1]], countsGr$counts[[2]], countsGr$gr,
      procs, verbose, eps, iterations, minP))
})
#' @rdname normR-enrichR
#'
#' @export
setMethod("enrichR", signature("character", "character", "character"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           procs=1L, verbose=TRUE, eps=1e-5, iterations=10,
           minP=5e-2) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome, verbose)
    return(enrichR(treatment, control, genome, countConfig,
      procs, verbose, eps, iterations, minP))
})

#' Difference Calling between conditional ChIP-seq data in normR with diffR
#'
#' Difference calling between \code{treatment} (ChIP-seq 1) and \code{control}
#' (ChIP-seq 2) in normR is done by fitting background and two conditional
#' enrichments simultaenously.  Therefore, a mixture of three binomials is fit
#' to the data with Expectation Maximization (EM). After convergence of the EM,
#' the fitted background component is used to calculate significance for
#' treatment and control count pair. Based on this statistic, user can extract
#' significantly enriched/depleted regions in a condition with a desired
#' significance level.  These regions can be further analyzed within R or
#' exported (see \code{\link{NormRFit-class}}).  Furthermore, diffR
#' calculates a standardized conditional-specific enrichment given the
#' fitted background component. See also Details
#'
#' Supplied count vectors for treatment and control should be of same length
#' and of type \code{integer}.
#'
#' For convenience, read count vectors can be obtained directly from bam files.
#' In this case, please specify a bam file for treatment and control each and a
#' \code{genome}. Bam files should be indexed using samtools (\emph{i.e.}
#' samtools index file file.bai). Furthermore, bam files should contain a valid
#' header with given chromosome names. If \code{genome == NULL}(default),
#' chromosome names will be read from treatment bamheader. Please be aware that
#' bamheader might contain irregular contigs and chrM which influence the fit.
#' Also be sure that treatment and control contain the same chromosomes.
#' Otherwise an error will be thrown. If \code{genome} is a \code{character},
#' \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} is used to
#' resolve this to a valid UCSC genome identifier (see
#' \url{https://genome.ucsc.edu/cgi-bin/hgGateway} for available genomes). In
#' this case, only assembled molecules will be considered (no circular). Please
#' check if your bam files obey this annotation. If \code{genome} is a
#' \code{data.frame}, it represents the chromosome specification. The first
#' column will be used as chromosome ID and the second column will be used as
#' the chromosome lengths. If \code{genome} is a \code{GenomicRanges}, it
#' should contain the equally sized genomic loci to count in, e.g. promoters.
#' The binsize in the supplied NormRCountConfig is ignore in this case.
#'
#' \code{bamCountConfig} is an instance of class \code{\link{NormRCountConfig}}
#' specifying settings for read counting on bam files. You can specify the
#' binsize, minimum mapping quality, shifting of read ends etc.. Please refer
#' to \code{\link{NormRFit-class}} for details.
#'
#' @param treatment An \code{integer} vector of treatment counts or a
#' \code{character} pointing to the treatment bam file. In the latter case an
#' "\code{treatment}.bai" index file should exist in the same folder.
#' @param control An \code{integer} vector of control counts or a
#' \link{character} pointing to the control bam file. In the latter case an
#' "\code{control}.bai" index file should exist in the same folder.
#' @param genome Either \code{NULL} (default), a \code{character} specifying a
#' USCS genome identifier, a \link{data.frame} consisting of two columns or a
#' \link{GenomicRanges} specifying the genomic regions (see Details).
#' @param countConfig A \code{\link{NormRCountConfig}} object specifying bam
#' counting parameters for read count retrieval. See Details.
#' @param procs An \code{integer} giving the number of parallel threads to
#' use.
#' @param verbose A \code{logical} indicating whether verbose output is
#' desired.
#' @param eps A \code{numeric} specifying the T Filter threshold and the
#' threshold for EM convergence, \emph{i.e.} the minimal difference in
#' log-likelihood in two consecutive steps.
#' @param iterations An \code{integer} specifying how many times the EM is
#' initialized with random model parameters.
#' @param minP An \code{integer} controlling the threshold for the T
#' method when filtering low power regions, i.e. regions with low counts.
#' @param ... Optional arguments for the respective implementations of
#' \code{\link{diffR}}.
#'
#' @return A \code{\link{NormRFit}} container holding results of the fit
#' with type \code{diffR}.
#'
#' @seealso \code{\link{NormRFit-class}} for functions on accessing and
#' exporting the diffR fit. \code{\link{NormRCountConfig-class}} for
#' configuration of the read counting procedure (binsize, mapping quality,...).
#'
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#'
#' @example inst/examples/methods_example.R
#'
#' @rdname normR-diffR
#' @aliases diffR-generic diffR diffr differenceCall DifferenceCalling
#' @export
setGeneric("diffR", function(treatment, control, genome, ...)
  standardGeneric("diffR"))
#' @rdname normR-diffR
#'
#' @export
setMethod("diffR", signature("integer", "integer", "GenomicRanges"),
  function(treatment, control, genome, procs=1L, verbose=TRUE, eps=1e-5,
           iterations=10, minP=5e-2) {
    if (length(treatment) != length(control)) {
      stop("incompatible treatment and control count vectors")
    }
    fit <- normr:::normr_core(control, treatment, 3L, eps, iterations, minP,
                              1, TRUE, verbose, procs)

    #T filter is computed for label-switched fit as well
    fit2 <- normr:::normr_core(treatment, control, 3L, eps, iterations, minP,
                               1, TRUE, FALSE, procs)
    fit$filteredT <- intersect(fit$filteredT,
                           which(colSums(fit$map$values) >= fit2$Tthreshold)
    )

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
             amount=as.integer(fit$map$amount),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             thresholdT=max(fit$Tthreshold, fit2$Tthreshold),
             lnqvals=lnqvals, classes=as.integer(classes))

    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=TRUE)
    }
    return(o)
})
#' @rdname normR-diffR
#'
#' @export
setMethod("diffR", signature("character", "character", "GenomicRanges"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
            procs=1L, verbose=TRUE, eps=1e-5, iterations=10, minP=5e-2) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharGR(treatment, control, genome, countConfig,
                                 procs, verbose)
    return(diffR(countsGr$counts[[1]], countsGr$counts[[2]], genome, procs,
       verbose, eps, iterations, minP))
})
#' @rdname normR-diffR
#'
#' @export
setMethod("diffR", signature("character", "character", "data.frame"),
  function(treatment, control, genome, countConfig=countConfigSingleEnd(),
           procs=1L, verbose=TRUE, eps=1e-5, iterations=10, minP=5e-2) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig,
                                 procs, verbose)
    return(diffR(countsGr$counts[[1]], countsGr$counts[[2]], countsGr$gr,
       procs, verbose, eps, iterations, minP))
})
#' @rdname normR-diffR
#'
#' @export
setMethod("diffR", signature("character", "character", "character"),
  function(treatment, control, genome="", countConfig=countConfigSingleEnd(),
           procs=1L, verbose=TRUE, eps=1e-5, iterations=10, minP=5e-2) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome, verbose)
    return(diffR(treatment, control, genome, countConfig, procs, verbose, eps,
      iterations, minP))
})

#' Regime Enrichment Calling for ChIP-seq data in normR with regimeR
#'
#' Regime enrichment calling between \code{treatment} (ChIP-seq) and
#' \code{control} (Input) in normR is done by fitting background and multiple
#' enrichment regimes simultaenously.  Therefore, a mixture of \code{models}
#' binomials is fit to the data with Expectation Maximization (EM). After
#' convergence of the EM, the fitted background component is used to calculate
#' significance for treatment and control count pair. Based on this statistic,
#' user can extract significantly enriched regions with a desired significance
#' level. Regime assignments are done by Maximum A Posteriori. Regions can be
#' further analyzed within R or exported (see \code{\link{NormRFit-class}}).
#' Furthermore, regimeR calculates a standardized enrichment given the fitted
#' background component. For example, 3 regimes discriminate background, broad
#' and peak enrichment. See also Details.
#'
#' Supplied count vectors for treatment and control should be of same length
#' and of type \code{integer}.
#'
#' For convenience, read count vectors can be obtained directly from bam files.
#' In this case, please specify a bam file for treatment and control each and a
#' \code{genome}. Bam files should be indexed using samtools (\emph{i.e.}
#' samtools index file file.bai). Furthermore, bam files should contain a valid
#' header with given chromosome names. If \code{genome == NULL}(default),
#' chromosome names will be read from treatment bamheader. Please be aware that
#' bamheader might contain irregular contigs and chrM which influence the fit.
#' Also be sure that treatment and control contain the same chromosomes.
#' Otherwise an error will be thrown. If \code{genome} is a \code{character},
#' \code{\link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}} is used to
#' resolve this to a valid UCSC genome identifier (see
#' \url{https://genome.ucsc.edu/cgi-bin/hgGateway} for available genomes). In
#' this case, only assembled molecules will be considered (no circular). Please
#' check if your bam files obey this annotation. If \code{genome} is a
#' \code{data.frame}, it represents the chromosome specification. The first
#' column will be used as chromosome ID and the second column will be used as
#' the chromosome lengths. If \code{genome} is a \code{GenomicRanges}, it
#' should contain the equally sized genomic loci to count in, e.g. promoters.
#' The binsize in the supplied NormRCountConfig is ignore in this case.
#'
#' \code{bamCountConfig} is an instance of class \code{\link{NormRCountConfig}}
#' specifying settings for read counting on bam files. You can specify the
#' binsize, minimum mapping quality, shifting of read ends etc.. Please refer
#' to \code{\link{NormRFit-class}} for details.
#'
#' @param treatment An \code{integer} vector of treatment counts or a
#' \code{character} pointing to the treatment bam file. In the latter case an
#' "\code{treatment}.bai" index file should exist in the same folder.
#' @param control An \code{integer} vector of control counts or a
#' \link{character} pointing to the control bam file. In the latter case an
#' "\code{control}.bai" index file should exist in the same folder.
#' @param genome Either \code{NULL} (default), a \code{character} specifying a
#' USCS genome identifier, a \link{data.frame} consisting of two columns or a
#' \link{GenomicRanges} specifying the genomic regions (see Details).
#' @param models An \code{integer} specifying the number of mixture
#' components to fit [\code{\link{regimeR}} only]. Default is \code{3}.
#' @param countConfig A \code{\link{NormRCountConfig}} object specifying bam
#' counting parameters for read count retrieval. See Details.
#' @param procs An \code{integer} giving the number of parallel threads to
#' use.
#' @param verbose A \code{logical} indicating whether verbose output is
#' desired.
#' @param eps A \code{numeric} specifying the T Filter threshold and the
#' threshold for EM convergence, \emph{i.e.} the minimal difference in
#' log-likelihood in two consecutive steps.
#' @param iterations An \code{integer} specifying how many times the EM is
#' initialized with random model parameters.
#' @param minP An \code{integer} controlling the threshold for the T
#' method when filtering low power regions, i.e. regions with low counts.
#' @param ... Optional arguments for the respective implementations of
#' \code{\link{regimeR}}.
#'
#' @return A \code{\link{NormRFit}} container holding results of the fit
#' with type \code{regimeR}.
#'
#' @seealso \code{\link{NormRFit-class}} for functions on accessing and
#' exporting the regimeR fit. \code{\link{NormRCountConfig-class}} for
#' configuration of the read counting procedure (binsize, mapping quality,...).
#'
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#'
#' @example inst/examples/methods_example.R
#'
#' @rdname normR-regimeR
#' @aliases regimeR-generic regimeR regimer regimeCall RegimeCalling
#' @export
setGeneric("regimeR", function(treatment, control, genome, models, ...)
  standardGeneric("regimeR"))

#' @rdname normR-regimeR
#'
#' @export
setMethod("regimeR",
          signature("integer", "integer", "GenomicRanges", "numeric"),
  function(treatment, control, genome, models=3, procs=1L, verbose=TRUE,
           eps=1e-5, iterations=10, minP=5e-2) {
    if (models <= 2) stop("invalid number of models specified")
    if (length(treatment) != length(control)) {
      stop("incompatible treatment and control count vectors")
    }
    models = as.integer(models)
    fit <- normr:::normr_core(control, treatment, models, eps, iterations,
                              minP, 0, FALSE, verbose, procs)

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
             amount=as.integer(fit$map$amount),
             names=c("treatment", "control"), thetastar=fit$qstar,
             theta=exp(fit$lntheta), mixtures=exp(fit$lnprior), lnL=fit$lnL,
             eps=eps, lnposteriors=fit$lnpost, lnenrichment=fit$lnenrichment,
             lnpvals=fit$lnpvals, filteredT=fit$filteredT,
             thresholdT=fit$Tthreshold, lnqvals=lnqvals,
             classes=as.integer(classes))

    #Print logging information
    if (verbose) {
      message("\n\n+++ OVERALL RESULT ++++\n")
      summary(o, print=TRUE)
    }
    return(o)
})
#' @rdname normR-regimeR
#'
#' @export
setMethod("regimeR",
          signature("character", "character", "GenomicRanges", "numeric"),
  function(treatment, control, genome, models=3,
           countConfig=countConfigSingleEnd(), procs=1L, verbose=TRUE,
           eps=1e-5, iterations=10, minP=5e-2) {
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharGR(treatment, control, genome, countConfig,
                                 procs, verbose)
    return(regimeR(countsGr$counts[[1]], countsGr$counts[[2]], genome, models,
      procs, verbose, eps, iterations, minP))
})
#' @rdname normR-regimeR
#' @export
setMethod("regimeR",
          signature("character", "character", "data.frame", "numeric"),
  function(treatment, control, genome, models=3,
           countConfig=countConfigSingleEnd(), procs=1L, verbose=TRUE,
           eps=1e-5, iterations=10, minP=5e-2) {
    if (models <= 2) stop("invalid number of models specified")
    treatment <- path.expand(treatment); control <- path.expand(control)
    countsGr <- handleCharCharDf(treatment, control, genome, countConfig,
                                 procs, verbose)
    return(regimeR(countsGr$counts[[1]], countsGr$counts[[2]], countsGr$gr,
      models, procs, verbose, eps, iterations, minP))
})
#' @rdname normR-regimeR
#' @export
setMethod("regimeR",
          signature("character", "character", "character", "numeric"),
  function(treatment, control, genome="", models=3,
           countConfig=countConfigSingleEnd(), procs=1L, verbose=TRUE,
           eps=1e-5, iterations=10, minP=5e-2) {
    if (models <= 2) stop("invalid number of models specified")
    treatment <- path.expand(treatment); control <- path.expand(control)
    genome <- handleCharCharChar(treatment, control, genome, verbose)
    return(regimeR(treatment, control, genome, models, countConfig, procs,
      verbose, eps, iterations, minP))
})
