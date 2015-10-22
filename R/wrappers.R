# Copyright (C) 2015 Johannes Helmuth
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
#' Normalization and difference calling for Next Generation Sequencing (NGS)
#' experiments via joint multinomial modeling
#'
#' Functions for normalization and difference calling in NGS experiment setups.
#' A binomial mixture model with a given number of components is fit and used 
#' for identifying enriched or depleted regions in two given data tracks.
#' Log-space multinomial model is fit by Expectation maximization in C++.
#'
#' @name diffr
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import IRanges
#' @import Rcpp
#' @import Rsamtools
#' @import parallel
#' @import bamsignals
#' @docType package
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#' @useDynLib diffr, .registration=TRUE
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
#' @param p.values Flag for P value computation (not implemented)
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
#' @examples
#' \dontrun{
#' #bam file input
#' diffr("H3K4me3.bam", "Input.bam")
#' #count vector input
#' diffr(H3K4me3.counts, Input.counts)
#' #Normalize a H3K4me3 experiment from bam input with 2 enrichment regimes
#' diffr("H3K4me3.bam", "Input.bam", bin.size=150, models=3)
#' #Difference Calling
#' diffr("H3K4me3.bam", "H3K4me2.bam", bin.size=150, models=3)
#' }
#' 
#' @export
diffR <- function(treatment, 
                  control,  
                  genome=NULL, 
                  bin.size=300, 
                  models=2, 
                  eps=.001,
                  p.values=T,
                  procs=1, 
                  mapqual=20, 
                  shift=0,
                  paired.end="ignore",
                  verbose=T) {

	# Check if treatment and control give bampaths or counts
	counts = NULL
	bam.files = NULL
	if (class(treatment) == "character") {
		if(!file.exists(paste(treatment, ".bai", sep=""))) 
			stop("No index file for", treatment, "found.\n")
		bam.files = c(bam.files, treatment)
	} else if (class(treatment) == "numeric" | class(treatment) == "integer") {
		counts[[1]] = treatment
	}
	if (class(control) == "character") {
		if(!file.exists(paste(control, ".bai", sep=""))) stop("No index file for", control, "found.\n")
		bam.files = c(bam.files, control)
	} else if (class(control) == "numeric" | class(control) == "integer") {
		counts[[2]] = control
	}

	# count in GRanges with bamsignals::count if necessary
	if (is.null(counts)) {
		if (!class(genome) %in% c("data.frame", "matrix")) { #no chrom lengths given
			if (is.null(genome)) { #read from bamheader
				header <- Rsamtools::scanBamHeader(treatment)[[1]][[1]]
				genome <- cbind(names(header), header)
			} else { #read from UCSC
				genome <- fetchExtendedChromInfoFromUCSC(genome)
				genome <- genome[which(!genome$circular & genome$SequenceRole == "assembled-molecule"),1:2]
			}
		}
		gr <- bin.genome(genome, bin.size)
		counts <- processByChromosome(bam.files=c(treatment, control), 
		                              gr=gr, 
		                              procs=procs, 
		                              bamsignals.function=bamCount,
		                              mapqual=mapqual,
		                              shift=shift,
		                              paired.end=paired.end,
		                              verbose=verbose)
	}

	if (length(counts[[1]]) != length(counts[[2]])) 
		stop("Incompatible lengths of treatment and control. Please provide compatible numeric arrays.\n")

	# fit multinomial model
	result <- list("ranges"=gr, "control"=counts[[2]], "treatment"=counts[[1]])
	result <- c(result, em(counts[[2]], counts[[1]], models, eps, verbose, procs))

	# P-values
	#if (p.values) {
	#  if (verbose) {
	#    cat( "... computing P-values\n" )
	#  }
	#  result$log10.p     <- matrix(0, nrow=nrow(result$posterior), ncol=models)
	#  result$log10.p[,1] <- -logSum( cbind(pbinom( counts[[1]], counts[[1]]+counts[[2]], result$fit$theta[1], lower.tail=F, log.p=T), dbinom( counts[[1]], counts[[1]]+counts[[2]], result$fit$theta[1], log=T)) )/log(10)
	#  #result$log10.adjp
	#} else {
	#  result$log10.p <- NULL
	#}

	# some logging
	if (verbose) {
		cat(models , "-component multinomial mixture model for treatment with control:\n",
		    "LogLik =", tail(result$lnL, 1), ", runs =", 10*length(result$lnL), "\n", 
		    "\tq*=",   format( result$qstar, 2, 2), "\n",
		    "\tq =",   format( result$theta, 2, 2), "\n",
		    "\ttransitions =", format( result$prior, 2, 2), "\n")
		#if(p.values) {
		#  cat("\tenriched (P-value     <= 0.05) =", length( which(result$log10.p[,1] > -log10(0.05))), "\n",
		#      "\tenriched (adj P-value <= 0.05) =", length( which( p.adjust(10^(-result$log10.p.values[,1]), method="BH") < 0.05)), "\n")
		#}
	}

	(result)
}

#' Create bins for a given genome annotation.
#' 
#' Computes non-overlapping bins for given chromosome names and lengths.
#'
#' @param genome A \link{data.frame} consisting of two columns. First column 
#' gives sequence names. Second column gives chromosome lengths. This 
#' information can be retrieved via the UCSC chromSizes table.
#' @param bin.size The size of the bins in bp.
#' @return \link{GenomicRanges}-object specifying bins for genome
#' 
#' @export
bin.genome <- function(genome, bin.size=300) {
	lens <- as.numeric(genome[,2])
	n <- floor(lens / bin.size)
	names(n) <- genome[,1]
	bin.sizes <- rep(bin.size, dim(genome)[1])
	names(bin.sizes) <- genome[,1]
	gr <- GRanges()
	suppressWarnings( {
		for (ch in genome[,1]) {
			gr <- c(gr, GRanges(seqnames=ch, IRanges(start=0:(n[ch] - 1) * bin.sizes[ch] + 1, width=bin.size)))
		}
	})
	gr <- sort( gr )
	GenomeInfoDb::seqlengths(gr) <- lens
	(gr) 
}

#' Process a list of bam files with a number of processes.
#'
#'
#' @param bam.files A list of strings specifying file system locations for bam 
#' files to count.
#' @param gr A\link{GenomicRanges}-object giving the genomic intervals to count
#' in.
#' @param procs Number of threads to use in parallel::mclapply()
#' @param bamsignals.function The corresponding counting function in bamsignals
#' to call.
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
#' @return a list of length(bam.files) with bamsignals.function results
processByChromosome <- function(bam.files, gr, procs, bamsignals.function, mapqual, shift=0, paired.end="ignore", verbose=F) {
	bam.files <- path.expand(bam.files)
	x <- safe_mclapply(as.character(unique(seqnames(gr))), function(chunk) {
	                   if (verbose)
	                   	cat("[", paste(Sys.time()),"] Counting on chromosome", chunk, "\n")
	                   gr.sub <- gr[ seqnames(gr) %in% chunk]
	                   lapply(bam.files, bamsignals.function, gr=gr.sub, mapqual=mapqual, shift=shift, paired.end=paired.end, verbose=F)
	}, mc.cores=procs)
	lapply(1:length(bam.files), function( i ) {
	       unlist(lapply(x, "[[", i)) 
	}
	)
}

#stolen from epicseg
#make sure that errors occurring in mclapply get propagated
propagateErrors <- function(l){
	for (el in l) if (inherits(el, "try-error")) stop(el)
	l
}
safe_mclapply <- function(...) propagateErrors(mclapply(...))
