#' Normalization and difference calling for Next Generation Sequencing (NGS)
#' experiments via joint multinomial modeling
#'
#' Functions for normalization and difference calling in NGS experiment setups.
#' A binomial mixture model with a given number of components is fit and used 
#' for identifying enriched or depleted regions in two given data tracks.
#' Log-space multinomial model is fit by Expectation maximization in C/C++.
#'
#' @name diffr
#' @import Rcpp
#' @import IRanges
#' @import GenomicRanges
#' @import parallel
#' @import bamsignals
#' @docType package
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#' @useDynLib diffr
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
#' @param treatment A \link{numeric} vector of treatment counts or a 
#' \link{character} pointing to the treatment bam file. In the latter case an
#' \code{treatment}.bai index file should exist in the same folder. Should be
#' consistent with control.
#' @param control A \link{numeric} vector of control counts or a 
#' \link{character} pointing to the control bam file. In the latter case an
#' \code{control}.bai index file should exist in the same folder. Should be
#' consistent with treatment.
#' @param genome A \link{data.frame} consisting of two columns. First column 
#' gives sequence names. Second column gives chromosome lengths. This 
#' information can be retrieved via the UCSC chromSizes table.
#' @param bin.size Width of genomic bins in bp.
#' @param models Number of model components.
#' @param eps Threshold for EM convergence.
#' @param p.values Flag for P value computation.
#' @param procs Number of threads to use in parallel::mclapply()
#' @param mapqual discard reads with mapping quality strictly lower than this 
#' parameter. The value 0 ensures that no read will be discarded, the value 254
#' that only reads with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs. 
#' This can be handy in the analysis of chip-seq data.
#' @param paired.end A logical value indicating whether the bampath contains 
#' paired-end sequencing output. In this case, only first reads in proper mapped
#' pairs are considered (FLAG 66).
#' @param paired.end.midpoint A logical value indicating whether the fragment 
#' middle points of each fragment should be counted. Therefore it uses the tlen
#' information from the given bam file (MidPoint = fragment_start + int( 
#' abs(tlen) / 2) )). For even tlen, the given midpoint is the round half down 
#' real midpoint.
#' @param verbose A logical value indicating whether verbose output is desired
#'
#' @return a \link{list} with the following elements:
#' 	\item{posterior}{a matrix containing posteriors for model components.}
#' 	\item{foldchange}{a matrix containing foldchanges for model components
#'  (i=2...models) to background component (i=1).}
#'  \item{fit}{
#'    Result of the multinomial fit. \code{qstar} is naive estimate of 
#'    background intensity. \code{theta} gives binomial mixture model 
#'    parameters. \code{prior} gives binomial mixture model priors.\code{lnL}
#'    gives lLn Likelihood trace.
#'  }
#'  \item{p.values}{-Log10 P-values for enrichment over background for model 
#'  components (i=2..models).}
#'
#' @examples
#' \dontrun{
#' #bam file input
#' normalize( "H3K4me3.bam", "Input.bam" )
#' #count vector input
#' normalize( H3K4me3.counts, Input.counts )
#' #Normalize a H3K4me3 experiment from bam input with 2 enrichment regimes
#' normalize( "H3K4me3.bam", "Input.bam", bin.size=150, models=3 )
#' }
#' 
#' @export
normalize <- function( treatment, 
 					   control,  
					   genome, 
					   bin.size=300, 
					   models=2, 
					   eps=.001,
						 p.values=T,
					   procs=1, 
					   mapqual=20, 
					   shift=0,
					   paired.end=F,
					   paired.end.midpoint=T,
					   verbose=T) {
	# construct GRanges
	if( is.null(genome) ) stop( "No genome specification given. Please provide chromSizes data frame.\n")
	gr <- bin.genome(genome, bin.size)

	# check if treatment and control give bamfiles or counts
	counts = NULL
	bam.files = NULL
	if ( class(treatment) == "character") {
		#counting only possible on indexed bamfiles
		if( !file.exists(paste(treatment, ".bai", sep="")) ) stop( "No index file for", treatment, "found.\n")
		bam.files = c(bam.files, treatment)
	} else if ( class(treatment) == "numeric" | class(treatment) == "integer" ) {
		counts[[1]] = treatment
	}
	if ( class(control) == "character") {
		#counting only possible on indexed bamfiles
		if( !file.exists(paste(control, ".bai", sep="")) ) stop( "No index file for", control, "found.\n")
		bam.files = c(bam.files, control)
	} else if ( class(control) == "numeric" | class(control) == "integer" ) {
		counts[[2]] = control
	}

	# count in GRanges with bamsignals::count if necessary
	if ( is.null(counts) ) {
		counts <- processByChromosome( bam.files=c(treatment, control), gr=gr, procs=procs, bamsignals.function=count, mapqual=mapqual, shift=shift, paired.end=paired.end, paired.end.midpoint=paired.end.midpoint, verbose=verbose)
	}

	if ( length(counts[[1]]) != length(counts[[2]]) ) stop( "Incompatible lengths of treatment and control. Please provide compatible numeric arrays.\n" )

	# fit multinomial model
	result <- em(counts[[1]], counts[[2]], models, eps, verbose)

	# P-values
	if (p.values) {
		if (verbose) {
			cat( "... computing P-values\n" )
		}
		result$log10.p.values     <- matrix(0, nrow=nrow(result$posterior), ncol=(models-1))
		result$log10.p.values[,1] <- -logSum( cbind(pbinom( counts[[1]], counts[[1]]+counts[[2]], result$fit$theta[1], lower.tail=F, log.p=T), dbinom( counts[[1]], counts[[1]]+counts[[2]], result$fit$theta[1], log=T)) )/log(10)
	} else {
		result$log10.p.values <- NULL
	}

	# some logging
	if (verbose) {
		cat(models , "-component multinomial mixture model for treatment with control:\n",
				"LogLik =", tail(result$fit$lnL, 1), ", runs =", 10*length(result$fit$lnL), "\n", 
				"\tq*=",   format( result$fit$qstar, 2, 2), "\n",
				"\tq =",   format( result$fit$theta, 2, 2), "\n",
				"\ttransitions =", format( result$fit$prior, 2, 2), "\n")
		if(p.values) {
			cat("\tenriched (P-value     <= 0.05) =", length( which(result$log10.p.values[,1] > -log10(0.05))), "\n",
					"\tenriched (adj P-value <= 0.05) =", length( which( p.adjust(10^(-result$log10.p.values[,1]), method="BH") < 0.05)), "\n")
		}
	}

	( result )
}
#' Call differences in two NGS experiments.
#'
#' This function implements difference calling of two NGS experiments by joint
#' multinomial modeling. A given number of model components are fit 
#' simultaenously via an efficient Expectation Maximization implementation in 
#' C++. Convergence is achieved by a threshold on the minimum change in model 
#' ln likelihood. The background component is assumed to be the component 
#' with lowest mean and is used to compute two treatments fold change values 
#' and p-values for statistical significance of enrichment.
#'
#' @param treatment.1 A \link{numeric} vector of treatment.1 counts or a 
#' \link{character} pointing to the treatment.1 bam file. In the latter case an
#' \code{treatment.1}.bai index file should exist in the same folder. Should be
#' consistent with treatment.2.
#' @param treatment.2 A \link{numeric} vector of treatment.2 counts or a 
#' \link{character} pointing to the treatment.2 bam file. In the latter case an
#' \code{treatment.2}.bai index file should exist in the same folder. Should be
#' consistent with treatment.1.
#' @param genome A \link{data.frame} consisting of two columns. First column 
#' gives sequence names. Second column gives chromosome lengths. This 
#' information can be retrieved via the UCSC chromSizes table.
#' @param bin.size Width of genomic bins in bp.
#' @param models Number of model components. Note that this should assemble to
#' 1 background components and at least 2 treatment components. Default: 4
#' @param eps Threshold for EM convergence.
#' @param p.values Flag for P value computation.
#' @param procs Number of threads to use in parallel::mclapply()
#' @param mapqual discard reads with mapping quality strictly lower than this 
#' parameter. The value 0 ensures that no read will be discarded, the value 254
#' that only reads with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs. 
#' This can be handy in the analysis of chip-seq data.
#' @param paired.end A logical value indicating whether the bampath contains 
#' paired-end sequencing output. In this case, only first reads in proper mapped
#' pairs are considered (FLAG 66).
#' @param paired.end.midpoint A logical value indicating whether the fragment 
#' middle points of each fragment should be counted. Therefore it uses the tlen
#' information from the given bam file (MidPoint = fragment_start + int( 
#' abs(tlen) / 2) )). For even tlen, the given midpoint is the round half down 
#' real midpoint.
#' @param verbose A logical value indicating whether verbose output is desired
#'
#' @return a \link{list} with the following elements:
#' 	\item{posterior}{a matrix containing posteriors for model components.}
#' 	\item{foldchange}{a matrix containing foldchanges for model components
#'  (i=2...models) to background component (i=1).}
#'  \item{fit}{
#'    Result of the multinomial fit. \code{qstar} is naive estimate of 
#'    background intensity. \code{theta} gives binomial mixture model 
#'    parameters. \code{prior} gives binomial mixture model priors.\code{lnL}
#'    gives lLn Likelihood trace.
#'  }
#'  \item{p.values}{-Log10 P-values for enrichment over background for model 
#'  components (i=2..models).}
#'
#' @examples
#' \dontrun{
#' #bam file input
#' normalize( "H3K4me3.bam", "H3K27me3.bam" )
#' #count vector input
#' normalize( H3K4me3.counts, H3K27me3.counts )
#' #Differences of H3K4me3 and H3K27me3 from bam input with 2 enrichment regimes
#' normalize( "H3K4me3.bam", "H3K27me3.bam", bin.size=150, models=5 )
#' }
#' 
#' @export
diffcall <- function( treatment.1, 
					  treatment.2, 
					  genome, 
					  bin.size=300, 
					  models=4, 
					  eps=.001,
						p.values=T,
					  procs=1, 
					  mapqual=20, 
					  shift=0,
					  paired.end=F,
					  paired.end.midpoint=T,
					  verbose=T) {
	# construct GRanges
	if( is.null(genome) ) stop( "No genome specification given. Please provide chromSizes data frame.\n")
	gr <- bin.genome(genome, bin.size)

	# check if treatment.1 and treatment.2 give bamfiles or counts
	counts = NULL
	bam.files = NULL
	if ( class(treatment.1) == "character") {
		#counting only possible on indexed bamfiles
		if( !file.exists(paste(treatment.1, ".bai", sep="")) ) stop( "No index file for", treatment.1, "found.\n")
		bam.files = c(bam.files, treatment.1)
	} else if ( class(treatment.1) == "numeric" | class(treatment.1) == "integer" ) {
		counts[[1]] = treatment.1
	}
	if ( class(treatment.2) == "character") {
		#counting only possible on indexed bamfiles
		if( !file.exists(paste(treatment.2, ".bai", sep="")) ) stop( "No index file for", treatment.2, "found.\n")
		bam.files = c(bam.files, treatment.2)
	} else if ( class(treatment.2) == "numeric" | class(treatment.2) == "integer" ) {
		counts[[2]] = treatment.2
	}

	# count in GRanges with bamsignals::count if necessary
	if ( is.null(counts) ) {
		counts <- processByChromosome( bam.files=c(treatment.1, treatment.2), gr=gr, procs=procs, bamsignals.function=count, mapqual=mapqual, shift=shift, paired.end=paired.end, paired.end.midpoint=paired.end.midpoint, verbose=verbose)
	}

	if ( length(counts[[1]]) != length(counts[[2]]) ) stop( "Incompatible lengths of treatment.1 and treatment.2 Please provide compatible numeric arrays.\n" )

	# fit multinomial model
	result <- em(counts[[1]], counts[[2]], models, eps, verbose)

	# P-values
	if (p.values) {
		if (verbose) {
			cat( "... computing P-values\n" )
		}
		result$log10.p.values     <- matrix(0, nrow=nrow(result$posterior), ncol=(models-1))
		result$log10.p.values[,1] <- -logSum( cbind(pbinom( counts[[1]], counts[[1]]+counts[[2]], result$fit$theta[1], lower.tail=F, log.p=T), dbinom( counts[[1]], counts[[1]]+counts[[2]], result$fit$theta[1], log=T)) )/log(10)
		result$log10.p.values[,2] <- -logSum( cbind(pbinom( counts[[2]], counts[[1]]+counts[[2]], result$fit$theta[1], lower.tail=F, log.p=T), dbinom( counts[[2]], counts[[1]]+counts[[2]], result$fit$theta[1], log=T)) )/log(10)
	} else {
		result$log10.p.values <- NULL
	}

 	# some logging
	if (verbose) {
		cat(models , "-component multinomial mixture model for treatment with control:\n",
				"LogLik =", tail(result$fit$lnL, 1), ", runs =", 10*length(result$fit$lnL), "\n", 
				"\tq*=",   format( result$fit$qstar, 2, 2), "\n",
				"\tq =",   format( result$fit$theta, 2, 2), "\n",
				"\ttransitions =", format( result$fit$prior, 2, 2), "\n")
		if(p.values) {
			cat("\t'treatment.1' enriched (P-value     <= 0.05) =", length( which(result$log10.p.values[,1] > -log10(0.05))), "\n",
					"\t'treatment.1' enriched (adj P-value <= 0.05) =", length( which( p.adjust( 10^(-result$log10.p.values[,1]), method="BH") < 0.05)), "\n",
					"\t'treatment.2' enriched (P-value     <= 0.05) =", length( which(result$log10.p.values[,1] > -log10(0.05))), "\n",
					"\t'treatment.2' enriched (adj P-value <= 0.05) =", length( which( p.adjust( 10^(-result$log10.p.values[,1]), method="BH") < 0.05)), "\n")
		}
	}

	( result )
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
	n <- floor(genome[,2] / bin.size)
	names(n) <- genome[,1]
	bin.sizes <- rep(bin.size, dim(genome)[1])
	names(bin.sizes) <- genome[,1]
	gr <- GRanges()
	for (ch in genome[,1]) {
		gr <- c(gr, GRanges(seqnames=ch, IRanges(start=0:(n[ch] - 1) * bin.sizes[ch] + 1, width=bin.size)))
	}
	gr <- sort( gr )
	GenomeInfoDb::seqlengths(gr) <- genome[,2]
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
#' @param paired.end A logical value indicating whether the bampath contains 
#' paired-end sequencing output. In this case, only first reads in proper mapped
#' pairs are considered (FLAG 66).
#' @param paired.end.midpoint A logical value indicating whether the fragment 
#' middle points of each fragment should be counted. Therefore it uses the tlen
#' information from the given bam file (MidPoint = fragment_start + int( 
#' abs(tlen) / 2) )). For even tlen, the given midpoint is the round half down 
#' real midpoint.
#' @param verbose A logical value indicating whether verbose output is desired
#' @return a list of length(bam.files) with bamsignals.function results
processByChromosome <- function(bam.files, gr, procs, bamsignals.function, mapqual, shift=0, paired.end=F, paired.end.midpoint=T, verbose=F) {
 	x <- mclapply(as.character(unique(seqnames(gr))), function(chunk) {
			 	  cat("[", paste(Sys.time()),"] Counting on chromosome", chunk, "\n")
				  gr.sub <- gr[ seqnames(gr) %in% chunk]
				  lapply( bam.files, bamsignals.function, gr=gr.sub, mapqual=mapqual, shift=shift, paired.end=paired.end, paired.end.midpoint=paired.end.midpoint, verbose=verbose)
 				}, mc.cores=procs)
	lapply( 1:length(bam.files), function( i ) {
				unlist(lapply(x, "[[", i)) 
			}
	)
}

#EOF
