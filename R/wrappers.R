#' Normalization and difference calling for Next Generation Sequencing (NGS) experiments via joint multinomial
#'
#' Functions for normalization and difference calling in NGS experiment setups. A binomial mixture model with a given number
#' of components is fit and used for identifying enriched or depleted regions in two given data tracks.
#' Log-space multinomial model is fit by Expectation maximization in C/C++.
#'
#' @name diffr
#' @imports Rcpp
#' @imports IRanges
#' @imports GenomicRanges
#' @imports parallel
#' @imports bamsignals
#' @docType package
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#' @useDynLib diffr
NULL
#' Normalize a NGS experiment with given background data.
#'
#' Compute read density in the regions specified by a GenomicRanges object.
#' A read position is always specified by its 5' end, so a read mapping to the reference strand
#' is positioned at its leftmost coordinate, a read mapping to the alternative strand
#' is positioned at its rightmost coordinate. To change that use the \code{shift} parameter
#' or the \code{coverage} function.
#' @param gr GenomicRanges object used to specify the regions
#' @param bampath path to the bam file storing the read. The file must be indexed. 
#' If a range is on the negative strand the profile will be reverse-complemented.
#' @param If the value is set to 1, the method will return basepair-resolution read densities,
#' for bigger values the density profiles will be binned (and the memory requirements
#' will scale accordingly). 
#' @param mapqual discard reads with mapping quality strictly lower than this parameter.
#' The value 0 ensures that no read will be discarded, the value 254 that only reads
#' with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs. This can
#' be handy in the analysis of chip-seq data.
#' @param ss produce a strand-specific profile or ignore the strand of the read. This option
#' does not have any effect on the strand of the region. Strand-specific profiles are
#' twice as long then strand-independent profiles.
#' @param format attempts to find a suitable matrix/array format for the count vector. 
#' if the profile is strand-specific one dimension will correspond to sense
#' and anti-sense strand, if the ranges have all the same width one dimension
#' will correspond to the range number.
#' @param paired.end a logical value indicating whether the bampath contains paired-end 
#' sequencing output. In this case, only first reads in proper mapped pairs are considered 
#' (FLAG 66).
#' @param paired.end.midpoint a logical value indicating whether the fragment middle 
#' points of each fragment should be counted. Therefore it uses the tlen information from
#' the given bam file (MidPoint = fragment_start + int( abs(tlen) / 2) )). For even tlen, 
#' the given midpoint is the round half down real midpoint.
#' @param paired.end.max.frag.length an integer indicating which fragments should be 
#' considered in paired-end sequencing data. Default value of 1,000 bases is generally
#' a good pick.
#' @return a list with the following arguments:
#' 	\item{counts}{the vector containing the read counts. This will be formatted
#' 	into a matrix or an array depending on whether the profile is strand-specific
#' 	and whether the ranges have all the same length.}
#' 	\item{starts, ends}{Vectors defining the boundaries of the count vector. 
#' 	To extract counts relative to the i-th range, use 
#'		\code{as.numeric(counts)[starts[i]:ends[i]]}, 
#' 	or the \code{getSignal} function to preserve the formatting.}
#'		\item{format}{This element is present if pu$counts is formatted
#' 	differently than a simple vector and it describes the formatting.}
#' @export
normalize <- function( treatment, 
					   control, 
					   genome, 
					   bin.size=300, 
					   models=2, 
					   eps=.001,
					   procs=1, 
					   mapqual=20, 
					   shift=0,
					   paired.end=F,
					   paired.end.midpoint=T,
					   verbose=T) {
	# construct GRanges
	gr <- bin.genome(genome, bin.size)

	# check if treatment and control give bamfiles or counts
	counts = NULL
	bam.files = NULL
	if ( class(treatment) == "character") {
		#counting only possible on indexed bamfiles
		if( !file.exists(paste(treatment, ".bai", sep="")) ) stop( "No index file for", treatment, "found.\n")
		bam.files = c(bam.files, treatment)
	} else if ( class(treatment) == "numeric" ) {
		counts[[1]] = treatment
	}
	if ( class(control) == "character") {
		#counting only possible on indexed bamfiles
		if( !file.exists(paste(control, ".bai", sep="")) ) stop( "No index file for", control, "found.\n")
		bam.files = c(bam.files, control)
	} else if ( class(control) == "numeric" ) {
		counts[[2]] = control
	}

	# count in GRanges with bamsignals::count if necessary
	if ( is.null(counts) ) {
		counts <- processByChromosome( bam.files=c(treatment, control), gr=gr, procs=procs, bamsignals.function=count, mapqual=mapqual, shift=shift, paired.end=paired.end, paired.end.midpoint=paired.end.midpoint, verbose=verbose)
	}

	if ( length(counts[[1]]) != length(counts[[2]]) ) stop( "Incompatible lengths of treatment and control. Please provide compatible numeric arrays.\n"

	# fit multinomial model
	result            <- em(counts[[1]], counts[[2]], models, eps, verbose)
	result$p.vals     <- matrix(0, nrow=nrow(post), ncol=(models-1))
	result$p.vals[,1] <- -logSum( cbind(pbinom( counts[[1]], counts[[1]]+counts[[2]], result$theta[1], lower.tail=F, log.p=T), dbinom( counts[[1]], counts[[1]]+counts[[2]], result$theta[1], log=T)) )/log(10)

	# some logging
	if (verbose) {
		cat( model, "-component multinomial mixture model for treatment='", treatment, "' with control='", control, "':\n",
			 "LogLik =", tail(result$fit$lnL, 1), ", runs =", 10*length(result$fit$lnL), "\n", 
			 "\tq*=",   format(sum(s[idx])/(sum(s[idx]+r[idx])),2,2), "\n",
			 "\tq =",   format( result$fit$theta, 2, 2), "\n",
			 "\tmix =", format( result$fit$prior, 2, 2), "\n",
			 "\tenriched (P-value     <= 0.05) =", length( which(result$p.vals[,1] > -log10(0.05))), "\n",
			 "\tenriched (adj P-value <= 0.05) =", length( which( p.adjust( 10^(-result$p.pvals[,1]), method="BH") < 0.05)), "\n"
		)
	}

	( list( cbind("t"=fc[,1], "p"=p.vals[,1], "p.adj"=p.adj[[1]], "post"=post[,2]), "theta"=theta, "mix"=prior ) )
	( result )
}

#' @export
diffcall <- function( treatment, 
					   control, 
					   genome, 
					   bin.size=300, 
					   models=3, 
					   eps=.001,
					   procs=1, 
					   mapqual=20, 
					   shift=0,
					   paired.end=F,
					   paired.end.midpoint=T,
					   verbose=T) {
	# construct GRanges
	gr <- bin.genome(genome, bin.size)

	# check if treatment and control give bamfiles or counts
	counts = NULL
	bam.files = NULL
	if ( class(treatment) == "character") {
		#counting only possible on indexed bamfiles
		if( !file.exists(paste(treatment, ".bai", sep="")) ) stop( "No index file for", treatment, "found.\n")
		bam.files = c(bam.files, treatment)
	} else if ( class(treatment) == "numeric" ) {
		counts[[1]] = treatment
	}
	if ( class(control) == "character") {
		#counting only possible on indexed bamfiles
		if( !file.exists(paste(control, ".bai", sep="")) ) stop( "No index file for", control, "found.\n")
		bam.files = c(bam.files, control)
	} else if ( class(control) == "numeric" ) {
		counts[[2]] = control
	}

	# count in GRanges with bamsignals::count if necessary
	if ( is.null(counts) ) {
		counts <- processByChromosome( bam.files=c(treatment, control), gr=gr, procs=procs, bamsignals.function=count, mapqual=mapqual, shift=shift, paired.end=paired.end, paired.end.midpoint=paired.end.midpoint, verbose=verbose)
	}

	if ( length(counts[[1]]) != length(counts[[2]]) ) stop( "Incompatible lengths of treatment and control. Please provide compatible numeric arrays.\n"

	# fit multinomial model
	result            <- em(counts[[1]], counts[[2]], models, eps, verbose)
	result$p.vals     <- matrix(0, nrow=nrow(post), ncol=(models-1))
	result$p.vals[,1] <- -logSum( cbind(pbinom( counts[[1]], counts[[1]]+counts[[2]], result$theta[1], lower.tail=F, log.p=T), dbinom( counts[[1]], counts[[1]]+counts[[2]], result$theta[1], log=T)) )/log(10)
	result$p.vals[,2] <- -logSum( cbind(pbinom( counts[[2]], counts[[1]]+counts[[2]], result$theta[1], lower.tail=F, log.p=T), dbinom( counts[[2]], counts[[1]]+counts[[2]], result$theta[1], log=T)) )/log(10)

	# some logging
	if (verbose) {
		cat( model, "-component multinomial mixture model for treatment='", treatment, "' with control='", control, "':\n",
			 "LogLik =", tail(result$fit$lnL, 1), ", runs =", 10*length(result$fit$lnL), "\n", 
			 "\tq*=",   format(sum(s[idx])/(sum(s[idx]+r[idx])),2,2), "\n",
			 "\tq =",   format( result$fit$theta, 2, 2), "\n",
			 "\tmix =", format( result$fit$prior, 2, 2), "\n",
			 "\t'treatment' enriched (P-value     <= 0.05) =", length( which(result$p.vals[,1] > -log10(0.05))), "\n",
			 "\t'treatment' enriched (adj P-value <= 0.05) =", length( which( p.adjust( 10^(-result$p.pvals[,1]), method="BH") < 0.05)), "\n"
			 "\t'control' enriched (P-value     <= 0.05) =", length( which(result$p.vals[,1] > -log10(0.05))), "\n",
			 "\t'control' enriched (adj P-value <= 0.05) =", length( which( p.adjust( 10^(-result$p.pvals[,1]), method="BH") < 0.05)), "\n"
		)
	}

	( list( cbind("t"=fc[,1], "p"=p.vals[,1], "p.adj"=p.adj[[1]], "post"=post[,2]), "theta"=theta, "mix"=prior ) )
	( result )

			( list( cbind("t.s"=fc[,1], "p.s"=p.vals[,1], "p.s.adj"=p.adj[[1]], "post.s"=post[,2], "t.r"=fc[,2], "p.r"=p.vals[,2], "p.r.adj"=p.adj[[2]], "post.r"=post[,3]), "theta"=theta, "mix"=prior ) )
}

#' Create bins for a given genome annotation.
#' 
#' Computes non-overlapping bins for given chromosome names and lengths.
#'
#' @param genome A \link{data.frame} consisting of two columns. First column 
#' gives sequence names. Second column gives chromosome lengths.
#' @param bin.size The size of the bins in bp.
#' @return a GenomicRanges object specifying bins
#'
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
	seqlengths(gr) <- genome[,2]
	(gr) 
}

#' Process a list of bam files with a number of processes.
#'
#'
#' @param bam.files A list of strings specifying file system locations for 
#' bam files to count.
#' @param gr \code{GenomicRanges}-object giving the genomic intervals to count in.
#' @param bamsignals.function The corresponding counting function in bamsignals to call.
#' @param mapqual discard reads with mapping quality strictly lower than this parameter.
#' The value 0 ensures that no read will be discarded, the value 254 that only reads
#' with the highest possible mapping quality will be considered.
#' @param procs Number of processors to use.
processByChromosome <- function(bam.files, gr, procs, bamsignals.function, mapqual, shift=0, paired.end=F, paired.end.midpoint=T, verbose=F) {
	x <- mclapply(as.character(unique(seqnames(gr))), function(chunk) {
			 	  cat("[", paste(Sys.time()),"] Counting on chromosome", chunk, "\n")
				  gr.sub <- gr[ seqnames(gr) %in% chunk]
				  lapply( bam.files, get( bamsignals.function ), gr=gr.sub, mapqual=mapqual, shift=shift, paired.end=paired.end, paired.end.midpoint=paired.end.midpoint, verbose=verbose)
				}, mc.cores=procs)
	list("s"=unlist(lapply(x, "[[",1)),
		 "r" =unlist(lapply(x, "[[",2)) 
		 )
}

#EOF
