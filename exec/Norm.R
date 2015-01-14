#!/usr/local/package/bin/Rscript --vanilla

# Argument parsing
library(optparse, quietly=T)

description <- "
Norm.R was tested with R version 3.1. 
Necessary libraries might be installed via CRAN with invoking the R command:

	install.packages( c(\"optparse\", \"GenomicRanges\", \"parallel\", \"rtracklayer\", \"devtools\") );
	install_github( \"lamortenera/bamsignals\" );
	"
epilogue <- "Johannes Helmuth <helmuth@molgen.mpg.de>"
option_list <- list(
 		    make_option(c("-t", "--treatment"), 
				help="Path to treatment bamfile for normalization. Please ensure that an index file with .bai suffix is found in the same directory [needed]"),
		    make_option(c("-c", "--control"), 
				help="Path to control bamfile for normalization. Please ensure that an index file with .bai suffix is found in the same directory [needed]"),
		    make_option(c("-g", "--genome"), 
				help="Genome chromosome sizes file, e.g. chromSize table from UCSC."),
		    make_option(c("-b", "--bin.size"), type="integer", default=1000,
				help="Size of bins for enrichment calculation in bp. If no bin.size is set, a best fitting size (under BIC) for bins (multiple of 50 bp) will be assigned. [default=1000]"),
		    make_option(c("-m", "--models"), type="integer", default="2",
				help="How many binomial model components are desired [default=2]"),
		    make_option(c("-o", "--output.dir"), 
				help="Path to output directory [default=$(treatment.basename)Signal_$(control.basename)Background_$(bin.size)bp]"),
		    make_option(c("-j", "--jobs"), type="integer", default="1",
				help="Number of jobs used for normalization. Ideally this should be set to the number of chromosomes, if enough processors are available [default=1]"),
		    make_option(c("-f", "--mapqual"), type="integer", default="30",
				help="Consider only reads with mapping quality >= mapqual. Note, in paired end data only the first mate in proper pair is counted [default=30]"),
		    make_option(c("-s", "--shift"), type="integer", default="0",
				help="Integer specifying shift of counted tag relative to 5'-fragment end [default=0]"),
		    make_option(c("-p", "--paired.end"), default="FALSE", action="store_true",
				help="Consider data to be paired end. Fragment middle points will be counted. [default=FALSE]"),
	  	    make_option(c("-q", "--quiet"), default="FALSE", action="store_false",
				help="Flag indicating if quiet output [unset by default].")
	 	    )
args <- parse_args(OptionParser(option_list = option_list, description=description, epilogue=epilogue))
if ( is.null( args$treatment ) ) {
    cat( "ERROR: No treatment file supplied. Exiting...\n" )
    q( "no", status=1 )
}
if ( is.null( args$control ) ) {
    cat( "ERROR: No control file supplied. Exiting...\n" )
    q( "no", status=1 )
}


# create output directory
default.prefix <- F
if ( is.null(prefix) ) {
 	prefix <- paste( basename(args$treatment), "Signal_", basename(args$control), "Background", sep="") 
	default.prefix <- T
}
dir.create(prefix, recursive=T, showWarnings=F)


# diffr normalization
if (!args$quiet)
	cat( "[Starting Run @", paste(Sys.time()),"] ", prefix , "\n" )
try( library( diffr ), silent=T)
result <- normalize( treatment = args$treatment, 
					 control   = args$control, 
					 genome    = args$genome, 
					 bin.size  = args$bin.size, 
					 models    = args$models,
					 eps       = 0.001, 
					 procs     = args$jobs, 
					 mapqual   = args$mapqual,
					 shift     = args$shift, 
					 paired.end= args$paired.end, 
					 paired.end.midpoint= args$paired.end,
					 verbose   = !args$quiet)


# write result
if ( !args$quiet )
	cat( "[Started @", paste(Sys.time()),"] Saving bigwig file...\n" )
gr <- bin.genome( args$genome, args$bin.size )
try( library(rtracklayer), silent=T )
gr$score <- 1 - result[[1]][,1] # gives probability of NOT being background. Useful for more than 2 models
gr$score[ is.na(gr$score) ] <- 0
export.bw( gr, paste(prefix, "/", basename(args$treatment), "_normalized.bw", sep="") )

if ( !args$quiet )
	cat( "[Started @", paste(Sys.time()),"] Saving normalized value bed file...\n" )

#CONTINUE HERE

mat <- cbind( as.data.frame(gr)[,c(1,2,3)], -log10(result[[1]][,3]), result[[1]][,1:2],  counts[[1]], counts[[2]], result[[1]][,4])
write.table( x     = mat, 
			 file  = paste(prefix, "/", basename(args$treatment), "_normalized.bed", sep=""), 
			 quote = F, 
			 sep   = "\t", 
			 row.names=F, 
			 col.names=F)


if ( !args$quiet)
	cat( "[Finishing Run @", paste(Sys.time()),"] ", prefix , "\n" ) 

#EOF
