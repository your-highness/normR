## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ------------------------------------------------------------------------
#e <- enrichR("ChIP.bam", "Control.bam", "hg19")

## ------------------------------------------------------------------------
#de <- diffR("ChIP1.bam", "ChIP2.bam", "hg19")

## ------------------------------------------------------------------------
#re <- regimeR("ChIP.bam","Control.bam", "hg19", k)

## ------------------------------------------------------------------------
#Loading required package
library(normr)

inputBamfile <- system.file("extdata", "K562_Input.bam", package="normr")
k4me3Bamfile <- system.file("extdata", "K562_H3K4me3.bam", package="normr")
k27me3Bamfile <- system.file("extdata", "K562_H3K27me3.bam", package="normr")

## ------------------------------------------------------------------------
#Fetch chromosome information
#genome <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC("hg19")
#
##Filter out irregular chromosomes and delete unnecessary columns
#idx <- which(!genome$circular & genome$SequenceRole=="assembled-molecule")
#genome <- genome[idx,1:2]
#
##Toy data has only "chr1"
#genome <- genome[genome[,1] == "chr1",]

genome = data.frame("seqname"="chr1", "length"=25e6)
genome

## ------------------------------------------------------------------------
#Enrichment Calling for H3K4me3 and H3K27me3
k4me3Fit <- enrichR(treatment = k4me3Bamfile, control = inputBamfile, genome = genome, verbose = F)
k27me3Fit <- enrichR(treatment = k27me3Bamfile, control = inputBamfile, genome = genome, verbose = F)

## ------------------------------------------------------------------------
k4me3Fit

## ------------------------------------------------------------------------
k27me3Fit

## ------------------------------------------------------------------------
summary(k4me3Fit)

## ---- fig.show='hold', fig.cap="Density of enrichment", fig.width=5, fig.height=3----
par(mfrow=c(1,2), bty="n")
hist(getEnrichment(k4me3Fit), breaks=20, main="H3K4me3 Enrichment")
hist(getEnrichment(k4me3Fit), breaks=20, main="H3K4me3 Enrichment")

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
#get path to executable
inputBamfile <- system.file("extdata", "K562_Input.bam", package="normr")


