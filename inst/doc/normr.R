## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----eval=FALSE----------------------------------------------------------
#  e <- enrichR(treatment = "ChIP.bam",
#               control   = "Control.bam",
#               genome    = "hg19")

## ----eval=FALSE----------------------------------------------------------
#  de <- diffR(treatment = "ChIP1.bam",
#              control   = "ChIP2.bam",
#              genome    = "hg19")

## ----eval=FALSE----------------------------------------------------------
#  re <- regimeR(treatment = "ChIP.bam",
#                control   = "Control.bam",
#                genome    = "hg19",
#                models    = k)

## ----eval=FALSE----------------------------------------------------------
#  #export enriched regions with FDR<=10% for downstream analysis
#  exportR(obj      = e,
#          filename = "enriched.bed",
#          type     = "bed",
#          fdr      = 0.1)
#  #or
#  #write normalized differential enrichment to bigWig for genome browser display
#  exportR(obj      = de,
#          filename = "diffEnrichment.bw",
#          type     = "bigWig")

## ------------------------------------------------------------------------
#Loading required package
library(normr)

inputBamfile <- system.file("extdata", "K562_Input.bam", package="normr")
k4me3Bamfile <- system.file("extdata", "K562_H3K4me3.bam", package="normr")
k27me3Bamfile <- system.file("extdata", "K562_H3K27me3.bam", package="normr")

## ----warning=FALSE-------------------------------------------------------
#Fetch chromosome information
genome <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC("hg19")

#Filter out irregular chromosomes and delete unnecessary columns
idx <- which(!genome$circular & genome$SequenceRole=="assembled-molecule")
genome <- genome[idx,1:2]
genome

#Toy data has only "chr1"
genome <- genome[genome[,1] == "chr1",]
genome

## ----warning=FALSE-------------------------------------------------------
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
#plot(k4me3Fit)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
#get path to executable
inputBamfile <- system.file("extdata", "K562_Input.bam", package="normr")


