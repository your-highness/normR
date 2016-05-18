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

## ------------------------------------------------------------------------
summary(k27me3Fit)

## ------------------------------------------------------------------------
#integer vector with <NA> set to non-significant regions
k4me3Classes <- getClasses(k4me3Fit, fdr = 0.05)
k4me3Ranges <- getRanges(k4me3Fit)[!is.na(k4me3Classes)]
#Alternatively you can extract ranges without storing the class vector
k4me3Ranges <- getRanges(k4me3Fit, fdr = 0.05)

#as expected we get 577 regions
length(k4me3Ranges)

## ----warning=FALSE-------------------------------------------------------
#example gene annotation for representative region (chr1:22500000-25000000)
genes <- read.delim(file = system.file("extdata", "genes.bed", package="normr"),
                    header = F,
                    stringsAsFactors = F)
library(GenomicRanges)
genes <- GRanges(seqnames = genes[, 1],
                 ranges = IRanges(start = genes[, 2], end = genes[, 3]),
                 strand = genes[, 6],
                 ENSTID = genes[, 4])
genes <- unique(genes)

#Fisher-test provides significance of overlap
#(total specifies number of bins in representative region)
overlapOdds <- function(query, subject, total = 10000) {
  subject <- reduce(subject, ignore.strand = T)
  ov1 <- countOverlaps(query, subject)
  m <- matrix(c(sum(ov1 != 0), sum(ov1 == 0),
              ceiling(sum(width(subject))/width(query)[1]-sum(ov1 != 0)), 0),
              ncol = 2)
  m[2,2] <- total - sum(m)
  fisher.test(m)
}

#Overlap of H3K4me3-enriched with genes
overlapOdds(k4me3Ranges, genes)
#Overlap of H3K4me3-enriched with promoters
overlapOdds(k4me3Ranges, promoters(genes))

## ------------------------------------------------------------------------
#Overlap of H3K4me3 with H3K27me3
k27me3Ranges <- getRanges(k27me3Fit, fdr = 0.05)
overlapOdds(k4me3Ranges, k27me3Ranges)

## ------------------------------------------------------------------------
#Overlap of H3K27me3-enriched with genes
overlapOdds(k27me3Ranges, genes)
#Overlap of H3K27me3-enriched with promoters
overlapOdds(k27me3Ranges, promoters(genes))

## ----warning=FALSE-------------------------------------------------------
#export coordinates of significantly (FDR <= 0.05) enriched regions
exportR(k4me3Fit, filename = "k4me3Fit.bed", type = "bed", fdr = 0.05)
exportR(k27me3Fit, filename = "k27me3Fit.bed", type = "bed", fdr = 0.05)

#export background-normalized enrichment
exportR(k4me3Fit, filename = "k4me3Fit.bw", type = "bigWig")
exportR(k27me3Fit, filename = "k27me3Fit.bw", type = "bigWig")

## ----warning=FALSE-------------------------------------------------------
#We could use read counts from above NormRFit objects
k4k27Dif <- diffR(treatment = getCounts(k4me3Fit)$treatment,
                  control   = getCounts(k27me3Fit)$treatment,
                  genome    = getRanges(k4me3Fit),
                  verbose   = F)
#or
#We could count again which is unnecessary
#k4k27Dif <- diffR(treatment = k4me3Bamfile, control = k27me3Bamfile, genome = genome, verbose = F)

#summary statistics
summary(k4k27Dif)

## ----warning=FALSE-------------------------------------------------------
k27me3Regimes <- regimeR(treatment = getCounts(k27me3Fit)$treatment,
                         control   = getCounts(k27me3Fit)$control,
                         genome    = getRanges(k27me3Fit),
                         models    = 3,
                         verbose   = F)

## ----eval=FALSE----------------------------------------------------------
#  #TODO

