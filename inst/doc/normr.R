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
citation("normr")

## ------------------------------------------------------------------------
#Loading required package
library(normr)

inputBamfile <- system.file("extdata", "K562_Input.bam", package="normr")
k4me3Bamfile <- system.file("extdata", "K562_H3K4me3.bam", package="normr")
k36me3Bamfile <- system.file("extdata", "K562_H3K36me3.bam", package="normr")

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
#Enrichment Calling for H3K4me3 and H3K36me3
k4me3Fit <- enrichR(treatment = k4me3Bamfile, control = inputBamfile,
                    genome = genome, verbose = F)
k36me3Fit <- enrichR(treatment = k36me3Bamfile, control = inputBamfile,
                     genome = genome, verbose = F)

## ------------------------------------------------------------------------
k4me3Fit

## ------------------------------------------------------------------------
k36me3Fit

## ------------------------------------------------------------------------
summary(k4me3Fit)

## ------------------------------------------------------------------------
summary(k36me3Fit)

## ------------------------------------------------------------------------
#background normalized enrichment
k4me3Enr <- getEnrichment(k4me3Fit)

#restrict to regions with non-zero counts
idx <- which(rowSums(do.call(cbind, getCounts(k4me3Fit))) != 0)
summary(k4me3Enr[idx])

## ----fig.width=10, fig.height=10-----------------------------------------
x <- k4me3Enr[idx]
y <- getEnrichment(k36me3Fit)[idx]
d.x <- density(x); d.y <- density(y)
limits <- range(x,y)
layout( matrix( c(0,2,2,1,3,3,1,3,3),ncol=3) )
plot(d.x$x, d.x$y, xlim=limits, type='l',
     main="H3K36me3 normalized Enrichment", xlab="", ylab="Density")
abline(v=0, lty=3, lwd=2, col=4)
plot(d.y$y, d.y$x, ylim=limits, xlim=rev(range(d.y$y)), type='l',
     main="H3K4me3 normalized Enrichment", xlab="Density", ylab="")
abline(h=0, lty=3, lwd=2, col=4)
color <- rep("grey10", length(idx))
plot(x, y, xlim=limits, ylim=limits, pch=20, xlab="", ylab="",
     col=adjustcolor(color, alpha.f=.2))
abline(0, 1, lty=2, lwd=3, col=2)
abline(v=0, lty=3, lwd=2, col=4)
abline(h=0, lty=3, lwd=2, col=4)

## ----fig.width=10, fig.height=10-----------------------------------------
#integer vector with <NA> set to non-significant regions
k4me3Classes <- getClasses(k4me3Fit, fdr = 0.05)
k36me3Classes <- getClasses(k36me3Fit, fdr = 0.05)

#Color scatter plot based on enrichment
color[!is.na(k4me3Classes[idx])] <- "#2C9500"
color[!is.na(k36me3Classes[idx])] <- "#990099"
color[!is.na(k4me3Classes+k36me3Classes)[idx]] <- "#971621"
plot(x, y, xlim=limits, ylim=limits, pch=20,
     col=adjustcolor(color, alpha.f=.5), xlab="H3K4me3 normalized Enrichment",
     ylab="H3K36me3 normalized Enrichment")
legend("topright", pch=20, col=unique(color), cex=.6, bg="white",
  legend=c("Background", "H3K36me3 enriched", "H3K4me3 enriched",
           "H3K4me3/K36me3 enriched")
  )

## ------------------------------------------------------------------------
k4me3Ranges <- getRanges(k4me3Fit)[!is.na(k4me3Classes)]
#Alternatively you can extract ranges without storing the class vector
k4me3Ranges <- getRanges(k4me3Fit, fdr = 0.05)

#as expected we get 380 regions
length(k4me3Ranges)

## ----warning=FALSE-------------------------------------------------------
#example gene annotation for representative region (chr1:22500000-25000000)
genes <- read.delim(file = system.file("extdata", "genes.bed", package="normr"),
                    header = F,
                    stringsAsFactors = F)
library(GenomicRanges, verbose = F)
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
  fisher.test(m, alternative="greater")
}

#Overlap of H3K4me3 with genes
overlapOdds(k4me3Ranges, genes)

#Overlap of H3K4me3 with promoters
promoters <- promoters(genes, upstream = 2000, downstream = 2000)
overlapOdds(k4me3Ranges, promoters)

## ------------------------------------------------------------------------
#Overlap of H3K36me3 with H3K4me3
k36me3Ranges <- getRanges(k36me3Fit, fdr = 0.05)
overlapOdds(k36me3Ranges, k4me3Ranges)

#Overlap of H3K36me3 with H3K4me3 at promoter regions
overlapOdds(k36me3Ranges[countOverlaps(k36me3Ranges, promoters) != 0],
            k4me3Ranges[countOverlaps(k4me3Ranges, promoters) != 0])

#Overlap of H3K36me3 with H3K4me3 in genes
overlapOdds(k36me3Ranges[countOverlaps(k36me3Ranges, genes) != 0],
            k4me3Ranges[countOverlaps(k4me3Ranges, genes) != 0])

## ------------------------------------------------------------------------
#Overlap of H3K36me3 in genes
overlapOdds(k36me3Ranges, genes)

#Overlap of H3K36me3 with promoters
overlapOdds(k36me3Ranges, promoters(genes, 1500, 1500))

## ----warning=FALSE-------------------------------------------------------
#export coordinates of significantly (FDR <= 0.05) enriched regions
exportR(k4me3Fit, filename = "k4me3Fit.bed", type = "bed", fdr = 0.05)
exportR(k36me3Fit, filename = "k36me3Fit.bed", type = "bed", fdr = 0.05)

#export background-normalized enrichment
exportR(k4me3Fit, filename = "k4me3Fit.bw", type = "bigWig")
exportR(k36me3Fit, filename = "k36me3Fit.bw", type = "bigWig")

## ----warning=FALSE-------------------------------------------------------
#We could use read counts from above NormRFit objects
k4k36Dif <- diffR(treatment = getCounts(k4me3Fit)$treatment,
                  control   = getCounts(k36me3Fit)$treatment,
                  genome    = getRanges(k4me3Fit),
                  verbose   = F)
#<or> (unnecessarily) count again
#k4k36Dif <- diffR(treatment = k4me3Bamfile, control = k36me3Bamfile,
#                  genome = genome, verbose = F)

#summary statistics
summary(k4k36Dif)

## ----warning=FALSE-------------------------------------------------------
exportR(k4k36Dif, filename = "k4k36Dif.bed", type = "bed", fdr = 0.05)
exportR(k4k36Dif, filename = "k4k36Dif.bw", type = "bigWig")

## ----warning=FALSE-------------------------------------------------------
k4me3Regimes <- regimeR(treatment = getCounts(k4me3Fit)$treatment,
                         control   = getCounts(k4me3Fit)$control,
                         genome    = getRanges(k4me3Fit),
                         models    = 3,
                         verbose   = F)
summary(k4me3Regimes)

## ----warning=FALSE-------------------------------------------------------
k36me3Regimes <- regimeR(treatment = getCounts(k36me3Fit)$treatment,
                         control   = getCounts(k36me3Fit)$control,
                         genome    = getRanges(k36me3Fit),
                         models    = 3,
                         verbose   = F)
summary(k36me3Regimes)

## ----warning=FALSE-------------------------------------------------------
exportR(k4me3Regimes, filename = "k4me3Regimes.bed", type = "bed", fdr = 0.05)
exportR(k36me3Regimes, filename = "k36me3Regimes.bed", type = "bed", fdr = 0.05)

## ----eval=FALSE----------------------------------------------------------
#  #Single End:
#  # Count in 500bp bins.
#  # Consider only reads with Mapping Quality >= 20.
#  # Filter reads for marked duplicates (e.g. with picard mark-duplicates)
#  # Shift the counting position for a read 100 bp downstream.
#  countConfigSE <- countConfigSingleEnd(binsize = 500, mapq = 20,
#                                        filteredFlag = 1024, shift = 100)
#  
#  #Paired End:
#  # Count in 500bp bins.
#  # Consider only reads with Mapping Quality >= 30.
#  # Count the midpoint of the aligned fragment instead of 5' ends.
#  # Consider only reads corresponding to fragments with size from 100 to 300bp
#  countConfigPE <- countConfigPairdEnd(binsize = 500, mapq = 30, midpoint = T,
#                                       tlenfilter = c(100, 300))
#  
#  #Plug in the counting configuration into normR, e.g. in enrichR()
#  fit <- enrichR(treatment   = k4me3Bamfile,
#                 control     = inputBamfile,
#                 genome      = genome,
#                 countConfig = countConfigPE)

## ----warning=FALSE-------------------------------------------------------
#regions have identical size?
all(width(promoters) == 4000)

## ----eval=FALSE----------------------------------------------------------
#  #count in predefined regions
#  library(bamsignals)
#  k4Counts <- bamCount(bampath = k4me3Bamfile, gr = promoters, verbose = F)
#  inputCounts <- bamCount(bampath = inputBamfile, gr = promoters, verbose = F)
#  
#  #Fit only on promoters
#  promotersFit <- enrichR(treatment = k4Counts, control = inputCounts,
#                          genome = promoters, verbose = F)

## ----eval=FALSE----------------------------------------------------------
#  cnvs <- diffR(treatment   = treatmentInputBamfile,
#                control     = controlInputBamfile,
#                genome      = genome,
#                countConfig = countConfigSingleEnd(binsize = 2.5e4))
#  
#  #export the CNV calls
#  exportR(cnvs, "CNVs.bed")
#  
#  #Filter previous ChIP-seq difference calls for CNVs
#  ov <- countOverlaps(getRanges(diffFit, fdr = .05), getRanges(cnvs, fdr = .05))
#  idx <- which(ov == 0)
#  cnvCleanedGR <- getRanges(diffFit, fdr = .05)[idx]

