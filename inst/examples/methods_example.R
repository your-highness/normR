require(GenomicRanges)

### enrichR(): Calling Enrichment over Input
#load some example bamfiles
input <- system.file("extdata", "K562_Input.bam", package="normr")
chipK4 <- system.file("extdata", "K562_H3K4me3.bam", package="normr")
#region to count in (example files contain information only in this region)
gr <- GRanges("chr1", IRanges(seq(22500001, 25000000, 1000), width = 1000))
#configure your counting strategy (see BamCountConfig-class)
countConfiguration <- countConfigSingleEnd(binsize = 1000,
                                           mapq = 30, shift = 100)
#invoke enrichR to call enrichment
enrich <- enrichR(treatment = chipK4, control = input,
                  genome = gr,  countConfig = countConfiguration,
                  iterations = 10, procs = 1, verbose = TRUE)
#inspect the fit
enrich
summary(enrich)

## Not run:
#write significant regions to bed
#exportR(enrich, filename = "enrich.bed", fdr = 0.01)
#write normalized enrichment to bigWig
#exportR(enrich, filename = "enrich.bw")
## End(**Not run**)

### diffR(): Calling differences between two conditions
chipK36 <- system.file("extdata", "K562_H3K36me3.bam", package="normr")
diff <- diffR(treatment = chipK36, control = chipK4,
              genome = gr,  countConfig = countConfiguration,
              iterations = 10, procs = 1, verbose = TRUE)
summary(diff)

### regimeR(): Identification of broad and peak enrichment
regime <- regimeR(treatment = chipK36, control = input, models = 3,
                  genome = gr,  countConfig = countConfiguration,
                  iterations = 10, procs = 1, verbose = TRUE)
summary(regime)
