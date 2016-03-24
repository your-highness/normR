inputfile <- system.file("extdata", "K562_Input.bam", package="normr")
chipfiles <- system.file("extdata", "K562_H3K27me3.bam"), package="normr")
gr <- GRanges("chr1", IRanges(22500000, 25000000))#chr1:1-249250621
fit <- enrichR(chipfile, inputfile, genome.df, countConfigPairedEnd(), verbose=F)
summary(fit)
