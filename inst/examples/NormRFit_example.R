require(GenomicRanges)

#Create a toy instance of type 'enrichR'
fit <- new("NormRFit",
           type="enrichR", n=10L,
           ranges=GRanges("chr1", IRanges(seq(1,100,10), width=10)),
           k=2L, B=1L, map=rep(1:5,2), counts=list(1:5, 1:5),
           amount=rep(2L,5), names=c("chip", "input"), thetastar=.35,
           theta=c(.15,.55), mixtures=c(.9,.1), lnL=seq(-50,-1,10), eps=.001,
           lnposteriors=log(matrix(runif(10), ncol=2)),
           lnenrichment=log(runif(5,0,.2)), lnpvals=log(runif(5)),
           filteredT=2:5, thresholdT=1L, lnqvals=log(runif(5,0,.2)),
           classes=sample(1:2,5,TRUE))

#print some statistics on fits
fit
summary(fit)

## Not run:
#write significant regions to bed
#exportR(fit, filename = "enrich.bed", fdr = 0.1)
#write normalized enrichment to bigWig
#exportR(fit, filename = "enrich.bw")
## End(**Not run**)

###AccessorMethods
#get original counts
getCounts(fit)
#get genomic coordinates for significant ranges as a GenomicRanges instance
getRanges(fit, fdr = .1)
getPosteriors(fit)
getEnrichment(fit)
getPvalues(fit)
getQvalues(fit)
getClasses(fit)
