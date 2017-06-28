source("utils.R")

#some stub data for object generation
getData <- function(n, type) {
  if (type == "enrichR") {
    k <- 2L
    B <- 1L
  } else if (type == "diffR") {
    k <- 3L
    B <- 2L
  } else if (type == "regimeR") {
    k <- 3L
    B <- 1L
  }
  r <- rbinom (n, 10, .2)
  s <- r + rbinom(n, 40, .01)
  map <- map2unique(cbind(r,s))
  theta <- seq(.2,.9, length.out=k)
  mixtures <- seq(.9,.2, length.out=k)
  mixtures <- mixtures/sum(mixtures)
  lnpost <- sapply(theta, dbinom, x=s, size=r+s)
  lnpost <- log(lnpost)
  enr <- (log(s + .1) - log(r +.05)) / log(.05/.1)
  pvals <- pbinom(s, r+s, theta[1], F) +
           dbinom(s, r+s, theta[1])
  filteredT <- which((r+s) > 3)
  qvals <- rep(NA, length(r))
  qvals[filteredT] <- pvals[filteredT] * pvals[filteredT]#stub
  classes <- as.integer(rep(NA,length(qvals)))
  classes[qvals < 5e-2] <- 1L
  list("k"=k,"B"=B,"r"=r,"s"=s,"map"=map,"theta"=theta,"mixtures"=mixtures,
    "lnpost"=lnpost,"lnenr"=enr, "lnpvals"=log(pvals), "filteredT"=filteredT,
    "thresholdT"=3L, "lnqvals"= log(qvals), "classes"=classes)
}

runs <- function(expr){
  res <- try(force(expr), TRUE)
  msg1 <- "code did not generate an error"
  msg2 <- "code generated an error"
  expectation("error", msg2, msg1)
}

expect_runs <- function(expr){
  expect_that(expr, runs, label=testthat:::find_expr("expr"))
}

context("NormRFit-class")
test_that("Test NormRFit construction, validity and accessors", {
  n <- 1000
  gr <- GRanges("chr1", IRanges(seq(1,n*100-1,100), width=100))
  seqinfo(gr) <- Seqinfo("chr1", n*100)
  for(type in c("enrichR", "diffR", "regimeR")) {
    d <- getData(n, type)
    o <- new("NormRFit", type=type, n=as.integer(n),
             ranges=gr, k=d$k, B=d$B, map=d$map$map,
             counts=list(d$map$values[1,],d$map$values[2,]),
             amount=as.integer(table(d$map$map)),
             names=c("treatment", "control"),
             thetastar=sum(d$s)/sum(d$r+d$s), theta=d$theta,
             mixtures=d$mixtures, lnL=seq(-n,-1,10)+rnorm(n,0,.5), eps=.001,
             lnposteriors=applyMap(d$lnpost, d$map),
             lnenrichment=applyMap(d$lnenr, d$map),
             lnpvals=applyMap(d$lnpvals, d$map),
             filteredT=sort(unique(d$map$map[d$filteredT])),
             thresholdT=d$thresholdT,
             lnqvals=applyMap(d$lnqvals, d$map),
             classes=applyMap(d$classes, d$map))

    expect_that(o, is_a("NormRFit"))

    ###
    #output methods
    expect_error(exportR(o, "x.x"))
    #bed
    expect_silent(exportR(o, tempfile(), type = "bed"))
    expect_silent(exportR(o, tempfile(), type = "bed", fdr = NA))
    expect_error(exportR(o, tempfile(), type = "bed", color = rep("black", 4)))
    expect_silent(exportR(o, tempfile(fileext = ".bed")))
    #bigWig
    expect_silent(exportR(o, tempfile(), type = "bigWig"))
    expect_silent(exportR(o, tempfile(fileext=".bigWig")))
    expect_silent(exportR(o, tempfile(fileext=".bw")))
    #bedGraph
    expect_silent(exportR(o, tempfile(), type = "bedGraph"))
    expect_silent(exportR(o, tempfile(fileext=".bedGraph")))
    expect_silent(exportR(o, tempfile(fileext=".bg")))

    ###
    #plotting - not implemented yet
    expect_error(plot(o))

    ###
    #accessors to NormRFit internal map
    #counts
    expect_that(str(getCounts(o)), prints_text("List of 2"))
    expect_equal(getCounts(o)$control, d$r)
    expect_equal(getCounts(o)$treatment, d$s)
    #ranges
    gr <- GRanges(gr, "component"=d$classes)
    expect_equal(getRanges(o), gr)
    expect_equal(getRanges(o), gr)
    expect_equal(getRanges(o, k = 1), gr[which(gr$component == 1)])
    expect_equal(getRanges(o, fdr = .05), gr[which(d$lnqvals <= log(.05))])
    expect_equal(getRanges(o, fdr = .05, k = 1),
                 gr[which(d$lnqvals <= log(.05))])
    expect_error(getRanges(o, k = -1))
    expect_error(getRanges(o, k = 4))
    expect_error(getRanges(o, fdr = .05, k = 0))
    #posteriors
    expect_equal(getPosteriors(o), exp(d$lnpost))
    #enrichment
    expect_equal(getEnrichment(o), d$lnenr)
    if (type == "enrichR") {
      enr <- RgetEnrichment(exp(d$lnpost), d$r, d$s, d$theta)
      expect_equal(getEnrichment(o, B = 1), enr)
      expect_equal(getEnrichment(o, B = 1, F = 2), enr)
    } else if (type == "diffR") {
      enr <- RgetEnrichmentDiff(exp(d$lnpost), d$r, d$s, d$theta)
      expect_equal(getEnrichment(o, B = 2), enr)
      expect_equal(getEnrichment(o, B = 2, F = 1), enr)
    } else if (type == "regimeR") {
      enr <- RgetEnrichment(exp(d$lnpost), d$r, d$s, d$theta)
      expect_equal(getEnrichment(o, B = 1, F = 2), enr)
      enr <- RgetEnrichment(exp(d$lnpost), d$r, d$s, d$theta, fgIdx=3)
      expect_equal(getEnrichment(o, B = 1, F = 3), enr)
    }
    expect_error(getEnrichment(o, B = -1))
    expect_error(getEnrichment(o, B = 9))
    expect_error(getEnrichment(o, F = -1))
    expect_error(getEnrichment(o, F = 9))
    #pvalues
    pvals <- exp(d$lnpvals)
    expect_equal(getPvalues(o), pvals)
    idx <- which(rowSums(do.call(cbind, getCounts(o))) > d$thresholdT)
    expect_equal(getPvalues(o, filtered = TRUE), pvals[idx])
    #qvalues
    expect_equal(getQvalues(o), exp(d$lnqvals))
    expect_equal(sum(is.na(getQvalues(o))), n-length(idx))
    #classes
    classes <- d$classes
    expect_equal(getClasses(o), classes)
    idx <- which(d$lnqvals <= exp(0.05))
    expect_equal(getClasses(o, fdr = .05), classes)

    ###
    #printing functions
    #length
    expect_equal(length(o), n)
    #print
    expect_runs(capture.output(print(o)))
    o2 <- o; o2@theta <- as.numeric(NULL)
    expect_runs(capture.output(print(o2)))
    #show
    expect_runs(capture.output(show(o)))
    #summary
    expect_runs(capture.output(summary(o)))
    expect_runs(capture.output(summary(o2)))

  }
})

