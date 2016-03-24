context("NormRFit class and methods testing")

source("utils.R")

#some stub data for object generation
getData <- function(n) {
  r <- rbinom(n, 10, .2)
  s <- r + rbinom(n, 40, .01)
  map <- map2unique(cbind(r,s))
  theta <- c(.2,.4)
  mixtures <- c(.8,.2)
  lnpost <- cbind(0,dbinom(s, r+s, theta[2]))
  lnpost[,1] <- 1-lnpost[,2]
  lnpost <- log(lnpost)
  enr <- (log(s + .1) - log(r +.05)) / log(.05/.1)
  pvals <- pbinom(s, r+s, theta[1], F) +
           dbinom(s, r+s, theta[1])
  filteredT <- which((r+s) > 3)
  qvals <- rep(NA, length(r))
  qvals[filteredT] <- pvals[filteredT] * pvals[filteredT]#stub
  classes <- as.integer(rep(NA,length(qvals)))
  classes[qvals < 5e-2] <- 1L
  list("r"=r,"s"=s,"map"=map,"theta"=theta,"mixtures"=mixtures,"lnpost"=lnpost,
   "lnenr"=enr, "lnpvals"=log(pvals), "filteredT"=filteredT,
   "lnqvals"=log(qvals), "classes"=classes)
}

runs <- function(expr){
  res <- try(force(expr), TRUE)
  msg1 <- "code did not generate an error"
  msg2 <- "code generated an error"
  expectation(!inherits(res, "try-error"), msg2, msg1)
}

expect_runs <- function(expr){
  expect_that(expr, runs, label=testthat:::find_expr("expr"))
}

test_that("Test NormRFit construction, validity and accessors", {
  n <- 1000
  for(type in c("enrichR", "diffR", "regimeR")) {
    d <- getData(n)
    o <- new("NormRFit", type=type, n=as.integer(n),
             ranges=GRanges("chr1", IRanges(seq(0,n*100-1,100), width=100)),
             k=2L, B=1L, map=d$map$map,
             counts=list(d$map$values[1,],d$map$values[2,]),
             thetastar=sum(d$s)/sum(d$r+d$s), theta=d$theta,
             mixtures=d$mixtures, lnL=seq(-n,-1,10)+rnorm(n,0,.5), eps=.001,
             lnposteriors=applyMap(d$lnpost, d$map),
             lnenrichment=applyMap(d$lnenr, d$map),
             lnpvals=applyMap(d$lnpvals, d$map),
             filteredT=sort(unique(d$map$map[d$filteredT])),
             lnqvals=applyMap(d$lnqvals, d$map),
             classes=applyMap(d$classes, d$map))

    expect_that(o, is_a("NormRFit"))

    #length
    expect_equal(length(o), n)

    #object printing
    expect_runs(capture.output(print(o)))
    expect_runs(capture.output(show(o)))
    expect_runs(capture.output(summary(o)))

    #plotting - not implemented yet TODO
    expect_error(plot(o))

    #accessors to NormRFit internal map
    expect_that(str(getCounts(o)), prints_text("List of 2"))
    expect_equal(getCounts(o)$control, d$r)
    expect_equal(getCounts(o)$treatment, d$s)
    expect_equal(getEnrichment(o), d$lnenr)
    expect_equal(getPosteriors(o), exp(d$lnpost))
    expect_equal(getPvalues(o), exp(d$lnpvals))
    expect_equal(getQvalues(o), exp(d$lnqvals))
    expect_equal(getClasses(o), d$classes)
    expect_equal(sum(is.na(getQvalues(o))), n-length(d$filteredT))

    #output methods - not implemented yet
    expect_error(writeEnrichment(o))
    expect_error(writeBed(o))
  }
})

