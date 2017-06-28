#normR - testing utilitiy functions

###
# R MAP2UNIQUE PAIRS IMPLEMENTATION
###
#rle for a matrix
rle.matrix <- function(x) {
        n <- dim(x)[1]
        if (n == 0L) return( list(values = numeric(), map = integer()) )
        y <- !sapply(2:n, function(i) { all(x[i,] == x[(i-1),]) } )
        i <- c(which(y | is.na(y)), n)
        structure(list(lengths = diff(c(0L, i)), values = x[i,]))
}
# MAP2UNIQUE function
#returns a list with two elements:
#1: values, sorted unique values of the counts vector
#2: map, reference to the values vector such that:
# for vectors: all(values[map] == counts)
# for matrices: all(values[map,] == counts)
map2unique <- function(counts){
    if (is.matrix(counts)) {
        o <- do.call(order, lapply(1:NCOL(counts), function(i) counts[, i]))
        uval <- rle.matrix( counts[o,] )
        n <- dim(counts)[1]
        m <- dim(uval$values)[1]
    } else if (is.vector(counts)) {
        o <- order(counts)
        uval <- rle(counts[o])
        n <- length(counts)
        m <- length(uval$values)
    }
    values <- uval$values
    uval$values <- 1:m
    map <- integer(n)
    map[o] <- inverse.rle( uval )
    return(list(values=t(values), map=map))
}
#maps a vector v to unique value space
applyMap <- function(v, map) {
  idx <- which(!duplicated(map$map))
  idx <- idx[order(map$map[idx])]
  if (class(v) == "matrix") v[idx,]
  else v[idx]
}

###
# R normR implementation
###
#the workhorse implemented in R
RnormR <- function(s, r, nmodels=2, eps=1e-5){
  idx <- which((s + r) > 0)
  N <- length(idx); n <- s + r
  mixtures <- runif(nmodels)
  thetastar <- sum(s[idx]) / sum(n[idx])
  theta <- rep(thetastar, nmodels) - runif(nmodels, 0, (thetastar - eps)) # q*
  theta <- sort(theta)
  #helper functions
  logRowSum <- function(x) {
    apply(x,1,function(r) {
          m <- max(r)
          tmp <- sum(exp(r-m))
          return(m + log(tmp))
    })
  }
  ## EM
  runs <- 0; lnL <- -Inf; not.converged <- T
  while (runs < 30 | not.converged) { #Ensure burn in
    lnmixtures <- log(mixtures)
    ## Expectation:
    likelihood <- sapply(theta, function(p) log(p) * s[idx] + log(1 - p) *
                         r[idx])
    likelihood <- sapply(1:nmodels, function(i) likelihood[, i] + lnmixtures[i])
    lnZ <- logRowSum(likelihood)
    posteriors <- exp(likelihood - lnZ)
    mixtures <- colSums(posteriors, na.rm=TRUE)
    mixtures <- mixtures / sum(mixtures)
    theta <- colSums( posteriors * s[idx], na.rm=T) /
               colSums( posteriors * n[idx], na.rm=T)
    o <- order(theta)
    theta <- theta[o]
    mixtures <- mixtures[o]

    ## Convergence
    lnL.new <- sum(lnZ, na.rm=T)
    if (runs > 30 & abs(lnL.new - lnL) < eps)
      not.converged <- F
    lnL <- lnL.new
    runs <- runs + 1
  }
  ##posteriorserio and Pvalue computation for whole data set
  likelihood <- sapply(theta, function(p) log(p) * s + log(1 - p) * r)
  likelihood <- sapply(1:nmodels, function(i) likelihood[,i]+log(mixtures[i]))
  lnZ <- logRowSum(likelihood)
  posteriors <- exp(likelihood - lnZ)
  list(control=r, treatment=s, idx=idx, thetastar=thetastar, theta=theta,
       mixtures=mixtures, lnL=log(sum(exp(lnZ))), eps=eps,
       posteriors=posteriors)
}

###
# R enrichR methods
###
RgetP <- function(s, r, p) {
  return(pbinom(s, r+s, p, lower.tail=F) + dbinom(s, r+s, p))
}
Rtfilter <- function(fit, thresh=5e-2, bgIdx=1) {
    marg = 0
    r = 0
    s = 0
    run=T
    border = 0
    while (run) {
        p <- RgetP(0, marg, fit$theta[bgIdx])
        if (p <= thresh) {
            border = marg
            break
        }
        if ( (marg - 1) > 0 ) {
            for (i in (marg-1):1) {
                p <- RgetP(marg-i, i, fit$theta[bgIdx])
                if (p <= thresh) {
                    border = marg-i
                    run=F
                    break
                }
                if (marg-i != i) {
                    p <- RgetP(i, marg-i, fit$theta[bgIdx])
                    if (p <= thresh) {
                        border = i
                        run=F
                        break
                    }
                }
            }
        }
        marg = marg + 1
    }
    return(which((fit$treatment + fit$control) >= marg))
}
RgetEnrichment <- function(post, r, s, theta, bgIdx=1, fgIdx=2) {
  p <- post[,bgIdx]
  pseu_r <- sum(p * r) / sum(p)
  pseu_s <- sum(p * s) / sum(p)
  foldchange <- log((s+pseu_s)/(r+pseu_r))
  regularization <- log(pseu_r / pseu_s)
  standardization <-
    log(theta[fgIdx]/(1-theta[fgIdx])*(1-theta[bgIdx])/theta[bgIdx])
  return((foldchange + regularization)/standardization)
}
RenrichR <- function(s, r, eps=1e-5, bgIdx=1) {
  fit <- RnormR(s,r, eps = eps)

  #get enrichment
  fit$lnenrichment <- RgetEnrichment(fit$posteriors,r,s,fit$theta)

  #calculate pvalues
  fit$pvals <- RgetP(s,r,fit$theta[bgIdx])
  fit$pvals[fit$pvals > 1] <- 1
  fit$pvals[fit$pvals < 0] <- 0

  #Apply Rtfilter filter
  fit$filteredT <- Rtfilter(fit, bgIdx=bgIdx)

  #Q values
  fit$qvals <- rep(NA, length(r))
  p <- fit$pvals
  fit$qvals[fit$filteredT] <- 
    qvalue::qvalue(fit$pvals[fit$filteredT], 
      lambda=seq(min(p), min(0.99, max(p)), .05))$qvalues

  #classes
  fit$classes <- as.integer(rep(NA, length(r)))
  fit$classes <- apply(fit$posteriors,1,which.max)
  fit$classes[fit$classes == 1] <- NA
  fit$classes <- fit$classes - 1

  return(fit)
}

###
# R diffR methods
###
RgetPDiff <- function(s, r, p) {
  return(sapply(1:length(s), function(i) {
    if ((s[i]+r[i]) == 0) 1
    else stats::binom.test(s[i], r[i]+s[i], p, alternative="two.sided")$p.value
  }))
}
RtfilterDiff <- function(fit, thresh=1e-2, bgIdx=2) {
  marg = 0
  r = 0
  s = 0
  run=T
  border = 0
  while (run) {
    p <- RgetPDiff(0, marg, fit$theta[bgIdx])
    if (p <= thresh) {
      border = marg
      break
    }
    if ( (marg - 1) > 0 ) {
      for (i  in (marg-1):1) {
        p <-  RgetPDiff(marg-i, i, fit$theta[bgIdx])
        if (p <= thresh) {
          border = marg-i
          run=F
          break
        }
        if (marg-i != i) {
          p <- RgetPDiff(i, marg-i, fit$theta[bgIdx])
          if (p <= thresh) {
            border = i
            run=F
            break
          }
        }
      }
    }
    marg = marg + 1
  }
  return(which((fit$treatment + fit$control) >= marg))
}
RgetEnrichmentDiff <- function(post, r, s, theta) {
  p <- post[,2]
  pseu_r <- sum(p * r) / sum(p)
  pseu_s <- sum(p * s) / sum(p)
  foldchange <- log((s+pseu_s)/(r+pseu_r))
  regularization <- log(pseu_r / pseu_s)

  #Standardized foldchange dependent on algebraic sign
  foldchange <- foldchange + regularization
  standardizationT <-
    log(theta[3]/(1-theta[3])*(1-theta[2])/theta[2])
  foldchange[foldchange > 0] <- foldchange[foldchange > 0]/standardizationT
  standardizationC <-
    -log(theta[1]/(1-theta[1])*(1-theta[2])/theta[2])
  foldchange[foldchange < 0] <- foldchange[foldchange < 0]/standardizationC
  return(foldchange)
}
RdiffR <- function(s, r, eps=1e-5) {
  fit <- RnormR(s,r,3,eps)

  #get enrichment
  fit$lnenrichment <- RgetEnrichmentDiff(fit$posteriors,r,s,fit$theta)

  #calculate pvalues
  fit$pvals <- RgetPDiff(s,r,fit$theta[2])
  fit$pvals[fit$pvals > 1] <- 1
  fit$pvals[fit$pvals < 0] <- 0

  #Apply Rtfilter filter
  fit2 <- RnormR(r,s,3)
  fit$filteredT <- intersect(RtfilterDiff(fit, eps, 2),
                             RtfilterDiff(fit2, eps, 2))

  #Q values
  fit$qvals <- rep(NA, length(r))
  p <- fit$pvals
  fit$qvals[fit$filteredT] <- 
    qvalue::qvalue(fit$pvals[fit$filteredT], 
      lambda=seq(min(p), min(0.99, max(p)), .05))$qvalues

  #classes
  fit$classes <- as.integer(rep(NA, length(r)))
  fit$classes <- apply(fit$posteriors[,c(1,3,2)],1,which.max)
  fit$classes[fit$classes == 3] <- NA

  return(fit)
}

###
# R regimeR methods
###
RregimeR <- function(s, r, nmodels, eps=1e-5, bgIdx=1) {
  fit <- RnormR(s,r,nmodels, eps=eps)

  #get enrichment
  fit$lnenrichment <- RgetEnrichment(fit$posteriors,r,s,fit$theta,1,nmodels)

  #calculate pvalues
  fit$pvals <- RgetP(s,r,fit$theta[bgIdx])
  fit$pvals[fit$pvals > 1] <- 1
  fit$pvals[fit$pvals < 0] <- 0

  #Apply Rtfilter filter
  fit$filteredT <- Rtfilter(fit, bgIdx=bgIdx)

  #Q values
  fit$qvals <- rep(NA, length(r))
  p <- fit$pvals
  fit$qvals[fit$filteredT] <- 
    qvalue::qvalue(fit$pvals[fit$filteredT], 
      lambda=seq(min(p), min(0.99, max(p)), .05))$qvalues

  #classes
  fit$classes <- as.integer(rep(NA, length(r)))
  fit$classes <- apply(fit$posteriors[,1:nmodels],1,which.max)
  fit$classes[fit$classes == 1] <- NA
  fit$classes <- fit$classes - 1

  return(fit)
}

###
# R implementation of getting classes by FDR filter
###
RgetClasses <- function(Rfit, fdr=.1, bgIdx=1) {
  idx <- which(Rfit$qvals <= fdr)
  guys <- Rfit$classes[idx]
  na <- which(is.na(guys))
  if (length(na) > 0) {
    if (length(Rfit$theta) == 2) {
      guys[na] <- 1
    } else if (length(na) == 1) {
      guys[na] <- which.max(Rfit$posteriors[idx,][na,-bgIdx])
    } else{
      guys[na] <- apply(Rfit$posteriors[idx,][na,-bgIdx], 1, which.max)
    }
  }
  classes <- as.integer(rep(NA_integer_, length(Rfit$classes)))
  classes[idx] <- guys
  classes
}
