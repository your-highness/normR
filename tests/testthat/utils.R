# Some helper methods for test that routines

getEnrichmentR <- function(fit, bg.idx=1) {
  p <- fit$posterior[,bg.idx]
  count_ctrl <- fit$control
  count_treat <- fit$treatment

  pseu_ctrl <- sum(p * count_ctrl) / sum(p)
  pseu_treat <- sum(p * count_treat) / sum(p)

  return(log((count_treat + pseu_treat)/(count_ctrl + pseu_ctrl))/
   log(pseu_ctrl / pseu_treat))
}

#PVALUE COMPUTATION
T.filter <- function(fit, thresh=.0001, bg.idx=1) {
  marg = 0 
  r = 0
  s = 0
  run=T
  border = 0
  pval.comp <- function(fit, r, s) {
    pbinom(s, r+s, fit$theta[bg.idx], lower.tail=F) + dbinom(s, r+s,
                                                             fit$theta[bg.idx])
  }
  while (run) {
    p <- pval.comp(fit, marg, 0)
    message("r=", marg, "; s=", 0, "; pval=", p)
    if (p <= thresh) {
      border = marg
      break
    }
    if ( (marg - 1) > 0 ) {
      for (i in (marg-1):1) {
        p <- pval.comp(fit, i, marg-i)
        message("r=", i, "; s=", marg-i, "; pval=", p)
        if (p <= thresh) {
          border = marg-i
          run=F
          break
        }
        if (marg-i != i) {
          p <- pval.comp(fit, marg-i, i)
          message("r=", marg-i, "; s=", i, "; pval=", p)
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

  #pi_0 becomes very small here
  #return(which((fit$treatment + fit$control) >= marg & fit$control > 0 & fit$treatment > border))
  return(which((fit$treatment + fit$control) >= marg))
}
#l is a list of c(treatment, control, binsize)
processFiles <- function(l, pe=paired.end, models=2, shift=0, nthreads=procs, gen=genome, bg.idx=1) {
  require(diffr)
  nor <- lapply(names(l), function(b) {
                if (l[[b]]["paired"] == T) {
                  pe = "midpoint"
                } else {
                  pe = "ignore"
                  shift=100
                }
                n <- diffR(treatment = l[[b]]["treatment"],
                           control   = l[[b]]["control"], 
                           genome    = gen,
                           bin.size  = binsize[b],
                           models    = models,
                           procs     = nthreads,
                           mapqual   = mapq,
                           shift     = shift,
                           paired.end= pe, #"ignore" (SE), "midpoint" (pe)
                           verbose   = T)


                n$pval <- sapply(1:models, function(i) {
                                 p <- pbinom(n$treatment, (n$control + n$treatment), n$theta[i], lower.tail=F) + dbinom(n$treatment, (n$control + n$treatment), n$theta[i])
                                 p[which(p > 1)] <- 1
                                 p[which(p < 0)] <- 0
                                 (p)
                           })

                #Apply T filter
                n$idx <- T.filter(n, .001, bg.idx=bg.idx)

                #Compute FDR
                require(qvalue)
                fdr <- qvalue(n$pval[n$idx, bg.idx])
                n$fdr <- fdr$qvalues
                n$significant <- n$idx[which(n$fdr <= fdr.thresh[b])]
                if (models > 2) {
                  n$regime <- apply(n$posterior[n$significant,-bg.idx], 1, which.max)
                }
                n    
   })
  names(nor) <- names(l)
  return (nor) 
}


#ENRICHMENT CALC
getEnrichment <- function(fit, bg.idx=1) {
  p <- fit$posterior[,bg.idx]
  count_ctrl <- fit$control
  count_treat <- fit$treatment

  pseu_ctrl <- sum(p * count_ctrl) / sum(p)
  pseu_treat <- sum(p * count_treat) / sum(p)

  log( (count_treat + pseu_treat) / (count_ctrl + pseu_ctrl) ) / log(pseu_ctrl / pseu_treat)
}

#DIFFR CALL
T.filter.diff <- function(fit, thresh=.0001, bg.idx=1) {
  marg = 0  
  r = 0
  s = 0
  run=T
  border = 0
  pval.comp <- function(fit, r, s) {
    if ((r+s) == 0)
      1
    else
      binom.test(s, r+s, fit$theta[bg.idx], alternative="two.sided")$p.value
  }
  while (run) {
    p <- pval.comp(fit, marg, 0)
    message("r=", marg, "; s=", 0, "; pval=", p)
    if (p <= thresh) {
      border = marg
      break 
    }
    if ( (marg - 1) > 0 ) {
      for (i  in (marg-1):1) {
        p <-  pval.comp(fit, i, marg-i)
        message("r=", i, "; s=", marg-i, "; pval=", p)
        if (p <= thresh) {
          border = marg-i
          run=F
          break
        }
        if (marg-i != i) {
          p <- pval.comp(fit, marg-i, i)
          message("r=", marg-i, "; s=", i, "; pval=", p)
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

  #pi_0 becomes very small here
  #return(which((fit$treatment + fit$control) >= marg & fit$control > 0 & fit$treatment > border))
  return(which((fit$treatment + fit$control) >= marg))
}
recomputeSignificance <- function(n, bg.idx=2) {
  binom.test2 <- function(x, n, ...) {
    if ((x+n) == 0)
      1
    else 
      binom.test(x, n, ...)$p.value
  }
  n$pval <- mcmapply(binom.test2, x=n$treatment, n=(n$control + n$treatment), MoreArgs=list("p"=n$theta[bg.idx], "alternative"="two.sided"), mc.cores=procs)

  #Apply T filter
  n$idx <- T.filter.diff(n, .001, bg.idx=2)

  require(qvalue)
  fdr <- qvalue(n$pval[n$idx])
  n$fdr <- fdr$qvalue
  n$significant <- n$idx[which(n$fdr <= fdr.thresh["H3K27me3"])]
  n$regime <- apply(n$posterior[n$significant,-bg.idx], 1, function(r) {
                    reg <- which(r > .9)
                    if (length(reg) == 0)
                      NA
                    else
                      reg
   })
  n
}


#WRITE RESULTS
writekModelsBed <- function(filename, nor, trackname="", trackdescr="", post=post.thresh, k=3, cols=NULL, nthreads=procs) {
  if (is.null(cols))
    cols <- rainbow(k-1)

  #put names for list
  peaks <- mclapply(2:k, function(i) {
                    idx <- nor$significant[which(nor$regime == (i-1))]
                    res <- nor$ranges[idx]
                    start(res) <- start(res)-1 #adjust for bed display
                    res$score <- as.integer(nor$posterior[idx,i] * 1000)
                    res$name <- paste0("k", i, "_score:", res$score) 
                    #res$col <- paste(col2rgb(cols[i-1]), collapse=",")
                    #Shading of colors dependent on posteriors (skipped for now)
                    colRamp <- apply(col2rgb(colorRampPalette(c("grey",cols[i-1]))(6)), 2, paste, collapse=",")
                    res$col <- colRamp[ceiling(nor$posterior[idx,i]/.2)+1]
                    #colRamp <- apply(col2rgb(colorRampPalette(c("white",cols[i-1]))(6)), 2, paste, collapse=",")
                    #res$col <- colRamp[((res$score - post * 1000)/1000/((1-post)/5) + 1)] #don't ask!
                    (res)
   }, mc.cores=nthreads)

  peaks <- do.call(c, peaks)
  peaks <- sort(peaks)

  cat(paste0('track name="', trackname, '" description="', trackdescr, '" visibility=2 itemRgb="On" \n'),
      file=filename)
  write.table(file=filename,
              x=as.data.frame(peaks)[,c(1,2,3,7,6,5,2,3,8)],
              sep="\t",
              col.names=F,
              row.names=F,
              quote=F,
              append=T)
}

#write bedgraph with acceptable precision (smaller filesize)
writeBedGraph <- function(fit, filename, bg.idx=1, log=F) {
  require(rtracklayer)
  gr <- fit$ranges
  fit.e <- getEnrichment(fit, bg.idx)
  gr$score <- as.numeric(format(fit.e, 1,1))
  if (!log)
    gr$score <- exp(gr$score)
  #restrict to regions where Input or ChIP have read counts -> reduces file size immensly
  idx <- fit$treatment == 0 | fit$control == 0
  message("Removing ", length(which(idx)), " regions with zero control or treatment counts.")
  gr <- gr[!idx]
  export.bedGraph(gr, filename)
  #TODO export a useful TrackLine
  #gtl <- new("GraphTrackLine", ...)
  #export.bedGraph(gr, filename, trackLine=gtl)
}
#write bigwig with acceptable precision (smaller filesize)
writeBigWig <- function(fit, filename, bg.idx=1) {
  require(rtracklayer) 
  gr <- fit$ranges
  fit.e <- getEnrichment(fit, bg.idx)
  gr$score <- exp(as.numeric(format(fit.e, 1,1)))
  #restrict to regions where Input or ChIP have read counts -> reduces file size immensly
  idx <- fit$treatment == 0 | fit$control == 0
  message("Removing ", length(which(idx)), " regions with zero control or treatment counts.")
  gr <- gr[!idx]
  export.bw(gr, filename)
}

writeEachDiffBed <- function(filename.prefix, nor, trackname="", trackdescr="", cols=NULL, nthreads=procs, bg.idx=2) {
  #dmp <- lapply(1:length(nor), function(i) {
  dmp <- mclapply(1:length(nor), function(i) {
                  idx <- nor[[i]]$significant#index and significance was previously computed in processFiles()
                  res <- nor[[i]]$ranges[idx]
                  res$score <- NA
                  res$col <- NA
                  for (j in 1:(NCOL(nor[[i]]$posterior)-1)) {
                    idx.2 <-which(nor[[i]]$regime == j) 
                    res$score[idx.2] <- as.integer(nor[[i]]$posterior[idx,-bg.idx][idx.2,j] * 1000)
                    colRamp <- apply(col2rgb(colorRampPalette(c("grey",cols[j]))(6)), 2, paste, collapse=",")
                    res$col[idx.2]   <- colRamp[ceiling(nor[[i]]$posterior[idx,-bg.idx][idx.2,j]/.2)+1]
                  }
                  res <- res[which(!is.na(res$score))]#remove the ones where no regime could be assigned by posterior threshold
                  res$name <- paste0(names(nor)[i], "_Cond",nor[[i]]$regime[!is.na(nor[[i]]$regime)] ,"_score:", res$score) 
                  res <- sort(res)

                  #writing
                  name <- names(nor)[i]
                  filename <- paste0(filename.prefix, "_", name, ".bed")
                  cat(paste0('track name=', paste0(trackname,"_", name), ' description="', paste(trackdescr, "enriched for", name), '" visibility=2 itemRgb="On"\n'),    file=filename, append=F)
                  write.table(file=filename,
                              x=as.data.frame(res)[,c(1,2,3,8,6,5,2,3,7)],
                              sep="\t",
                              col.names=F,
                              row.names=F,
                              quote=F,
                              append=T)
                  #                })
              }, mc.cores=nthreads)
}
