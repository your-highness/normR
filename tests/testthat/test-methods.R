source("utils.R")

context("methods arguments are checked correctly")
test_that("Function arguments are checked correctly", {
  for (func.str in c("enrichR", "diffR", "regimeR")) {
    func <- get(func.str)

    #input files
    expect_error(func("dmp.bam", "dmp2.bam"), info="Non existing files")
    if (file.create(c("dmp.bam", "dmp2.bam"), showWarnings=F)[1]) {
      expect_error(func("dmp.bam", "dmp2.bam", data.frame()),
        info="No index")
      file.remove(c("dmp.bam", "dmp2.bam"))
    }
    expect_error(func("dmp.bam", c(1,1,1,1)),
      info="file and vector input not allowed")

    #input vectors
    expect_error(func(c(0.1,0.1,0.1), c(0.1,0.1,0.1), NULL),
      info="Numeric type")
    expect_error(func(c(1,1,1), c(1,1,1,1), NULL),
      info="Different vector length")

    if (func.str %in% c("diffR", "regimeR")) {
      expect_error(func(c(1,1,1), c(1,1,1), NULL, k=0),
        info="models=0")
      expect_error(func(c(1,1,1), c(1,1,1), NULL, k=-1),
        info="models<0")
      }
  }

  #check data handling methods
  expect_error(recomputeP(new("NormRFit"), 0), info="B=0")
  expect_error(recomputeP(new("NormRFit"), -1), info="B<0")

  expect_error(exportR(new("NormRFit"), "dmp", type = "xxx"),
               info="Invalid format")
  expect_error(exportR(new("NormRFit"), "dmp", type = "bed"),
               info="Invalid format")
  expect_error(exportR(new("NormRFit"), "dmp", type = "bigWig"),
               info="Invalid format")
})

#comparing the fit computed by normR with the R implementation
testFit <- function(Rfit, Cfit, label, tolerance=1e-4, bgIdx=1) {
  #data vectors
  expect_equal(label=paste0(label, " - length"), length(Rfit$treatment),
               length(Cfit), tolerance=tolerance)
  expect_equal(label=paste0(label, " - counts treatment"), Rfit$treatment,
               getCounts(Cfit)[["treatment"]], tolerance=tolerance)
  expect_equal(label=paste0(label, " - counts control"), Rfit$control,
               getCounts(Cfit)[["control"]], tolerance=tolerance)
  #model fit
  expect_equal(label=paste0(label, " - thetastar"), Rfit$thetastar,
               Cfit@thetastar, tolerance=tolerance)
  expect_equal(label=paste0(label, " - theta"), Rfit$theta,
               Cfit@theta, tolerance=tolerance)
  expect_equal(label=paste0(label, " - mixtures"), Rfit$mixtures,
               Cfit@mixtures, tolerance=tolerance)
  expect_equal(label=paste0(label, " - posteriors"), Rfit$posteriors,
               getPosteriors(Cfit), tolerance=tolerance)
  #model inferred statistics
  expect_equal(label=paste0(label, " - lnenrichment"), Rfit$lnenrichment,
               getEnrichment(Cfit), tolerance=tolerance)
  expect_equal(label=paste0(label, " - Pvalues"), Rfit$pvals,
               getPvalues(Cfit), tolerance=tolerance)
  expect_equal(label=paste0(label, " - Tfiltered"), Rfit$filteredT,
               which(Cfit@map %in% Cfit@filteredT), tolerance=tolerance)
  expect_equal(label=paste0(label, " - Qvalues"), Rfit$qvals,
               getQvalues(Cfit), tolerance=1e-3)
  expect_equal(label=paste0(label, " - Classes"), Rfit$classes,
               getClasses(Cfit), tolerance=tolerance)
  expect_equal(label=paste0(label, " - Number with Qvalues<=0.1"),
               sum(Rfit$qvals <= .1), sum(getQvalues(Cfit) <= .1))
  if (length(Rfit$theta) > 2) {#test correctness of maximum in regime posteriors
    expect_equal(label=paste0(label, " - Maximum A Posteriori"),
                 apply(Rfit$posteriors[,-bgIdx],1,which.max),
                 apply(getPosteriors(Cfit)[,-bgIdx], 1, which.max))
  }
  expect_equal(label=paste0(label, " - Classes with Qvalues<=0.1"),
               RgetClasses(Rfit, fdr=.1, bgIdx=bgIdx), getClasses(Cfit, fdr=.1))
}

context("enrichR() gives correct results")
test_that("enrichR() works correctly", {
  gr <- GRanges("chr1", IRanges(22500001, 25000000))#chr1:22500000-25000000
  inputfile <- system.file("extdata", "K562_Input.bam", package="normr")
  chipfiles <- system.file("extdata", paste0("K562_", c("H3K36me3", "H3K4me3"),
                           ".bam"), package="normr")

  #check files first
  sapply(c(inputfile, chipfiles), function(b)
         expect_true(file.exists(b)))
  sapply(c(inputfile, chipfiles), function(b)
         expect_true(file.exists(paste0(b, ".bai"))))

  #input:integer,interger,GRanges
  for (chipfile in chipfiles) {
    for (binsize in c(500, 1000)) {
      counts_ctrl <- suppressWarnings(bamProfile(inputfile,gr,binsize))[1]
      counts_treat <- suppressWarnings(bamProfile(chipfile,gr,binsize))[1]
      testFit(
        RenrichR(counts_treat, counts_ctrl),
        enrichR(counts_treat, counts_ctrl, unlist(tile(gr, width=binsize)),
                verbose=F),
        paste0("enrichR-integer,integer,GRanges{", "region=",
               seqnames(gr), ":", start(gr),"-",end(gr),
               ", binsize=", binsize, ",chipfile=", chipfile, ",inputfile=",
               inputfile, "}")
      )
    }
  }

  #input:character,character,data.frame
  genome.df <- data.frame(c("chr1", "chr2"), c(25e6, 1e3))
  binsize <- 1e3
  for (chipfile in chipfiles) {
    gr2 <- GRanges(genome.df[,1], IRanges(1, genome.df[,2]))
    testFit(
      RenrichR(unlist(as.list(suppressWarnings(
                 bamProfile(chipfile,gr2,500,20,0,F,"ignore",verbose=F)))),
               unlist(as.list(suppressWarnings(
                 bamProfile(inputfile,gr2,500,20,0,F,"ignore",verbose=F))))),
      enrichR(chipfile, inputfile, genome.df,
              countConfigSingleEnd(500,20,-1,0), verbose=F),
      paste0("enrichR-character,character,data.frame-SingleEnd{region=",
             seqnames(gr), ":", start(gr),"-",end(gr),
             ", binsize=500,chipfile=", chipfile, ",inputfile=",
              inputfile, "}")
    )
    fit <- enrichR(chipfile, inputfile, genome.df,
                   countConfigPairedEnd(500,20,-1,0,T,c(0,300)), verbose=F)
    testFit(
      RenrichR(unlist(as.list(suppressWarnings(
                 bamProfile(chipfile,gr2,500,20,0,F,"midpoint",c(0,300),1024,
                            verbose=F)))),
               unlist(as.list(suppressWarnings(
                 bamProfile(inputfile,gr2,500,20,0,F,"midpoint",c(0,300),1024,
                            verbose=F))))),
      fit,
      paste0("enrichR-character,character,data.frame-PairedEnd{region=",
             seqnames(gr), ":", start(gr),"-",end(gr),
             ", binsize=500,chipfile=", chipfile, ",inputfile=",
             inputfile, "}")
    )

    #exporting
    expect_silent(exportR(fit, tempfile(), type = "bed"))
    expect_silent(exportR(fit, tempfile(), type = "bedGraph"))
    expect_silent(exportR(fit, tempfile(), type = "bigWig"))
  }
})

context("diffR() gives correct results")
test_that("diffR() works correctly", {
  gr <- GRanges("chr1", IRanges(22500001, 25000000))#chr1:1-249250621
  chipfiles <- system.file("extdata", paste0("K562_", c("H3K36me3", "H3K4me3"),
                           ".bam"), package="normr")

  #input:integer,interger,GRanges
  for (binsize in c(500, 1000)) {
    counts_ctrl <- suppressWarnings(bamProfile(chipfiles[1],gr,binsize))[1]
    counts_treat <- suppressWarnings(bamProfile(chipfiles[2],gr,binsize))[1]
    testFit(
      RdiffR(counts_treat, counts_ctrl, eps=5e-2),
      diffR(counts_treat, counts_ctrl, unlist(tile(gr, width=binsize)),
            verbose=F, eps=5e-2, minP=5e-2),
      paste0("diffR-integer,integer,GRanges{", "region=",
             seqnames(gr), ":", start(gr),"-",end(gr),
             ", binsize=", binsize, ",chipfile=", chipfiles[2], ",inputfile=",
             chipfiles[1], "}"),
      bgIdx=2,
      tolerance=5e-3
    )
  }

  #input:character,character,data.frame
  genome.df <- data.frame("chr1", 25000000)
  gr <- GRanges("chr1", IRanges(1, 25000000))
  #testFit(
  #  RdiffR(suppressWarnings(
  #           bamProfile( chipfiles[2],gr,500,20,0,F,"ignore",verbose=F))[1],
  #         suppressWarnings(
  #           bamProfile(chipfiles[1],gr,500,20,0,F,"ignore",verbose=F))[1],
  #         eps=5e-2),
  #  diffR(chipfiles[2], chipfiles[1], genome.df,
  #        countConfigSingleEnd(500,20,-1,0), verbose=F, eps=5e-2, minP=5e-2),
  #  paste0("diffR-character,character,data.frame-SingleEnd{region=",
  #         seqnames(gr), ":", start(gr),"-",end(gr),
  #         ", binsize=500,chipfile=", chipfiles[2], ",inputfile=",
  #         chipfiles[1], "}"),
  #  bgIdx=2,
  #  tolerance=1e-1
  #)
  fit <- diffR(chipfiles[2], chipfiles[1], genome.df,
               countConfigPairedEnd(500,20,-1,0,T,c(0,300)), verbose=F,
               eps=5e-2, minP=5e-2)
  testFit(
    RdiffR(suppressWarnings(
             bamProfile(chipfiles[2],gr,500,20,0,F,"midpoint",c(0,300),1024,
                        verbose=F))[1],
           suppressWarnings(
             bamProfile(chipfiles[1],gr,500,20,0,F,"midpoint",c(0,300),1024,
                        verbose=F))[1],
           eps=5e-2),
    fit,
    paste0("diffR-character,character,data.frame-PairedEnd{region=",
           seqnames(gr), ":", start(gr),"-",end(gr),
           ", binsize=500,chipfile=", chipfiles[2], ",inputfile=",
            chipfiles[1], "}"),
    bgIdx=2,
    tolerance=5e-2
  )

  #exporting
  expect_silent(exportR(fit, tempfile(), type = "bed"))
  expect_silent(exportR(fit, tempfile(), type = "bedGraph"))
  expect_silent(exportR(fit, tempfile(), type = "bigWig"))
})

context("regimeR() gives correct results")
test_that("regimeR() works correctly", {
  gr <- GRanges("chr1", IRanges(22500001, 25000000))#chr1:22500000-25000000
  inputfile <- system.file("extdata", "K562_Input.bam", package="normr")
  chipfiles <- system.file("extdata", paste0("K562_", c("H3K36me3", "H3K4me3"),
                           ".bam"), package="normr")

  #input:integer,interger,GRanges
  for (chipfile in chipfiles) {
    for (binsize in c(500, 1000)) {
      k <- 3L
      counts_ctrl <- suppressWarnings(bamProfile(inputfile,gr,binsize))[1]
      counts_treat <- suppressWarnings(bamProfile(chipfile,gr,binsize))[1]
      testFit(
        RregimeR(counts_treat, counts_ctrl, k),
        regimeR(counts_treat, counts_ctrl, unlist(tile(gr, width=binsize)),
                k, verbose=F),
        paste0("regimeR-integer,integer,GRanges{", "region=",
               seqnames(gr), ":", start(gr),"-",end(gr),
               ", binsize=", binsize, ",chipfile=", chipfile, ",inputfile=",
               inputfile, ",models=", k, "}"),
        tolerance=5e-3
      )
    }
  }

  #input:character,character,data.frame
  genome.df <- data.frame("chr1", 25000000)
  for (chipfile in chipfiles) {
    k <- 3L
    gr <- GRanges("chr1", IRanges(1, 25000000))
    testFit(
      RregimeR(suppressWarnings(
                 bamProfile(chipfile,gr,500,20,0,F,"ignore",verbose=F))[1],
               suppressWarnings(
                 bamProfile(inputfile,gr,500,20,0,F,"ignore",verbose=F))[1],
               k),
      regimeR(chipfile, inputfile, genome.df, k,
              countConfigSingleEnd(500,20,-1,0), verbose=F),
      paste0("regimeR-character,character,data.frame-SingleEnd{region=",
             seqnames(gr), ":", start(gr),"-",end(gr),
             ", binsize=500,chipfile=", chipfile, ",inputfile=",
             inputfile, ",models=", k, "}"),
      tolerance=5e-3
    )
    Cfit <- regimeR(chipfile, inputfile, genome.df, k,
                    countConfigPairedEnd(500,20,-1,0,T,c(0,300)), verbose=F)
    testFit(
      RregimeR(suppressWarnings(
                 bamProfile(chipfile,gr,500,20,0,F,"midpoint",c(0,300),1024,
                            verbose=F))[1],
               suppressWarnings(
                 bamProfile(inputfile,gr,500,20,0,F,"midpoint",c(0,300),1024,
                            verbose=F))[1],
               k),
      Cfit,
      paste0("regimeR-character,character,data.frame-PairedEnd{region=",
             seqnames(gr), ":", start(gr),"-",end(gr),
             ", binsize=500,chipfile=", chipfile, ",inputfile=",
             inputfile, ",models=", k, "}"),
      tolerance=5e-3
    )

    #exporting
    expect_silent(exportR(Cfit, tempfile(), type = "bed"))
    expect_silent(exportR(Cfit, tempfile(), type = "bedGraph"))
    expect_silent(exportR(Cfit, tempfile(), type = "bigWig"))
  }
})
