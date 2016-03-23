context( "normr-methods test" )

source("utils.R")

getArtificalData <- function(n) {

}

test_that("Function arguments are checked correctly", {
  #for (func.str in c("enrichR", "diffR", "regimeR")) {
  for (func.str in c("enrichR")) {#TODO for now
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

  expect_error(exportR(new("NormRFit"), "dmp", "xxx"), info="Invalid format")
  expect_error(exportR(new("NormRFit"), "dmp", "bed"), info="Invalid format")
  expect_error(exportR(new("NormRFit"), "dmp", "bigWig"), info="Invalid format")
})

test_that("Fitting with enrichR() works correctly", {
  genome <- GRanges("chr1", IRanges(22500000, 25000000))#chr1:1-249250621
  inputfile <- system.file("extdata", 
                          "K562_Input.bam", 
                          package="normr")
  chipfiles <- system.file("extdata", 
                          paste0("K562_", c("H3K27me3", "H3K4me3"), ".bam"), 
                          package="normr")

  #check files first
  sapply(c(inputfile, chipfiles), function(b) expect_true(file.exists(b)))
  sapply(c(inputfile, chipfiles), function(b) expect_true(file.exists(paste0(b, ".bai"))))

  #test for multiple binsizes
  for (binsize in c(150, 500, 1000)) {
    for (chipfile in chipfiles) {
      label <- paste0("enrichR{",
        "region=", seqnames(genome),":",start(genome),"-",end(genome),
        ", binsize=", binsize, ",chipfile=", chipfile,
        ",inputfile=", inputfile, "}")

    Rfit <- RenrichR(suppressWarnings(bamProfile(chipfile,genome,binsize))[1],
                     suppressWarnings(bamProfile(inputfile,genome,binsize))[1])
    Cfit <- enrichR(Rfit$treatment, Rfit$control,
                    unlist(tile(genome, width=binsize)), verbose=F)

    expect_equal(label=paste0(label, " - Treatment counts"),
                 Rfit$treatment, getCounts(Cfit)[["treatment"]])
    expect_equal(label=paste0(label, " - Control counts"),
                 Rfit$control, getCounts(Cfit)[["control"]])
    }
  }
  #expect_equal(label=paste0("bamCount{", 
  #                          paste("shift", shift, "mapq", mapq, "ss", ss, "pe", pe, 
  #                                "tlen.filter", paste0(tlen.filter, collapse=","), 
  #                                sep="="),
  #                          "}"),
  #             countR(reads, regions, ss=ss, shift=shift, paired.end=pe, mapqual=mapq,
  #                    tlen.filter=tlen.filter), 
  #             bamCount(bampath, regions, ss=ss, shift=shift, paired.end=pe, mapqual=mapq, 
  #                      tlen.filter=tlen.filter, verbose=FALSE))
  #load_all()
  #require(GenomicRanges)
  #genome <- GRanges("chr1", IRanges(22500000, 25000000))#chr1:1-249250621
  #require(bamsignals)
  #binsize=1000
  #datasets <- paste0("K562_", c("H3K27me3", "H3K4me3", "Input"), ".bam")
  #bamfiles <- system.file("extdata", datasets, package="normr")
  #count_ctrl <- bamProfile(bamfiles[3],genome,binsize)[1]
  #count_treat <- bamProfile(bamfiles[1],genome,binsize)[1]
  #load_all()
  #gr <- unlist(tile(genome, width=binsize))
  #Cfit <- enrichR(count_ctrl, count_treat, gr)

})

test_that("Fitting with diffR() works correctly", {
})

test_that("Recomputation of P works", {
})

test_that("Exporting of results works", {
})
