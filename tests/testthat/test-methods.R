context( "normr-methods test" )

source("utils.R")



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

test_that("Fitting with counting from bam files works correctly", {
})

test_that("Fitting with counts files works correctly", {
})

test_that("Recomputation of P works", {
})

test_that("Exporting of results works", {
})
