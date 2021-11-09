context("GetBarcodes Dataframe")

library(automata)
library(testthat)

barcodes_CHOL_df <- getBarcodes("TCGA-CHOL")
barcodes_BRCA_df <- getBarcodes("TCGA-BRCA", microarray = TRUE)

test_that("getBarcodes admits lower case",{
  dummy <- getBarcodes("tCgA-cHoL")
  expect_equal(class(dummy), "data.frame")
})

test_that("CHOL-RNA dataframe is generated", {
  expect_equal(class(barcodes_CHOL_df), "data.frame")
})

test_that("CHOL-RNA dataframe has five columns", {
  expect_equal(ncol(barcodes_CHOL_df), 5)
})

test_that("CHOL-RNA dataframe is not empty", {
  expect(nrow(barcodes_CHOL_df) > 0, "The dataframe is empty" )
})

test_that("BRCA-MICROARRAY dataframe is generated", {
  expect_equal(class(barcodes_BRCA_df), "data.frame")
})

test_that("CHOL dataframe has five columns", {
  expect_equal(ncol(barcodes_BRCA_df), 5)
})

test_that("CHOL dataframe is not empty", {
  expect(nrow(barcodes_BRCA_df) > 0, "The dataframe is empty" )
})

