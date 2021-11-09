context("Download TCGA data")

library(automata)
library(testthat)

test_that("if automata_id not given the function does not crash",{
  downloadClinical("TCGA-CHOL", "TCGA-3X-AAV9")
  dummy <- read_csv("TCGA_results/TCGA-CHOL/data/")
  expect_equal(class(dummy), "data.frame")
})

