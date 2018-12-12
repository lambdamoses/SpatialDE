context("SpatialDE on matrix")
library(testthat)

test_that("RunSpatialDE io", {
  data("mouseob_small")
  data("locations")
  expect_error(RunSpatialDE(locations[,2:3], mouseob_small))
  expect_error(RunSpatialDE(locations, mouseob_small, row_sample = TRUE))
  expect_error(RunSpatialDE(locations[,2:3], t(mouseob_small), row_sample = TRUE))
  results <- RunSpatialDE(locations[,2:3], mouseob_small, row_sample = TRUE)
  expect_equal(nrow(results), ncol(mouseob_small))
  expect_equal(ncol(results), 18)
})
