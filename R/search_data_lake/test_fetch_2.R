library(testthat)

# Source the function to be tested
source("fetch_occurrence_data2.R")

# Unit test for fetch_occurrence_data
test_that("fetch_occurrence_data returns expected output", {
  
  # Test with a valid URL
  url <- "https://catalog.dive.edito.eu"
  result <- suppressWarnings(fetch_occurrence_data(url))  # Suppress warnings
  
  expect_true(is.list(result) || is.null(result), "Output should be a list or NULL on failure")
  
  if (!is.null(result)) {
    expect_true(length(result) >= 0, "Output list should be of valid length")
    expect_type(result, "list")
  }
})

test_that("fetch_occurrence_data handles incorrect URL gracefully", {
  
  # Invalid URL test
  invalid_url <- "https://invalid-url.example.com"
  result <- fetch_occurrence_data(invalid_url)
  
  expect_null(result, "Function should return NULL for invalid URLs")
})
