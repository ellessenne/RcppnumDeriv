testthat::test_that("genD returns same results of numDeriv", {
  func <- function(x) c(x[1], x[1], x[2]^2)
  xxx <- c(2, 2, 5)
  original <- numDeriv::genD(func, xxx)
  new <- RcppnumDeriv::genD(func, xxx)

  testthat::expect_equal(object = new$D, expected = original$D)
})

# library(benchmark)
# microbenchmark::microbenchmark(
#   "original" = numDeriv::genD(func, xxx),
#   "new" = RcppnumDeriv::genD(func, xxx),
#   times = 1000
# )
