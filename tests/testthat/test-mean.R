test_that('range checks throw correct warning', {
  library(PSIlence)
  data(PUMS5extract10000, package = "PSIlence")

  my_n <- 10000
  my_epsilon <- 0.1
  my_delta <- 10^-6

  dp.mean1 <- dpMean$new(mechanism='mechanismLaplace', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=c(100))
  expect_error(dp.mean1$release(PUMS5extract10000$age), "range argument in error: requires upper and lower values as vector of length 2.")


  dp.mean2 <- dpMean$new(mechanism='mechanismLaplace', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=c(-10,0,100))
  expect_warning(dp.mean2$release(PUMS5extract10000$age), "range argument supplied has more than two values.  Will proceed using min and max values as range.")

  dp.mean3 <- dpMean$new(mechanism='mechanismLaplace', var.type='numeric', n=my_n, epsilon=my_epsilon, delta=my_delta, rng=c(0,100))
  dp.mean3$release(PUMS5extract10000$age)
  expect_equal(length(dp.mean3$result$release), 1)
  expect_equal(dp.mean3$epsilon, my_epsilon)
  expect_equal(length(dp.mean3$result$interval), 2)
})