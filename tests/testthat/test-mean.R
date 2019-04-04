library(PSIlence)
context("mean")

test_that('range checks throw correct warning', {
  data(PUMS5extract10000)

  my_n <- 10000
  my_epsilon <- 0.1
  my_delta <- 10^-6

  expect_error(dpMean$new(mechanism='mechanismLaplace', var.type='numeric', variable='age', n=my_n, epsilon=my_epsilon, rng=c(100)), 
               "range argument in error: requires upper and lower values as vector of length 2.")

  expect_warning(dpMean$new(mechanism='mechanismLaplace', var.type='numeric', variable='age', n=my_n, epsilon=my_epsilon, rng=c(-10,0,100)),
                 "range argument supplied has more than two values.  Will proceed using min and max values as range.")

  dp.mean3 <- dpMean$new(mechanism='mechanismLaplace', var.type='numeric', variable='age', n=my_n, epsilon=my_epsilon, delta=my_delta, rng=c(0,100))
  dp.mean3$release(PUMS5extract10000)
  expect_equal(length(dp.mean3$result$release), 1)
  expect_equal(dp.mean3$epsilon, my_epsilon)
  expect_equal(length(dp.mean3$result$interval), 2)
})