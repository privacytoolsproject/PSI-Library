library(PSIlence)
context("mean")

test_that('range checks throw correct warning', {
  data(PUMS5extract10000, package = "PSIlence")

  my_n <- 10000
  my_epsilon <- 0.1
  my_delta <- 10^-6

  expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=c(100)), 
               "range argument in error: requires upper and lower values as vector of length 2.")
  expect_warning(dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=c(-10,0,100)), 
                 "range argument supplied has more than two values.  Will proceed using min and max values as range.")

  dp.mean <- dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, delta=my_delta, rng=c(0,100))
  dp.mean$release(PUMS5extract10000)
  expect_equal(length(dp.mean$result$release), 1)
  expect_equal(dp.mean$epsilon, my_epsilon)
  expect_equal(length(dp.mean$result$interval), 2)
})

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {
    my_epsilon <- 0.1
    my_delta <- 10^-6
    expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=-1, epsilon=my_epsilon, rng=c(100)),
                 "n must be a positive whole number")
    expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=0.5, epsilon=my_epsilon, rng=c(-10,0,100)),
                 "n must be a positive whole number")
})

# make sure you do not have to enter range for a logical variable
test_that('range checks throw correct warning', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    dp.mean <- dpMean$new(mechanism='mechanismLaplace', variable='sex', var.type='logical', n=my_n, epsilon=my_epsilon, delta=my_delta)
    dp.mean$release(PUMS5extract10000)
    
    expect_equal(length(dp.mean$result$release), 1)
    expect_equal(dp.mean$epsilon, my_epsilon)
    expect_equal(length(dp.mean$result$interval), 2)
})

# test accuracy and epsilon calculations 
test_that('histogram getAccuracy and getEpsilon return approximately correct values for laplace mechanism', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    my_rng <- c(18,93)
    my_acc <- 0.25
    
    dp.mean <- dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=my_rng)
    dp.mean$release(PUMS5extract10000)
    
    acc <- round(dp.mean$result$accuracy, digits = 5)
    
    expect_equal(acc, 0.22468)
    
    dp.mean2 <- dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, accuracy=my_acc, rng=my_rng)
    dp.mean2$release(PUMS5extract10000)
    
    ep <- round(dp.mean2$result$epsilon, digits = 2)
    
    expect_equal(ep, 0.09)
})
