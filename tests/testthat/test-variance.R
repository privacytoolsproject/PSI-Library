library(PSIlence)
context("covariance")

data(PUMS5extract10000)

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {
    my_epsilon <- 0.1
    my_delta <- 10^-6
    expect_error(dpVariance$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=-1, epsilon=my_epsilon, rng=c(18,93)),
                 "n must be a positive whole number")
    expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=0.5, epsilon=my_epsilon, rng=c(18,93)),
                 "n must be a positive whole number")
})

# make sure you do not have to enter range for a logical variable
test_that('range checks throw correct warning', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    
    dp.variance <- dpVariance$new(mechanism='mechanismLaplace', variable='sex', var.type='logical', n=my_n, epsilon=my_epsilon)
    dp.variance$release(PUMS5extract10000)
    
    expect_equal(length(dp.variance$result$release), 1)
    expect_equal(dp.variance$epsilon, my_epsilon)
})

# make sure error is thrown when dimension of range entered is incorrect
test_that('range checks throw correct warning', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    expect_error(dpVariance$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=c(100)), 
                 "range argument in error: requires upper and lower values as vector of length 2.")
    expect_warning(dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=c(-10,0,100)), 
                   "range argument supplied has more than two values.  Will proceed using min and max values as range.")
    
    dp.variance <- dpVariance$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=c(0,100))
    dp.variance$release(PUMS5extract10000)
    expect_equal(length(dp.variance$result$release), 1)
    expect_equal(dp.variance$epsilon, my_epsilon)
})

# test accuracy and epsilon calculations 
test_that('getAccuracy and getEpsilon return approximately correct values for laplace mechanism', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_rng <- c(18,93)
    my_sens <- (my_n - 1) / my_n^2 * diff(my_rng)^2
    
    val1 <- round(histogram.getAccuracy(mechanism = 'mechanismLaplace', n=my_n, epsilon=0.1, sensitivity=my_sens))
    val2 <- round(histogram.getEpsilon(mechanism = 'mechanismLaplace', n=my_n, accuracy=20, sensitivity=my_sens), digits = 1)
    expect_equal(val1, 17)
    expect_equal(val2, 0.1)
})

# check for correct errors when imputation range is outside of entered range
test_that('error messages when imputation range is outside of data range', {
    my_n <- 10000
    my_epsilon <- 0.1
    my_rng <- c(18,93)
    
    expect_warning(dpVariance$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=my_rng, impute.rng=c(0,93)),
                   'Lower bound of imputation range is outside of the data range. Setting lower bound of the imputation range to the lower bound of the data range.')
    
    expect_warning(dpVariance$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=my_rng, impute.rng=c(18,200)),
                   'Upper bound of imputation range is outside of the data range. Setting upper bound of the imputation range to the upper bound of the data range.')
    
    expect_warning(dpVariance$new(mechanism='mechanismLaplace', variable='sex', var.type='logical', n=my_n, epsilon=my_epsilon, impute.rng=c(2,3)),
                   'Imputation range entered for variable that is not of numeric or integer type. Setting imputation range to data range.')
    
    expect_warning(dpVariance$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=my_n, epsilon=my_epsilon, rng=my_rng, impute.rng=c('wrong','type')),
                   'Imputation range for a numeric variable must be numeric. Setting imputation range to data range.')
})
