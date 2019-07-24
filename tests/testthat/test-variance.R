library(PSIlence)
context("covariance")

data(PUMS5extract10000)

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {
    my_epsilon <- 0.1
    my_delta <- 10^-6
    expect_error(dpVariance$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=-1, epsilon=my_epsilon, rng=c(100)),
                 "n must be a positive whole number")
    expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', var.type='numeric', n=0.5, epsilon=my_epsilon, rng=c(-10,0,100)),
                 "n must be a positive whole number")
})
