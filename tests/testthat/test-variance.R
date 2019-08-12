library(PSIlence)
context("variance")

data(PUMS5extract10000)

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {
    myEpsilon <- 0.1
    expect_error(dpVariance$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=-1, epsilon=myEpsilon, rng=c(18,93)),
                 "n must be a positive whole number")
})

# make sure you do not have to enter range for a logical variable
test_that('range checks throw correct warning', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    myEpsilon <- 0.1
    
    dpVar <- dpVariance$new(mechanism='mechanismLaplace', variable='sex', varType='logical', n=my_n, epsilon=myEpsilon)
    dpVar$release(PUMS5extract10000)
    
    expect_equal(length(dpVar$result$release), 1)
    expect_equal(dpVar$epsilon, myEpsilon)
})

# make sure error is thrown when dimension of range entered is incorrect
test_that('range checks throw correct warning', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    myEpsilon <- 0.1
    
    expect_error(dpVariance$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=my_n, epsilon=myEpsilon, rng=c(100)), 
                 "range argument in error: requires upper and lower values as vector of length 2.")
    
    dpVar <- dpVariance$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=my_n, epsilon=myEpsilon, rng=c(0,100))
    dpVar$release(PUMS5extract10000)
    expect_equal(length(dpVar$result$release), 1)
    expect_equal(dpVar$epsilon, myEpsilon)
})

# check for correct errors when imputation range is outside of entered range
test_that('error messages when imputation range is outside of data range', {
    my_n <- 10000
    myEpsilon <- 0.1
    my_rng <- c(18,93)
    
    expect_warning(dpVariance$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=my_n, epsilon=myEpsilon, rng=my_rng, imputeRng=c(0,93)),
                   'Lower bound of imputation range is outside of the data range. Setting lower bound of the imputation range to the lower bound of the data range.')
    
    expect_warning(dpVariance$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=my_n, epsilon=myEpsilon, rng=my_rng, imputeRng=c(18,200)),
                   'Upper bound of imputation range is outside of the data range. Setting upper bound of the imputation range to the upper bound of the data range.')
    
    expect_warning(dpVariance$new(mechanism='mechanismLaplace', variable='sex', varType='logical', n=my_n, epsilon=myEpsilon, imputeRng=c(2,3)),
                   'Imputation range entered for variable that is not of numeric or integer type. Setting imputation range to data range.')
    
    expect_warning(dpVariance$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=my_n, epsilon=myEpsilon, rng=my_rng, imputeRng=c('wrong','type')),
                   'Imputation range for a numeric variable must be numeric. Setting imputation range to data range.')
})
