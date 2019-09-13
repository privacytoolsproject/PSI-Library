library(PSIlence)
context("variance")

data(PUMS5extract10000)

# test sensitivity function
test_that('sensitivity function is consistent with intended implementation', {
  expect_equal(varianceSensitivity(2, c(0,10)),50)
  expect_equal(varianceSensitivity(5, c(5,10)),5)
})

# test accuracy, epsilon, and sensitivity calculations 
test_that('variancee getAccuracy and getEpsilon return approximately correct values for laplace mechanism', {
    # test sensitivity and accuracy
    nTest <- 10000
    epsilonTest <- 0.1
    
    dpVar <- dpVariance$new(variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=c(0,100))
    dpVar$release(PUMS5extract10000)
    
    sens <- round((dpVar$result$epsilon * dpVar$result$accuracy) / log(1/0.05))
    acc <- round(dpVar$result$accuracy)

    expect_equal(sens, 1)
    expect_equal(acc, 30)
    
    # test accuracy
    accuracyTest <- 30
    
    dpVar2 <- dpVariance$new(variable='age', varType='numeric', n=nTest, accuracy=accuracyTest, rng=c(0,100))
    dpVar2$release(PUMS5extract10000)
    
    epsilon <- round(dpVar2$result$epsilon, digits = 1)
    expect_equal(epsilon, 0.1)
})

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    expect_error(dpVariance$new(variable='age', varType='numeric', n=-1, epsilon=epsilonTest, rng=c(18,93)),
                 "n must be a positive whole number")
})

# make sure you do not have to enter range for a logical variable
test_that('range checks throw correct warning', {
    data(PUMS5extract10000, package = "PSIlence")
    
    nTest <- 10000
    epsilonTest <- 0.1
    
    dpVar <- dpVariance$new(variable='sex', varType='logical', n=nTest, epsilon=epsilonTest)
    dpVar$release(PUMS5extract10000)
    
    expect_equal(length(dpVar$result$release), 1)
    expect_equal(dpVar$epsilon, epsilonTest)
})

# make sure error is thrown when dimension of range entered is incorrect
test_that('range checks throw correct warning', {
    data(PUMS5extract10000, package = "PSIlence")
    
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    expect_error(dpVariance$new(variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=c(100)), 
                 "range argument in error: requires upper and lower values as vector of length 2.")

    expect_warning(dpVariance$new(variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=c(-10,0,100)), 
                   "range argument supplied has more than two values.  Will proceed using min and max values as range.")
    
    dpVar <- dpVariance$new(variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=c(0,100))
    dpVar$release(PUMS5extract10000)
    expect_equal(length(dpVar$result$release), 1)
    expect_equal(dpVar$epsilon, epsilonTest)
})

# check for correct errors when imputation range is outside of entered range
test_that('error messages when imputation range is outside of data range', {
    nTest <- 10000
    epsilonTest <- 0.1
    rngTest <- c(18,93)
    
    expect_warning(dpVariance$new(variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=rngTest, imputeRng=c(0,93)),
                   'Lower bound of imputation range is outside of the data range. Setting lower bound of the imputation range to the lower bound of the data range.')
    
    expect_warning(dpVariance$new(variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=rngTest, imputeRng=c(18,200)),
                   'Upper bound of imputation range is outside of the data range. Setting upper bound of the imputation range to the upper bound of the data range.')
    
    expect_warning(dpVariance$new(variable='sex', varType='logical', n=nTest, epsilon=epsilonTest, imputeRng=c(2,3)),
                   'Imputation range entered for variable that is not of numeric or integer type. Setting imputation range to data range.')
    
    expect_warning(dpVariance$new(variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=rngTest, imputeRng=c('wrong','type')),
                   'Imputation range for a numeric variable must be numeric. Setting imputation range to data range.')
})

test_that('output has correct dimensions', {
  x <- c(0,1,5,9,3,10)
  n <- length(x)
  data <- data.frame(x)
  
  dpVar <- dpVariance$new(variable='x', varType='numeric', n=length(x), epsilon=1, rng=c(0,10))
  out <- dpVar$release(data)

  expect_true(is.numeric(out$release))
  expect_true(!is.null(out$accuracy))
  expect_equal(out$variable, 'x')
})
