context("mean")

data <- data(PUMS5extract10000, package = "PSIlence")

test_that('getAccuracy and getEpsilon return approximately correct values for laplace mechanism', {
    # test sensitivity and accuracy
    nTest <- 10000
    epsilonTest <- 0.1
    
    dpMeanTest <- dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=c(0,100))
    dpMeanTest$release(PUMS5extract10000)
    
    sens <- 100 / nTest
    acc <- round(dpMeanTest$result$accuracy, digits = 1)

    expect_equal(sens, 0.01)
    expect_equal(acc, 0.3)
    
    # test accuracy
    accuracyTest <- 0.3

    dpMeanTest2 <- dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, accuracy=accuracyTest, rng=c(0,100))
    dpMeanTest2$release(PUMS5extract10000)

    epsilon <- round(dpMeanTest2$result$epsilon, digits = 1)
    expect_equal(epsilon, 0.1)
})

test_that('range checks throw correct warning', {

  nTest <- 10000
  epsilonTest <- 0.1
  deltaTest <- 10^-6

  expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=c(100)))
  expect_warning(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=c(-10,0,100)))
  expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=NULL))
  
  dpMeanTest <- dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, delta=deltaTest, rng=c(0,100))
  dpMeanTest$release(PUMS5extract10000)
  expect_equal(length(dpMeanTest$result$release), 1)
  expect_equal(dpMeanTest$epsilon, epsilonTest)
  expect_equal(length(dpMeanTest$result$interval), 2)
})

test_that('error thrown when n not positive or whole number', {
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=-1, epsilon=epsilonTest, rng=c(0,100)),
                 "n must be a positive whole number")
    expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=0.5, epsilon=epsilonTest, rng=c(0,100)),
                 "n must be a positive whole number")
})

test_that('you do not have to enter range for a logical variable', {
    
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    dpMeanTest <- dpMean$new(mechanism='mechanismLaplace', variable='sex', varType='logical', n=nTest, epsilon=epsilonTest, delta=deltaTest)
    dpMeanTest$release(PUMS5extract10000)
    
    expect_equal(length(dpMeanTest$result$release), 1)
    expect_equal(dpMeanTest$epsilon, epsilonTest)
    expect_equal(length(dpMeanTest$result$interval), 2)
    expect_equal(length(dpMeanTest$result$histogram), 2)
})

test_that('error messages occur when imputation range is outside of data range', {
    nTest <- 10000
    epsilonTest <- 0.1
    rngTest <- c(18,93)
    
    expect_warning(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=rngTest, imputeRng=c(0,93)),
                   'Lower bound of imputation range is outside of the data range. Setting lower bound of the imputation range to the lower bound of the data range.')
    
    expect_warning(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=rngTest, imputeRng=c(18,200)),
                   'Upper bound of imputation range is outside of the data range. Setting upper bound of the imputation range to the upper bound of the data range.')
    
    expect_warning(dpMean$new(mechanism='mechanismLaplace', variable='sex', varType='logical', n=nTest, epsilon=epsilonTest, imputeRng=c(2,3)),
                   'Imputation range entered for variable that is not of numeric or integer type. Setting imputation range to data range.')
    
    expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest, rng=rngTest, imputeRng=c('wrong','type')))
    
    # make sure error thrown when range is null for variable that requires range
    expect_error(dpMean$new(mechanism='mechanismLaplace', variable='age', varType='numeric', n=nTest, epsilon=epsilonTest),
                 'Input range may not be empty.')
})
