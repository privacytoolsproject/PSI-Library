context("histogram")

library(datasets)
data(esoph)
data(PUMS5extract10000, package = "PSIlence")

# test accuracy and epsilon calculation for stability mechanism
test_that('histogram getAccuracy and getEpsilon return approximately correct values for stability mechanism', {
	val1 <- round(histogramGetAccuracy(mechanism = 'mechanismStability', epsilon=0.2, delta=10^-6, sensitivity=2))
	val2 <- round(histogramGetEpsilon(mechanism = 'mechanismStability', accuracy=2, delta=10^-6, sensitivity=2))
	expect_equal(val1, 176)
	expect_equal(val2, 35)
})

# test accuracy and epsilon calculation for laplace mechanism 
test_that('histogram getAccuracy and getEpsilon return approximately correct values for laplace mechanism', {
	val1 <- round(histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=0.2, sensitivity=2))
	val2 <- round(histogramGetEpsilon(mechanism = 'mechanismLaplace', accuracy=0.5, sensitivity=2))
	expect_equal(val1, 30)
	expect_equal(val2, 12)
})

# enter a nonsense variable type and expect error
test_that('expect stability mechanism for unknown variable type', {
    data(PUMS5extract10000, package = "PSIlence")
    
    nBinsTest <- 16
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-9
    
    expect_error(dpHistogram$new(varType='number data', variable="educ", n=nTest, epsilon=epsilonTest, 
                                 nBins=nBinsTest, delta=deltaTest, rng=c(0,16)), 
                 "Variable type number data should be one of numeric, integer, logical, character", fixed=TRUE)
})

# test determineMechanism
# 1) if bins are entered, should be Laplace
test_that('histogram with bins entered', {

    nTest <- 88
    epsilonTest <- 1
    binsTest <- c("0-9g/day", "10-19", "20-29", "30+")
    
    catHistogram <- dpHistogram(varType='character', variable='tobgp', n=nTest, epsilon=epsilonTest, bins=binsTest)
    catHistogram$release(esoph)
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(catHistogram$epsilon, epsilonTest)
    expect_equal(catHistogram$accuracy, askAccuracy)
})

# 2) if logical variable, should be Laplace
test_that('histogram release has expected dimensions and accuracy for logical variable with impute = false (laplace mechanism)', {
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    dpHist <- dpHistogram$new(varType='logical', variable="sex", n=nTest, epsilon=epsilonTest)
    dpHist$release(PUMS5extract10000)
    expect_equal(length(dpHist$result$release), 3) # there should be 3 bins when impute = FALSE: 0,1,NA
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(dpHist$epsilon, epsilonTest)
    expect_equal(dim(dpHist$result$interval), c(3,2))
    expect_equal(dpHist$accuracy, askAccuracy)
})

# 3) if character variable and no bins entered, should be Stability
test_that('histogram on categorical data', {
    
    nTest <- 88
    epsilonTest <- 1
    deltaTest <- 10^-4
    
    catHistogram <- dpHistogram(varType='character', variable='tobgp', n=nTest, epsilon=epsilonTest, delta=deltaTest)
    catHistogram$release(esoph)
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismStability', epsilon=epsilonTest, delta=deltaTest, sensitivity=2)
    expect_equal(catHistogram$epsilon, epsilonTest)
    expect_equal(catHistogram$accuracy, askAccuracy)
})

# 4) if numeric and number of bins and range entered, should be Laplace
test_that('histogram releases have expected dimensions for Laplace mechanism', {
    
    nBinsTest <- 16
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    dpHist <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, nBins=nBinsTest, rng=c(1,16))
    dpHist$release(PUMS5extract10000)
    expect_equal(length(dpHist$result$release), 16)
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(dpHist$epsilon, epsilonTest)
    expect_equal(dim(dpHist$result$interval), c(16,2))
    expect_equal(dpHist$accuracy, askAccuracy)
})

# 5) if numeric and number of bins entered without a range, should be Stability
test_that('histogram has expected accuracy for stability mechanism', {

    nBinsTest <- 16
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-9
    
    dpHist2 <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, nBins=nBinsTest, delta=deltaTest)
    dpHist2$release(PUMS5extract10000)
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismStability', epsilon=epsilonTest, delta=deltaTest, sensitivity=2)
    expect_equal(dpHist2$epsilon, epsilonTest)
    expect_equal(dpHist2$accuracy, askAccuracy)
})

# 6) If numeric and number of bins not entered, expect error
test_that('histogram releases have expected dimensions for Laplace mechanism', {
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    expect_error(dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, rng=c(1,16), 'number of bins or granularity must be specified'))
})

# test on the stability mechanism: should return error if delta < 1/n^2
test_that('stability mechanism returns error if delta is >= 1/n', {
	nBinsTest <- 16
	nTest <- 10000
	epsilonTest <- 0.1
	deltaTest <- 0.1 # set delta to > 1/n^2
	
	dpHist2 <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, nBins=nBinsTest, delta=deltaTest)
	expect_error(dpHist2$release(PUMS5extract10000), "A delta value on the order of 1/n\\^2 is a privacy risk, as it allows for additional data to leak beyond the privacy parameter epsilon. Choose a smaller value for delta to maintain your privacy guarantee.")
})

# tests for determineBins

# enter bins, check errors
# 1. enter bins that are correct
# 2. enter bins that are not of the correct variable type
# 3. enter both numeric bins and a bin range
# 4. get correct bins for logical variable with impute = true or false
# 5. get correct number of bins when numeric range and number of bins are entered, or granularity is entered

# 1. enter bins that are correct
# expect the code to run and check the output
# numeric
test_that('test on determineBins - ensure correct number of bins when bins are entered correctly', {
    binsTest <- c(0,20,30,40,50,60,70,80,90,100)
    
    expectedNumberOfBins <- 9
    
    nTest <- 10000
    epsilonTest <- 0.1
    
    dpHist <- dpHistogram$new(varType='numeric', variable="age", n=nTest, epsilon=epsilonTest, bins=binsTest, nBins = 9)
    dpHist$release(PUMS5extract10000)
    expect_equal(length(dpHist$result$release), expectedNumberOfBins)
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(dpHist$epsilon, epsilonTest)
    expect_equal(dim(dpHist$result$interval), c(expectedNumberOfBins,2))
    expect_equal(dpHist$accuracy, askAccuracy)
})

# character
test_that('histogram on categorical data with bins entered', {
    
    nTest <- 88
    epsilonTest <- 1
    binsTest <- c("0-9g/day", "10-19", "20-29", "30+")
    
    catHistogram <- dpHistogram(varType='character', variable='tobgp', n=nTest, epsilon=epsilonTest, bins=binsTest)
    catHistogram$release(esoph)
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(catHistogram$epsilon, epsilonTest)
    expect_equal(catHistogram$accuracy, askAccuracy)
})

# 2. enter bins that are not of the correct variable type
# expect error saying character bins cannot be entered for numeric variables
test_that('test on determineBins - get error when you enter character bins for numeric variable', {
    
    binsTest <- c("should", "not", "be", "character", "bins")
    
    expectedNumberOfBins <- 9
    
    nTest <- 10000
    epsilonTest <- 0.1
    
    expect_error(dpHistogram$new(varType='numeric', variable="age", n=nTest, epsilon=epsilonTest, bins=binsTest), 
                 'Bins must be numeric for a numeric variable')
})

# expect error saying numeric bins cannot be entered for character variables
test_that('test on determineBins - get error when you enter numeric bins for character variable', {
    
    binsTest <- c(1,2,3,4,5)
    
    nTest <- 88
    epsilonTest <- 1
    
    expect_error(dpHistogram(varType='character', variable='tobgp', n=nTest, epsilon=epsilonTest, bins=binsTest), 
                 'Bins must be of type `character` for a variable of type `character`')
})

# expect error saying logical bins must be 0,1,NA
test_that('test on determineBins - get error when you enter incorrect bins for logical variable', {

    binsTest <- c("wrong", "bins")
    
    nTest <- 10000
    epsilonTest <- 0.1
    
    expect_error(dpHistogram$new(varType='logical', variable="sex", n=nTest, epsilon=epsilonTest, bins=binsTest), 
                 'Histogram bins for a logical variable may only be 0, 1, or NA')
})

# 3. enter both bins and a range
# expect an error that says you entered both, and the code is defaulting to the bins entered
test_that('test on determineBins - get error when you enter both bins and a range', {
    binsTest <- c(0,20,30,40,50,60,70,80,90,100)
    
    expectedNumberOfBins <- 9
    
    nTest <- 10000
    epsilonTest <- 0.1
    
    expect_warning(dpHistogram$new(varType='numeric', variable="age", n=nTest, epsilon=epsilonTest, bins=binsTest, rng=c(0.5,16.5)), "You have entered both bins and a data range, when you do not need both. Default is to use the bins that have been entered. If you would like to use the range, please enter the range and the desired number of bins and omit the bins.")
})

# 4. get correct bins for logical variable with impute = true or false
# no imputation
test_that('histogram release has expected dimensions and accuracy for logical variable with impute = false (laplace mechanism)', {

    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    dpHist <- dpHistogram$new(varType='logical', variable="sex", n=nTest, epsilon=epsilonTest)
    dpHist$release(PUMS5extract10000)
    expect_equal(length(dpHist$result$release), 3) # there should be 3 bins when impute = FALSE: 0,1,NA
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(dpHist$epsilon, epsilonTest)
    expect_equal(dim(dpHist$result$interval), c(3,2))
    expect_equal(dpHist$accuracy, askAccuracy)
})

# with imputation
test_that('histogram release has expected dimensions and accuracy for logical variable with impute = true (laplace mechanism)', {

    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    dpHist <- dpHistogram$new(varType='logical', variable="sex", n=nTest, epsilon=epsilonTest, impute = TRUE)
    dpHist$release(PUMS5extract10000)
    expect_equal(length(dpHist$result$release), 2) # there should be 2 bins when impute = TRUE: 0,1
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(dpHist$epsilon, epsilonTest)
    expect_equal(dim(dpHist$result$interval), c(2,2))
    expect_equal(dpHist$accuracy, askAccuracy)
})

# with imputation and the variable has NA values
test_that('histogram release has expected dimensions and accuracy for manually created logical variable with impute = true (laplace mechanism)', {

    logicalVar_withNA <- c(1,0,1,1,1,0,1,0,0,NA,1,0,NA,1,0,0,1,NA,NA,1,1,0,1,0,1,0)
    dataLog <- data.frame(logicalVar_withNA)
    
    nTest <- 26
    epsilonTest <- 1
    deltaTest <- 10^-3
    
    dpHist <- dpHistogram$new(varType='logical', variable="logicalVar_withNA", n=nTest, epsilon=epsilonTest, impute = TRUE)
    dpHist$release(dataLog)
    expect_equal(length(dpHist$result$release), 2) # there should be 2 bins when impute = TRUE: 0,1
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(dpHist$epsilon, epsilonTest)
    expect_equal(dim(dpHist$result$interval), c(2,2))
    expect_equal(dpHist$accuracy, askAccuracy)
})

# 5. get correct number of bins when numeric range and number of bins are entered, or granularity is entered
# number of bins entered
test_that('histogram releases have expected number of bins', {

    nBinsTest <- 16
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    dpHist <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, nBins=nBinsTest, rng=c(0,16))
    dpHist$release(PUMS5extract10000)
    expect_equal(length(dpHist$result$release), 16)
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(dpHist$epsilon, epsilonTest)
    expect_equal(dim(dpHist$result$interval), c(16,2))
    expect_equal(dpHist$accuracy, askAccuracy)
})

# granularity entered
test_that('histogram releases have expected dimensions for Laplace mechanism', {

    granularityTest <- 1000
    nTest <- 10000
    epsilonTest <- 0.1
    
    dpHist <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, granularity=granularityTest, rng=c(0,16))
    dpHist$release(PUMS5extract10000)
    expect_equal(length(dpHist$result$release), 10)
    
    askAccuracy <- histogramGetAccuracy(mechanism = 'mechanismLaplace', epsilon=epsilonTest, sensitivity=2)
    expect_equal(dpHist$epsilon, epsilonTest)
    expect_equal(dim(dpHist$result$interval), c(10,2))
    expect_equal(dpHist$accuracy, askAccuracy)
})

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {

    granularityTest <- 1000
    epsilonTest <- 0.1
    expect_error(dpHistogram$new(varType='numeric', variable="educ", n=-1, epsilon=epsilonTest, granularity=granularityTest, rng=c(0,16)),
                 "n must be a positive whole number")
    expect_error(dpHistogram$new(varType='numeric', variable="educ", n=0.5, epsilon=epsilonTest, granularity=granularityTest, rng=c(0,16)),
                 "n must be a positive whole number")
})

# make sure correct errors are thrown with incorrect values of nBins
test_that('errors thrown for incorrect values of nBins', {
  
    nTest <- 10000
    epsilonTest <- 0.1
    deltaTest <- 10^-6
    
    # expect warning and number of bins set to next-highest integer if user enters non-integer value
    nBinsTest <- 16.5
    expect_warning(dpHist <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, nBins=nBinsTest, rng=c(0,16)), 'non-integer value for number of bins converted to next highest integer value')
    dpHist$release(PUMS5extract10000)
    expect_equal(length(dpHist$result$release), 17)
    
    # expect error if number of bins is less than 2
    nBinsTest <- 1
    expect_error(dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, nBins=nBinsTest, rng=c(0,16)), 'number of bins must be at least 2')
})

# make sure delta value is correct, or correct error is thrown
test_that('check delta', {
    data(PUMS5extract10000, package = "PSIlence")
    
    granularityTest <- 1000
    nTest <- 10000
    epsilonTest <- 0.1
    
    # check that delta is set to 0 for histogram that uses laplace mechanism
    dpHist <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, granularity=granularityTest, rng=c(0,16))
    dpHist$release(PUMS5extract10000)
    expect_equal(dpHist$result$delta, 0)
    
    # check that warning is thrown when user entere delta value for histogram that uses laplace mechanism
    expect_warning(dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, granularity=granularityTest, rng=c(0,16), delta=10^-5), 'A delta parameter has been entered, but a mechanism that uses a delta value is not being used. Setting delta to 0.')
    
    # check that default value of delta set when stability mechanism used and delta value not entered
    nBinsTest <- 16
    
    dpHist2 <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, nBins=nBinsTest)
    dpHist2$release(PUMS5extract10000)
    expect_equal(dpHist2$result$delta, 2^-30)
    
    # check that the entered delta value is set as delta when stbaility mechanism used and delta value is entered
    dpHist3 <- dpHistogram$new(varType='numeric', variable="educ", n=nTest, epsilon=epsilonTest, nBins=nBinsTest, delta=10^-10)
    dpHist3$release(PUMS5extract10000)
    expect_equal(dpHist3$result$delta, 10^-10)
})
