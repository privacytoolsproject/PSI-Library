library(PSIlence)
context("histogram")

# test accuracy and epsilon calculation for stability mechanism
test_that('histogram getAccuracy and getEpsilon return approximately correct values for stability mechanism', {
	val1 <- round(histogram.getAccuracy(mechanism = 'mechanismStability', epsilon=0.2, delta=10^-6, sensitivity=2))
	val2 <- round(histogram.getEpsilon(mechanism = 'mechanismStability', accuracy=2, delta=10^-6, sensitivity=2))
	expect_equal(val1, 176)
	expect_equal(val2, 35)
})

# test accuracy and epsilon calculation for laplace mechanism 
test_that('histogram getAccuracy and getEpsilon return approximately correct values for laplace mechanism', {
	val1 <- round(histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=0.2, sensitivity=2))
	val2 <- round(histogram.getEpsilon(mechanism = 'mechanismLaplace', accuracy=0.5, sensitivity=2))
	expect_equal(val1, 30)
	expect_equal(val2, 12)
})




# enter a nonsense variable type and expect error
test_that('expect stability mechanism for unknown variable type', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n.bins <- 16
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-9
    
    expect_error(dpHistogram$new(var.type='number data', variable="educ", n=my_n, epsilon=my_epsilon, 
                                 n.bins=my_n.bins, delta=my_delta, rng=c(0,16)), 
                 "Please enter a data type of 'numeric', 'integer', 'logical', or 'character'")
})





# test determineMechanism
# 1) if bines are entered, should be Laplace
test_that('histogram with bins entered', {
    library(datasets)
    data(esoph)
    
    my_n <- 88
    my_epsilon <- 1
    my_bins <- c("0-9g/day", "10-19", "20-29", "30+")
    
    catHistogram <- dpHistogram(var.type='character', variable='tobgp', n=my_n, epsilon=my_epsilon, bins=my_bins)
    catHistogram$release(esoph)
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(catHistogram$epsilon, my_epsilon)
    expect_equal(catHistogram$accuracy, askAccuracy)
})

# 2) if logical variable, should be Laplace
test_that('histogram release has expected dimensions and accuracy for logical variable with impute = false (laplace mechanism)', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    dp.histogram <- dpHistogram$new(var.type='logical', variable="sex", n=my_n, epsilon=my_epsilon)
    dp.histogram$release(PUMS5extract10000)
    expect_equal(length(dp.histogram$result$release), 3) # there should be 3 bins when impute = FALSE: 0,1,NA
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(dp.histogram$epsilon, my_epsilon)
    expect_equal(dim(dp.histogram$result$interval), c(3,2))
    expect_equal(dp.histogram$accuracy, askAccuracy)
})

# 3) if character variable and no bins entered, should be Stability
test_that('histogram on categorical data', {
    library(datasets)
    data(esoph)
    
    my_n <- 88
    my_epsilon <- 1
    my_delta <- 10^-4
    
    catHistogram <- dpHistogram(var.type='character', variable='tobgp', n=my_n, epsilon=my_epsilon, delta=my_delta)
    catHistogram$release(esoph)
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismStability', epsilon=my_epsilon, delta=my_delta, sensitivity=2)
    expect_equal(catHistogram$epsilon, my_epsilon)
    expect_equal(catHistogram$accuracy, askAccuracy)
})

# 4) if numeric and number of bins and range entered, should be Laplace
test_that('histogram releases have expected dimensions for Laplace mechanism', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n.bins <- 16
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    dp.histogram <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, rng=c(1,16))
    dp.histogram$release(PUMS5extract10000)
    expect_equal(length(dp.histogram$result$release), 16)
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(dp.histogram$epsilon, my_epsilon)
    expect_equal(dim(dp.histogram$result$interval), c(16,2))
    expect_equal(dp.histogram$accuracy, askAccuracy)
})

# 5) if numeric and number of bins entered without a range, should be Stability
test_that('histogram has expected accuracy for stability mechanism', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n.bins <- 16
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-9
    
    dp.histogram2 <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, delta=my_delta)
    dp.histogram2$release(PUMS5extract10000)
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismStability', epsilon=my_epsilon, delta=my_delta, sensitivity=2)
    expect_equal(dp.histogram2$epsilon, my_epsilon)
    expect_equal(dp.histogram2$accuracy, askAccuracy)
})

# 6) If numeric and number of bins not entered, expect error
test_that('histogram releases have expected dimensions for Laplace mechanism', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    expect_error(dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, rng=c(1,16), 'number of bins or granularity must be specified'))
})





# test on the stability mechanism: should return error if delta < 1/n^2
test_that('stability mechanism returns error if delta is >= 1/n', {
	data(PUMS5extract10000, package = "PSIlence")
	
	my_n.bins <- 16
	my_n <- 10000
	my_epsilon <- 0.1
	my_delta <- 0.1 # set delta to > 1/n^2
	
	dp.histogram2 <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, delta=my_delta)
	expect_error(dp.histogram2$release(PUMS5extract10000), "A delta value on the order of 1/n\\^2 is a privacy risk, as it allows for additional data to leak beyond the privacy parameter epsilon. Choose a smaller value for delta to maintain your privacy guarantee.")
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
    data(PUMS5extract10000, package = "PSIlence")
    
    my_bins <- c(0,20,30,40,50,60,70,80,90,100)
    
    expected_number_of_bins <- 9
    
    my_n <- 10000
    my_epsilon <- 0.1
    
    dp.histogram <- dpHistogram$new(var.type='numeric', variable="age", n=my_n, epsilon=my_epsilon, bins=my_bins, n.bins = 9)
    dp.histogram$release(PUMS5extract10000)
    expect_equal(length(dp.histogram$result$release), expected_number_of_bins)
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(dp.histogram$epsilon, my_epsilon)
    expect_equal(dim(dp.histogram$result$interval), c(expected_number_of_bins,2))
    expect_equal(dp.histogram$accuracy, askAccuracy)
})

# character
test_that('histogram on categorical data with bins entered', {
    library(datasets)
    data(esoph)
    
    my_n <- 88
    my_epsilon <- 1
    my_bins <- c("0-9g/day", "10-19", "20-29", "30+")
    
    catHistogram <- dpHistogram(var.type='character', variable='tobgp', n=my_n, epsilon=my_epsilon, bins=my_bins)
    catHistogram$release(esoph)
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(catHistogram$epsilon, my_epsilon)
    expect_equal(catHistogram$accuracy, askAccuracy)
})

# 2. enter bins that are not of the correct variable type
# expect error saying character bins cannot be entered for numeric variables
test_that('test on determineBins - get error when you enter character bins for numeric variable', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_bins <- c("should", "not", "be", "character", "bins")
    
    expected_number_of_bins <- 9
    
    my_n <- 10000
    my_epsilon <- 0.1
    
    expect_error(dpHistogram$new(var.type='numeric', variable="age", n=my_n, epsilon=my_epsilon, bins=my_bins), 
                 'Bins must be numeric for a numeric variable')
})

# expect error saying numeric bins cannot be entered for character variables
test_that('test on determineBins - get error when you enter numeric bins for character variable', {
    library(datasets)
    data(esoph)
    
    my_bins <- c(1,2,3,4,5)
    
    my_n <- 88
    my_epsilon <- 1
    
    expect_error(dpHistogram(var.type='character', variable='tobgp', n=my_n, epsilon=my_epsilon, bins=my_bins), 
                 'Bins must be of type `character` for a variable of type `character`')
})

# expect error saying logical bins must be 0,1,NA
test_that('test on determineBins - get error when you enter incorrect bins for logical variable', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_bins <- c("wrong", "bins")
    
    my_n <- 10000
    my_epsilon <- 0.1
    
    expect_error(dpHistogram$new(var.type='logical', variable="sex", n=my_n, epsilon=my_epsilon, bins=my_bins), 
                 'Histogram bins for a logical variable may only be 0, 1, or NA')
})

# 3. enter both bins and a range
# expect an error that says you entered both, and the code is defaulting to the bins entered
test_that('test on determineBins - get error when you enter both bins and a range', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_bins <- c(0,20,30,40,50,60,70,80,90,100)
    
    expected_number_of_bins <- 9
    
    my_n <- 10000
    my_epsilon <- 0.1
    
    expect_warning(dpHistogram$new(var.type='numeric', variable="age", n=my_n, epsilon=my_epsilon, bins=my_bins, rng=c(0.5,16.5)), "You have entered both bins and a data range, when you do not need both. Default is to use the bins that have been entered. If you would like to use the range, please enter the range and the desired number of bins and omit the bins.")
})

# 4. get correct bins for logical variable with impute = true or false
# no imputation
test_that('histogram release has expected dimensions and accuracy for logical variable with impute = false (laplace mechanism)', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    dp.histogram <- dpHistogram$new(var.type='logical', variable="sex", n=my_n, epsilon=my_epsilon)
    dp.histogram$release(PUMS5extract10000)
    expect_equal(length(dp.histogram$result$release), 3) # there should be 3 bins when impute = FALSE: 0,1,NA
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(dp.histogram$epsilon, my_epsilon)
    expect_equal(dim(dp.histogram$result$interval), c(3,2))
    expect_equal(dp.histogram$accuracy, askAccuracy)
})

# with imputation
test_that('histogram release has expected dimensions and accuracy for logical variable with impute = true (laplace mechanism)', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    dp.histogram <- dpHistogram$new(var.type='logical', variable="sex", n=my_n, epsilon=my_epsilon, impute = TRUE)
    dp.histogram$release(PUMS5extract10000)
    expect_equal(length(dp.histogram$result$release), 2) # there should be 2 bins when impute = TRUE: 0,1
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(dp.histogram$epsilon, my_epsilon)
    expect_equal(dim(dp.histogram$result$interval), c(2,2))
    expect_equal(dp.histogram$accuracy, askAccuracy)
})

# with imputation and the variable has NA values
test_that('histogram release has expected dimensions and accuracy for manually created logical variable with impute = true (laplace mechanism)', {
    logicalVar_withNA <- c(1,0,1,1,1,0,1,0,0,NA,1,0,NA,1,0,0,1,NA,NA,1,1,0,1,0,1,0)
    dataLog <- data.frame(logicalVar_withNA)
    
    my_n <- 26
    my_epsilon <- 1
    my_delta <- 10^-3
    
    dp.histogram <- dpHistogram$new(var.type='logical', variable="logicalVar_withNA", n=my_n, epsilon=my_epsilon, impute = TRUE)
    dp.histogram$release(dataLog)
    expect_equal(length(dp.histogram$result$release), 2) # there should be 2 bins when impute = TRUE: 0,1
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(dp.histogram$epsilon, my_epsilon)
    expect_equal(dim(dp.histogram$result$interval), c(2,2))
    expect_equal(dp.histogram$accuracy, askAccuracy)
})

# 5. get correct number of bins when numeric range and number of bins are entered, or granularity is entered
# number of bins entered
test_that('histogram releases have expected number of bins', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n.bins <- 16
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    dp.histogram <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, rng=c(0,16))
    dp.histogram$release(PUMS5extract10000)
    expect_equal(length(dp.histogram$result$release), 16)
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(dp.histogram$epsilon, my_epsilon)
    expect_equal(dim(dp.histogram$result$interval), c(16,2))
    expect_equal(dp.histogram$accuracy, askAccuracy)
})

# granularity entered
test_that('histogram releases have expected dimensions for Laplace mechanism', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_granularity <- 1000
    my_n <- 10000
    my_epsilon <- 0.1
    
    dp.histogram <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, granularity=my_granularity, rng=c(0,16))
    dp.histogram$release(PUMS5extract10000)
    expect_equal(length(dp.histogram$result$release), 10)
    
    askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', epsilon=my_epsilon, sensitivity=2)
    expect_equal(dp.histogram$epsilon, my_epsilon)
    expect_equal(dim(dp.histogram$result$interval), c(10,2))
    expect_equal(dp.histogram$accuracy, askAccuracy)
})

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {
    my_granularity <- 1000
    my_epsilon <- 0.1
    expect_error(dpHistogram$new(var.type='numeric', variable="educ", n=-1, epsilon=my_epsilon, granularity=my_granularity, rng=c(0,16)),
                 "n must be a positive whole number")
    expect_error(dpHistogram$new(var.type='numeric', variable="educ", n=0.5, epsilon=my_epsilon, granularity=my_granularity, rng=c(0,16)),
                 "n must be a positive whole number")
})

# make sure correct errors are thrown with incorrect values of nBins
test_that('errors thrown for incorrect values of nBins', {
    data(PUMS5extract10000, package = "PSIlence")
    
    my_n <- 10000
    my_epsilon <- 0.1
    my_delta <- 10^-6
    
    # expect warning and number of bins set to next-highest integer if user enters non-integer value
    my_n.bins <- 16.5
    expect_warning(dp.histogram <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, rng=c(0,16)), 'non-integer value for number of bins converted to next highest integer value')
    dp.histogram$release(PUMS5extract10000)
    expect_equal(length(dp.histogram$result$release), 17)
    
    # expect error if number of bins is less than 2
    my_n.bins <- 1
    expect_error(dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, rng=c(0,16)), 'number of bins must be at least 2')
})
