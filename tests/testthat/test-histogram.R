library(PSIlence)
context("histogram")

test_that('histogram getAccuracy and getParameters return approximately correct values for stability mechanism', {
	val1 <- round(histogram.getAccuracy(mechanism = 'mechanismStability', n.bins=10, n=2000, epsilon=0.2, delta=10^-6))
	val2 <- round(histogram.getParameters(mechanism = 'mechanismStability', n.bins=20, n=2000, accuracy=0.5, delta=10^-6))
	expect_equal(val1, 175)
	expect_equal(val2, 70)
})

test_that('histogram getAccuracy and getParameters return approximately correct values for laplace mechanism', {
	val1 <- round(histogram.getAccuracy(mechanism = 'mechanismLaplace', n.bins=10, n=2000, epsilon=0.2, delta=10^-6))
	val2 <- round(histogram.getParameters(mechanism = 'mechanismLaplace', n.bins=20, n=2000, accuracy=0.5, delta=10^-6))
	expect_equal(val1, 30)
	expect_equal(val2, 12)
})

test_that('histogram releases have expected dimensions for lapalce mechanism', {
	data(PUMS5extract10000, package = "PSIlence")
	
	my_n.bins <- 16
	my_n <- 10000
	my_epsilon <- 0.1
	my_delta <- 10^-6
	
	dp.histogram <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, delta=my_delta, rng=c(0.5,16.5))
	dp.histogram$release(PUMS5extract10000)
	expect_equal(length(dp.histogram$result$release), 16)
	
	askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', n.bins=my_n.bins, n=my_n, epsilon=my_epsilon, delta=my_delta)
	expect_equal(dp.histogram$epsilon, my_epsilon)
	expect_equal(dim(dp.histogram$result$interval), c(16,2))
	expect_equal(dp.histogram$accuracy, askAccuracy)
})

test_that('histogram has expected accuracy for stability mechanism', {
	data(PUMS5extract10000, package = "PSIlence")
	
	my_n.bins <- 16
	my_n <- 10000
	my_epsilon <- 0.1
	my_delta <- 10^-6
	
	dp.histogram2 <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, delta=my_delta)
	dp.histogram2$release(PUMS5extract10000)
	
	askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismStability', n.bins=my_n.bins, n=my_n, epsilon=my_epsilon, delta=my_delta)
	expect_equal(dp.histogram2$epsilon, my_epsilon)
	expect_equal(dp.histogram2$accuracy, askAccuracy)
})

test_that('histogram release has expected dimensions and accuracy for logical variable with impute = false (laplace mechanism)', {
	data(PUMS5extract10000, package = "PSIlence")
	
	my_n <- 10000
	my_epsilon <- 0.1
	my_delta <- 10^-6
	
	dp.histogram <- dpHistogram$new(var.type='logical', variable="sex", n=my_n, epsilon=my_epsilon, delta=my_delta)
	dp.histogram$release(PUMS5extract10000)
	expect_equal(length(dp.histogram$result$release), 3) # there should be 3 bins when impute = FALSE: 0,1,NA
	
	askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', n.bins=my_n.bins, n=my_n, epsilon=my_epsilon, delta=my_delta)
	expect_equal(dp.histogram$epsilon, my_epsilon)
	expect_equal(dim(dp.histogram$result$interval), c(3,2))
	expect_equal(dp.histogram$accuracy, askAccuracy)
})

test_that('histogram release has expected dimensions and accuracy for logical variable with impute = true (laplace mechanism)', {
	data(PUMS5extract10000, package = "PSIlence")
	
	my_n <- 10000
	my_epsilon <- 0.1
	my_delta <- 10^-6
	
	dp.histogram <- dpHistogram$new(var.type='logical', variable="sex", n=my_n, epsilon=my_epsilon, delta=my_delta, impute = TRUE)
	dp.histogram$release(PUMS5extract10000)
	expect_equal(length(dp.histogram$result$release), 2) # there should be 2 bins when impute = TRUE: 0,1
	
	askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', n.bins=my_n.bins, n=my_n, epsilon=my_epsilon, delta=my_delta)
	expect_equal(dp.histogram$epsilon, my_epsilon)
	expect_equal(dim(dp.histogram$result$interval), c(2,2))
	expect_equal(dp.histogram$accuracy, askAccuracy)
})


test_that('stability mechanism returns error if delta is >= 1/n', {
	data(PUMS5extract10000, package = "PSIlence")
	
	my_n.bins <- 16
	my_n <- 10000
	my_epsilon <- 0.1
	my_delta <- 0.1 # set delta to > 1/n
	
	dp.histogram2 <- dpHistogram$new(var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, delta=my_delta)
	expect_error(dp.histogram2$release(PUMS5extract10000), "Delta must be less than 1/n")
})

test_that('histogram on categorical data', {
	library(datasets)
	data(esoph)
	
	my_n <- 88
	my_epsilon <- 1
	my_delta <- 10^-3
	
	catHistogram <- dpHistogram(var.type='character', variable='tobgp', n=my_n, epsilon=my_epsilon, delta=my_delta)
	catHistogram$release(esoph)
	
	askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismStability', n.bins=length(catHistogram$result$release), n=my_n, epsilon=my_epsilon, delta=my_delta)
	expect_equal(catHistogram$epsilon, my_epsilon)
	expect_equal(catHistogram$accuracy, askAccuracy)
})

test_that('histogram release has expected dimensions and accuracy for manually created logical variable with impute = true (laplace mechanism)', {
	logicalVar_withNA <- c(1,0,1,1,1,0,1,0,0,NA,1,0,NA,1,0,0,1,NA,NA,1,1,0,1,0,1,0)
	dataLog <- data.frame(logicalVar_withNA)
	
	my_n <- 26
	my_epsilon <- 1
	my_delta <- 10^-2
	
	dp.histogram <- dpHistogram$new(var.type='logical', variable="logicalVar_withNA", n=my_n, epsilon=my_epsilon, delta=my_delta, impute = TRUE)
	dp.histogram$release(dataLog)
	expect_equal(length(dp.histogram$result$release), 2) # there should be 3 bins when impute = FALSE: 0,1,NA
	
	askAccuracy <- histogram.getAccuracy(mechanism = 'mechanismLaplace', n.bins=my_n.bins, n=my_n, epsilon=my_epsilon, delta=my_delta)
	expect_equal(dp.histogram$epsilon, my_epsilon)
	expect_equal(dim(dp.histogram$result$interval), c(2,2))
	expect_equal(dp.histogram$accuracy, askAccuracy)
})
