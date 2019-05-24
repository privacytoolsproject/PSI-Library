library(PSIlence)
context("histogram")

test_that('histogram getAccuracy and getParameters return approximately correct values', {
  val1 <- round(histogram.getAccuracy(n.bins=10, n=2000, epsilon=0.2, stability=TRUE, delta=10^-6))
  val2 <- round(histogram.getParameters(n.bins=20, n=2000, accuracy=0.5, stability=TRUE, delta=10^-6))
  expect_equal(val1, 175)
  expect_equal(val2, 70)
})

test_that('histogram releases have expected dimensions', {
  data(PUMS5extract10000, package = "PSIlence")
  
  my_n.bins <- 16
  my_n <- 10000
  my_epsilon <- 0.1
  my_delta <- 10^-6
  
  dp.histogram <- dpHistogram$new(mechanism='mechanismLaplace', var.type='numeric', variable="educ", n=my_n, epsilon=my_epsilon, n.bins=my_n.bins, delta=my_delta, rng=c(0.5,16.5))
  dp.histogram$release(PUMS5extract10000)
  expect_equal(length(dp.histogram$result$release), 16)
  
  askAccuracy <- histogram.getAccuracy(n.bins=my_n.bins, n=my_n, epsilon=my_epsilon, stability=FALSE, delta=my_delta)
  expect_equal(dp.histogram$epsilon, my_epsilon)
  expect_equal(dim(dp.histogram$result$interval), c(16,2))
  expect_equal(dp.histogram$accuracy, askAccuracy)
})
