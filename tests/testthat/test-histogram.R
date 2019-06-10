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
  # print(dp.histogram$result$release)
  # print(table(PUMS5extract10000$educ))
  # print(dp.histogram$result)
  
  my_n.bins2 <- 2
  my_n2 <- 10000
  my_epsilon2 <- 0.1
  my_delta2 <- 10^-6
  
  dp.histogram2 <- dpHistogram$new(mechanism='mechanismLaplace', var.type='logical', variable="sex", n=my_n2, epsilon=my_epsilon2, n.bins=my_n.bins2, delta=my_delta2, impute = TRUE, rng = c(0,1))
  dp.histogram2$release(PUMS5extract10000)
  expect_equal(length(dp.histogram2$result$release), 2)
  # print(dp.histogram2$result$release)
  
  my_n.bins3 <- 76
  my_n3 <- 10000
  my_epsilon3 <- 0.1
  my_delta3 <- 10^-6
  
  dp.histogram3 <- dpHistogram$new(mechanism='mechanismLaplace', var.type='numeric', variable="age", n=my_n3, epsilon=my_epsilon3, n.bins=my_n.bins3, delta=my_delta3, rng=c(17.5,93.5))
  dp.histogram3$release(PUMS5extract10000)
  expect_equal(length(dp.histogram3$result$release), 76)
  # print(dp.histogram$result$release)
  print(table(cut(sort(PUMS5extract10000$age), breaks = 76)))
  print(dp.histogram3$result$release)
  print(dp.histogram3$stability)
  
  askAccuracy <- histogram.getAccuracy(n.bins=my_n.bins, n=my_n, epsilon=my_epsilon, stability=FALSE, delta=my_delta)
  expect_equal(dp.histogram$epsilon, my_epsilon)
  expect_equal(dim(dp.histogram$result$interval), c(16,2))
  expect_equal(dp.histogram$accuracy, askAccuracy)
})
