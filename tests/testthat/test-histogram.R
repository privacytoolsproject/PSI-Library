test_that('histogram getAccuracy and getParameters return approximately correct values', {
  library(PSIlence)
  expect_equal(round(histogram.getAccuracy(n.bins=10, n=2000, epsilon=0.2, stability=TRUE, delta=10^-6)*100), 28)
  expect_equal(histogram.getParameters(n.bins=20, n=2000, accuracy=0.5, stability=TRUE, delta=10^-6)*100), 11)
})
