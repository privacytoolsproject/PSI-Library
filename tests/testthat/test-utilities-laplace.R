library(PSIlence)
context('utilities-Laplace')

test_that('Accuracy derivation works for inputs of varying dimensions',{
  expOut <- c(laplaceGetAccuracy(1,2),laplaceGetAccuracy(2,4))
  trueOut <- laplaceGetAccuracy(c(1,2), c(2,4))
  expect_equal(expOut, trueOut)
  
  expOut <- c(laplaceGetAccuracy(2,3), laplaceGetAccuracy(2,4))
  trueOut <- laplaceGetAccuracy(2, c(3,4))
  expect_equal(expOut, trueOut)
  
  expOut <- c(laplaceGetAccuracy(3,2), laplaceGetAccuracy(4,2))
  trueOut <- laplaceGetAccuracy(c(3,4),2)
  expect_equal(expOut, trueOut)
})

test_that('Epsilon derivation works for inputs of varying dimensions', {
  expOut <- c(laplaceGetEpsilon(1,2),laplaceGetEpsilon(2,4))
  trueOut <- laplaceGetEpsilon(c(1,2), c(2,4))
  expect_equal(expOut, trueOut)
  
  expOut <- c(laplaceGetEpsilon(2,3), laplaceGetEpsilon(2,4))
  trueOut <- laplaceGetEpsilon(2, c(3,4))
  expect_equal(expOut, trueOut)
  
  expOut <- c(laplaceGetEpsilon(3,2), laplaceGetEpsilon(4,2))
  trueOut <- laplaceGetEpsilon(c(3,4),2)
  expect_equal(expOut, trueOut)
})