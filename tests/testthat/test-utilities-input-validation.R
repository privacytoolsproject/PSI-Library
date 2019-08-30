context('utilities-input-validation')

test_that('checkNumeric raises proper warnings and errors', {
  
})

test_that('checkN raises proper warnings and errors',{
  expect_error(checkN(-1))
  expect_error(checkN(1.5))
  expect_error(checkN(-1.5))
  expect_error(checkN(c(1,2)))
  expect_error(checkN('foo'))
})

test_that('checkRange raises proper warnings and errors on 1D input', {
  rng1 <- c(0,1)
  rng2 <- c(0,1,2)
  rng3 <- c(1)
  
  expect_equal(checkRange(rng1, "numeric"),rng1)
  
  rng2Str <- paste('c(',toString(rng2),')')
  expectedWarning  <- paste('Range argument of', rng2Str, 'has more than two values.  Will proceed using min and max values as range.')
  expect_warning(checkRange(rng2, "numeric"), expectedWarning, fixed=TRUE)
  expect_warning(expect_equal(checkRange(rng2, "numeric"), c(0,2)))
  
  rng3Str <- paste('c(',toString(rng3),')')
  expectedError <- paste('Error in range argument provided,', rng3Str, ': requires upper and lower values as vector of length 2.')
  expect_error(checkRange(rng3, "numeric"), expectedError, fixed=TRUE)
})

test_that('checkRange raises proper warnings and errors 2D input', {
  rng <- matrix(c(1,2,3,4,5,6), nrow = 3)
  expect_equal(rng, checkRange(rng, 'numeric'))
  
  rng <- matrix(c(1,2,3,4,5,6), nrow = 2)
  outRng <- matrix(c(1,2,5,6), nrow=2) 
  expect_warning(expect_equal(checkRange(rng, 'numeric'),outRng))
})

test_that('checkEpsilon raises proper warnings and errors', {
  expect_equal(checkEpsilon(1), 1)
  expect_equal(checkEpsilon(c(1,0.1), multipleEps=TRUE, expectedLength=2), c(1,0.1))
  
  expect_error(checkEpsilon(c(1,0.1)))
  expect_error(checkEpsilon(-1))
  expect_error(checkEpsilon(c(1,-0.1), multipleEps=TRUE))
  expect_error(checkEpsilon(c(1,2), multupleEps=TRUE, expectedLength=5))
  expect_error(checkEpsilon('foo'))
  
  expect_warning(checkEpsilon(4))
  expect_warning(checkEpsilon(c(1,4), multipleEps=TRUE, expectedLength=2))
})

test_that('checkAccuracy raises proper warnings and errors', {
  expect_error(checkAccuracy(-1))
  expect_equal(checkAccuracy(5),5)
  expect_error(checkAccuracy('foo'))
})