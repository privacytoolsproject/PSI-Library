context('utilities-input-validation')

test_that('checkNumeric raises proper warnings and errors', {
  expect_error(checkNumeric('foo'))
  expect_equal(checkNumeric(15.5), 15.5)
})

test_that('checkN raises proper warnings and errors',{
  expect_error(checkN(-1))
  expect_error(checkN(1.5))
  expect_error(checkN(-1.5))
  expect_error(checkN(c(1,2)))
  expect_error(checkN('foo'))
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

test_that('checkVariableType raises proper warnings and errors', {
  expect_equal(checkVariableType('Numeric', c('Numeric', 'Factor')), 'numeric')
  expect_equal(checkVariableType('numeric', c('Numeric', 'Integer')), 'numeric')
  
  expect_error(checkVariableType('Numeric', c('Factor', 'Integer')))
})

test_that('checkImputationRange raises proper warnings and errors', {
  expect_equal(checkImputationRange(c(1,2),c(0,3),'numeric'), c(1,2))
  expect_equal(checkImputationRange(c(1,2),c(0,3), 'integer'), c(1,2))
  expect_equal(checkImputationRange(c(1,2), c(1,2), 'numeric'), c(1,2))
  
  #To Do: Discuss this case with Megan
  expect_warning(expect_equal(checkImputationRange(c('cat','dog'),c('cat','dog','mouse'), 'character'), 
                              c('cat','dog','mouse')))
  
  expect_error(checkImputationRange(c('cat','dog'),c(1,5),'numeric'))
  expect_error(checkImputationRange(1,c(1,5),'numeric'))
  expect_error(checkImputationRange(c(1,2,5), 'numeric'))
  
  expect_warning(expect_equal(checkImputationRange(c(0,3),c(1,2),'numeric'), c(1,2)))
  expect_warning(expect_equal(checkImputationRange(c(1.5,3), c(1,2),'numeric'), c(1.5,2)))
  expect_warning(expect_equal(checkImputationRange(c(0.5,1.5),c(1,2), 'numeric'),c(1,1.5)))
  
})