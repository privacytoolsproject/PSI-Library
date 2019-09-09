context('utilities-input-validation')

test_that('checkEmpty raises proper errors', {
  expect_equal(checkEmpty(15.5), 15.5)
  expect_equal(checkEmpty(c(1,4)), c(1,4))
  expect_equal(checkEmpty(0), 0)
  expect_equal(checkEmpty(NULL, emptyOkay=TRUE), NULL) 
  expect_equal(checkEmpty(NA, emptyOkay=TRUE), NA)
  
  expect_error(checkEmpty(NULL))
  expect_error(checkEmpty(NA))
})

test_that('checkNumeric raises proper errors', {
  expect_error(checkNumeric('foo'))
  expect_equal(checkNumeric(15.5), 15.5)
  expect_equal(checkNumeric(c(15.5, 17)), c(15.5, 17))
  
  expect_error(checkNumeric(NA))
  expect_error(checkNumeric(NULL))
  expect_error(checkNumeric(c(1, NA)))
  
  expect_equal(checkNumeric(NULL, emptyOkay=TRUE), NULL)
  expect_equal(checkNumeric(NA, emptyOkay=TRUE), NA)
  expect_equal(checkNumeric(c(2, NA), emptyOkay=TRUE), c(2,NA))
})

test_that('checkLength raises proper warnings and errors', {
  xs <- c(1,2,3)
  expect_equal(checkLength(xs,3), xs)
  expect_error(checkLength(xs,5))
  
  li <- list(1,"a", c(1,2))
  expect_equal(checkLength(li, 3), li)
  expect_error(checkLength(li, 5))
})

test_that('checkN raises proper warnings and errors',{
  expect_equal(checkN(8),8)
  expect_error(checkN(-1))
  expect_error(checkN(1.5))
  expect_error(checkN(-1.5))
  expect_error(checkN(c(1,2)))
  expect_error(checkN('foo'))

  expect_equal(checkN(c(1,2,3), expectedLength=3), c(1,2,3))
  expect_error(checkN(1, expectedLength=4))
  expect_error(checkN(c(1,-1,3), expectedLength=3))
  expect_error(checkN(c(1, 1.5, 4)))
  
  expect_equal(checkN(NA, emptyOkay=TRUE), NA)
  expect_equal(checkN(NULL, emptyOkay=TRUE), NULL)
  expect_error(checkN(NA), "Input n may not be NA or NULL.")
  expect_error(checkN(NULL, expectedLength=0), "Input n may not be NA or NULL.")
  
  expect_equal(checkN(c(1,NA,2), expectedLength=3, emptyOkay=TRUE), c(1,NA,2))
  expect_error(checkN(c(1,NA,2), expectedLength=3))
  expect_error(checkN(c(1,NA, 2.4), expectedLength=3))
})

test_that('checkEpsilon raises proper warnings and errors', {
  expect_equal(checkEpsilon(1), 1)
  expect_equal(checkEpsilon(c(1,0.1), expectedLength=2), c(1,0.1))
  
  expect_error(checkEpsilon(c(1,0.1)))
  expect_error(checkEpsilon(-1))
  expect_error(checkEpsilon(c(1,-0.1)))
  expect_error(checkEpsilon(c(1,2), expectedLength=5))
  expect_error(checkEpsilon('foo'))
  
  expect_warning(checkEpsilon(4))
  expect_warning(checkEpsilon(c(1,4), expectedLength=2))
  
  expect_error(checkEpsilon(NULL))
  expect_error(checkEpsilon(NA))
  expect_error(checkEpsilon(c(1,NA), expectedLength=2))
})

test_that('checkAccuracy raises proper warnings and errors', {
  expect_error(checkAccuracy(-1))
  expect_equal(checkAccuracy(5),5)
  expect_error(checkAccuracy('foo'))
  
  expect_equal(checkAccuracy(c(1,5), expectedLength=2), c(1,5))
  expect_error(checkAccuracy(c(1,5)))
})

test_that('checkRange raises proper warnings and errors on 1D input', {
  rng1 <- c(0,1)
  rng2 <- c(0,1,2)
  rng3 <- c(1)
  
  expect_equal(checkRange(rng1, "numeric", "vector"), rng1)
  expect_equal(checkRange(rng1, "numeric", "vector", emptyOkay=TRUE), rng1)
  
  rng2Str <- paste('c(',toString(rng2),')')
  expectedWarning  <- paste('Range argument of', rng2Str, 'has more than two values.  Will proceed using min and max values as range.')
  expect_warning(checkRange(rng2, "numeric", "vector"), expectedWarning, fixed=TRUE)
  expect_warning(expect_equal(checkRange(rng2, "numeric", "vector"), c(0,2)))
  
  rng3Str <- paste('c(',toString(rng3),')')
  expectedError <- paste('Error in range argument provided,', rng3Str, ': requires upper and lower values as vector of length 2.')
  expect_error(checkRange(rng3, "numeric", "vector"), expectedError, fixed=TRUE)
  
  expect_error(checkRange(NULL, "numeric", "vector"), "Input range may not be empty.")
  expect_error(checkRange(c(NA,1), "numeric", "vector"), "Input range may not contain NA values.")
  expect_error(checkRange(c(NULL, 1), "numeric", "vector"), "Error in range argument provided, c( 1 ) : requires upper and lower values as vector of length 2.", fixed=TRUE)
  
  expect_equal(checkRange(NULL, "numeric", "vector", emptyOkay=TRUE), NULL)
  expect_equal(checkRange(NA, "numeric", "vector", emptyOkay=TRUE), NULL)
  
  expect_warning(expect_equal(checkRange(c(NA,1), "numeric", "vector", emptyOkay=TRUE), NULL), "Range argument provided c( NA, 1 ) has NA value. Setting range to NULL.", fixed=TRUE)
  
  expect_equal(checkRange(NULL, "logical", "vector"), c(0,1))
  expect_equal(checkRange(5, "logical", "vector"), c(0,1))
})


test_that('checkRange raises proper warnings and errors on list and matrix inputs', {
  #skip("need to add expectedLength")
  li <- list(c(1,4),c(2,5),c(3,6))
  m <- matrix(c(1,2,3,4,5,6), nrow = 3)
  expect_equal(checkRange(li, 'numeric', 'list', 3), li)
  expect_equal(checkRange(m, 'numeric', 'matrix', 3), li)
  
  rng <- matrix(c(1,2,3,4,5,6), nrow = 2)
  outRng <- list(c(1,5),c(2,6)) 
  expect_warning(expect_equal(checkRange(rng, 'numeric', 'matrix', 2), outRng))
  
  rng <- matrix(rep(NA,3,3),3)
  expect_equal(checkRange(rng, "numeric",'matrix', 3, emptyOkay=TRUE), list(NULL, NULL, NULL))
  expect_error(checkRange(rng, "numeric", 'matrix', 3), "Input range may not be empty.")
  
  rng <- list(c(1,2), c(2,NA))
  expect_equal(expect_warning(checkRange(rng, "numeric", 'list', 2, emptyOkay=TRUE)), list(c(1,2), NULL))
  expect_equal(checkRange(rng, "logical", 'list', 2), list(c(0,1), c(0,1)))
  expect_error(checkRange(rng, "numeric",'list', 2))
  
  rng <- list(c(1,2), NA)
  expect_equal(checkRange(rng, "numeric",'list', 2, emptyOkay=TRUE), list(c(1,2), NULL))
  expect_equal(checkRange(rng, "logical",'list', 2), list(c(0,1), c(0,1)))
  expect_error(checkRange(rng, "numeric",'list', 2))
  
  rng <- list(c(1,2,3), 2)
  expect_error(expect_warning(checkRange(rng, "numeric",'list', 2)))
  
  rng <- list(c(1,3))
  expect_equal(checkRange(rng, "numeric", 'list',1), rng)
  expect_equal(checkRange(rng, "numeric",'list', 1, emptyOkay=TRUE), rng)
  
  rng <- list(NULL, c(5,6))
  expect_equal(checkRange(rng, "numeric",'list', 2, emptyOkay=TRUE), rng)
  expect_error(checkRange(rng, "numeric",'list'))
})


test_that('checkVariableType raises proper warnings and errors', {
  expect_equal(checkVariableType('Numeric', c('Numeric', 'Factor')), 'numeric')
  expect_equal(checkVariableType('numeric', c('Numeric', 'Integer')), 'numeric')
  
  expect_error(checkVariableType('Numeric', c('Factor', 'Integer')))
  
  expect_equal(checkVariableType(NULL, c('Factor', 'Integer'), emptyOkay=TRUE), NULL)
  expect_equal(checkVariableType(NULL, 'Factor', emptyOkay=TRUE), NULL)
  expect_error(checkVariableType(NULL, c('Factor', 'Integer')))
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

test_that('checkMechanism raises proper warnings and errors', {
  ls <- c('foo', 'bar')
  expect_equal(checkMechanism('foo', c('Foo', 'Bar')), 'Foo')
  expect_error(checkMechanism(NA, ls), "Input may not be NA or NULL.")
  expect_error(checkMechanism(NULL, ls))
  expect_error(checkMechanism(ls, ls))
  expect_error(checkMechanism('fooo', ls))
  
})

test_that("Matrix to list conversion works", {
  m <- matrix(c(1,2,3,4,5,6), nrow = 3)
  li <- list(c(1,4), c(2,5), c(3,6))
  
  expect_equal(matrixToList(m), li)
  
  m <- matrix(rep(NA,3),3)
  li <- list(NA, NA, NA)
  expect_equal(matrixToList(m), li)
  
  m <- matrix(c(1,NA,2,NA,3,4), nrow=3)
  li <- list(c(1,NA), c(NA,3), c(2,4))
  expect_equal(matrixToList(m), li)
})

test_that("isVector runs correctly", {
  expect_true(isVector(c(1,2,3)))
  expect_true(isVector(1))
  expect_false(isVector(list(1,2)))
  expect_false(isVector(matrix(c(1,2,3,4))))
})