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

test_that('checkRange is as expected', {
  rng1 = c(0,1)
  rng2 = c(0,1,2)
  rng3 = c(1)
  
  expect_equal(checkRange(rng1, "numeric"),rng1)
  expect_warning(checkRange(rng2, "numeric"),"range argument supplied has more than two values.  Will proceed using min and max values as range.")
  expect_warning(expect_equal(checkRange(rng2, "numeric"), c(0,2)))
  expect_error(checkRange(rng3, "numeric"),"range argument in error: requires upper and lower values as vector of length 2.")
})

test_that('censorData is as expected', {
  residence <- factor(c("WA", "OR", "OR", "OR", "WA","CA"))
  chars = c('a', 'b', 'c', 'c', 'd')
  nums = 1:10
  
  residenceOut = censorData(x=residence, varType='factor', levels=c("OR"))
  charsOut = censorData(x=chars, varType='character', levels=c('a', 'b', 'c'))
  numsOut = censorData(x=nums, varType='integer', rng=c(2.5, 7))
  
  expect_equal(residenceOut, factor(c(NA,"OR","OR","OR",NA,NA)))
  expect_equal(charsOut, factor(c('a','b','c','c',NA)))
  expect_equal(numsOut, c(2.5,2.5,3.0,4.0,5.0,6.0,7.0,7.0,7.0,7.0))
})