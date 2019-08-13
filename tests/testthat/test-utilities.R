library(PSIlence)
context("utilities")

test_that('Uniform distribution drawing functions as expected', {
  uniform_repeatable1 <- dpUnif(n=1, seed=75436)
  uniform_repeatable2 <- dpUnif(n=1, seed=75436)
  expect_equal(uniform_repeatable1, uniform_repeatable2)
  
  n <- 5
  uniform_list <- dpUnif(n=n)
  expect_length(uniform_list, n)
  expect_equal(sum(uniform_list<1), n)
  expect_equal(sum(uniform_list>0), n)
})

test_that('checkRange is as expected', {
  rng1 = c(0,1)
  rng2 = c(0,1,2)
  rng3 = c(1)
  
  expect_equal(checkRange(rng1, "numeric"),rng1)
  expect_warning(checkRange(rng2, "numeric"),"range argument supplied has more than two values.  Will proceed using min and max values as range.")
  expect_equal(checkRange(rng2, "numeric"), c(0,2))
  expect_error(checkRange(rng3, "numeric"),"range argument in error: requires upper and lower values as vector of length 2.")
})

test_that('censorData is as expected', {
  residence <- factor(c("WA", "OR", "OR", "OR", "WA","CA"))
  chars = c('a', 'b', 'c', 'c', 'd')
  nums = 1:10
  
  residence_out = censorData(x=residence, varType='factor', levels=c("OR"))
  chars_out = censorData(x=chars, varType='character', levels=c('a', 'b', 'c'))
  nums_out = censorData(x=nums, varType='integer', rng=c(2.5, 7))
  
  expect_equal(residence_out, factor(c(NA,"OR","OR","OR",NA,NA)))
  expect_equal(chars_out, factor(c('a','b','c','c',NA)))
  expect_equal(nums_out, c(2.5,2.5,3.0,4.0,5.0,6.0,7.0,7.0,7.0,7.0))
})

test_that('getFuncArgs performs as expected', {
  foo = function(thing1, var1){
    return('bar')
  }
  
  testList <- list(var1 = 1, var2 = 2, var3 = 3)
  
  testClass <- setRefClass(
    Class = 'test',
    fields = list(
      thing1 = 'numeric',
      thing2 = 'numeric',
      thing3 = 'numeric'
    )
  )
  
  testClass$methods(
    getFields = function() {
      f <- names(getRefClass()$fields())
      out <- setNames(vector('list', length(f)), f)
      for (fd in f) {
        out[[fd]] <- .self[[fd]]
      }
      return(out)
    })
  
  testObject <- testClass$new(thing1=1, thing2=2, thing3=3)
  
  expect_equal(getFuncArgs(foo, inputList=testList, inputObject=testObject), list(thing1=1,var1=1))
  expect_equal(getFuncArgs(foo, inputList=testList), list(var1=1))
  expect_equal(getFuncArgs(foo, inputObject=testObject), list(thing1=1))
  expect_equal(getFuncArgs(foo), list())
})


############################# fillMissing tests #######################################

test_that('fillMissing as expected', {
  x <- dpUnif(10)
  animals <- as.factor(c('dog', 'zebra', 'bird', 'hippo'))
  lowBound <- 1
  upBound <- 5
  
  scaledNumeric <- scaleValues(x, 'numeric', lower=lowBound, upper=upBound)
  scaledLogical <- scaleValues(x, 'logical')
  scaledInteger <- scaleValues(x, 'integer', lower=lowBound, upper=upBound)
  scaledFactor <- scaleValues(x, 'factor', categories=animals)
  
  expect_equal(length(scaledNumeric), length(x))
  expect_equal(length(scaledLogical), length(x))
  expect_equal(length(scaledInteger), length(x))
  expect_equal(length(scaledFactor), length(x))
  
  expect_equal(sum(scaledNumeric < lowBound), 0)
  expect_equal(sum(scaledNumeric > upBound), 0)
  expect_equal(sum(scaledInteger < lowBound), 0)
  expect_equal(sum(scaledInteger > upBound), 0)
  expect_equal(sum(scaledLogical < 0), 0)
  expect_equal(sum(scaledLogical > 1), 0)
  expect_equal(sum((scaledLogical > 0) & (scaledLogical < 1)), 0)
  
  y <- rnorm(100)
  y[sample(1:100, size=10)] <- NA
  y_imputed <- fillMissing(x=y, varType='numeric', imputeRng=c(-1,1))
  
  expect_equal(sum(is.na(y_imputed)), 0)
  
  s <- sample(animals, size=100, replace=TRUE)
  s[sample(1:100, size=10)] <- NA
  s_imputed <- fillMissing(x=s, varType='factor', categories=animals)
  
  expect_equal(sum(is.na(s_imputed)), 0)
  expect_true(is.factor(s_imputed))
  
  N <- 100
  x1 <- x2 <- rnorm(N)
  x1[sample(1:N, size=10)] <- NA
  x2[sample(1:N, size=10)] <- NA
  impRng <- matrix(c(-3, 3, -2, 2), ncol=2, byrow=TRUE)
  df <- data.frame(x1, x2)
  df_imputed <- fillMissing(x=df, varType='numeric', imputeRng=impRng)
  
  expect_equal(sum(is.na(df_imputed)), 0)
})