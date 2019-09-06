library(PSIlence)
context("utilities")

test_that('Uniform distribution drawing functions as expected', {
  uniformRepeatable1 <- dpUnif(n=1, seed=75436)
  uniformRepeatable2 <- dpUnif(n=1, seed=75436)
  expect_equal(uniformRepeatable1, uniformRepeatable2)
  
  n <- 5
  uniformList <- dpUnif(n=n)
  expect_length(uniformList, n)
  expect_equal(sum(uniformList<1), n)
  expect_equal(sum(uniformList>0), n)
})

test_that('checkRange is as expected', {
  rng1 = c(0,1)
  rng2 = c(0,1,2)
  rng3 = c(1)
  
  expect_equal(checkRange(rng1, "numeric"),rng1)
  expect_warning(checkRange(rng2, "numeric"),"range argument supplied has more than two values.  Will proceed using min and max values as range.")
  expect_error(checkRange(rng3, "numeric"),"range argument in error: requires upper and lower values as vector of length 2.")
  # make sure that after the warning is thrown for a range argument is too long, the correct value is returned
  rangeTest <- expect_warning(checkRange(rng2, "numeric"),"range argument supplied has more than two values.  Will proceed using min and max values as range.")
  expect_equal(rangeTest, c(0,2))
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
  yImputed <- fillMissing(x=y, varType='numeric', imputeRng=c(-1,1))
  
  expect_equal(sum(is.na(yImputed)), 0)
  
  s <- sample(animals, size=100, replace=TRUE)
  s[sample(1:100, size=10)] <- NA
  sImputed <- fillMissing(x=s, varType='factor', categories=animals)
  
  expect_equal(sum(is.na(sImputed)), 0)
  expect_true(is.factor(sImputed))
  
  N <- 100
  x1 <- x2 <- rnorm(N)
  x1[sample(1:N, size=10)] <- NA
  x2[sample(1:N, size=10)] <- NA
  impRng <- matrix(c(-3, 3, -2, 2), ncol=2, byrow=TRUE)
  df <- data.frame(x1, x2)
  dfImputed <- fillMissing(x=df, varType='numeric', imputeRng=impRng)
  
  expect_equal(sum(is.na(dfImputed)), 0)
})


# testing checkDelta() utility function
test_that('testing checkDelta() function', {
    expect_warning(checkDelta('mechanismLaplace', 10^-5), 'A delta parameter has been entered, but a mechanism that uses a delta value is not being used. Setting delta to 0.')

    delta2 <- checkDelta('mechanismLaplace')
    expect_equal(delta2, 0)

    delta3 <- checkDelta('mechanismExponential')
    expect_equal(delta3, 2^-30)

    delta4 <- checkDelta('mechanismStability')
    expect_equal(delta4, 2^-30)

    delta5 <- checkDelta('mechanismGaussian')
    expect_equal(delta5, 2^-30)

    delta6 <- checkDelta('mechanismExponential', 10^-5)
    expect_equal(delta6, 10^-5)

    delta7 <- checkDelta('mechanismStability', 10^-10)
    expect_equal(delta7, 10^-10)

    delta8 <- checkDelta('mechanismGaussian', 2^-15)
    expect_equal(delta8, 2^-15)
})
