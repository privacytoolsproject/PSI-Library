library(PSIlence)
context("utilities")

## checkrange tests##
rng1 = c(0,1)
rng2 = c(0,1,2)
rng3 = c(1)

expect_equal(checkrange(rng1),rng1)
expect_warning(checkrange(rng2),"range argument supplied has more than two values.  Will proceed using min and max values as range.")
expect_equal(checkrange(rng2), c(0,2))
expect_error(checkrange(rng3),"range argument in error: requires upper and lower values as vector of length 2.")

## censordata tests ##

# create test data
residence <- factor(c("WA", "OR", "OR", "OR", "WA","CA"))
chars = c('a', 'b', 'c', 'c', 'd')
nums = 1:10

residence_out = censordata(x=residence, var_type='factor', levels=c("OR"))
chars_out = censordata(x=chars, var_type='character', levels=c('a', 'b', 'c'))
nums_out = censordata(x=nums, var_type='integer', rng=c(2.5, 7))

expect_equal(residence_out, factor(c(NA,"OR","OR","OR",NA,NA)))
expect_equal(chars_out, factor(c('a','b','c','c',NA)))
expect_equal(nums_out, c(2.5,2.5,3.0,4.0,5.0,6.0,7.0,7.0,7.0,7.0))

## getFuncArgs tests ##

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
