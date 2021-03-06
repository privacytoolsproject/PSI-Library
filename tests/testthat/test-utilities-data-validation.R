library(PSIlence)
context('utilities-data-validation')

test_that('censorData runs as expected', {
    residence <- factor(c("WA", "OR", "OR", "OR", "WA","CA"))
    chars <- c('a', 'b', 'c', 'c', 'd')
    nums <- 1:10
    m <- matrix(1:10, ncol=2)
    
    residenceOut <- censorData(residence, varType='factor', levels=c("OR"))
    charsOut <- censorData(chars, varType='character', levels=c('a', 'b', 'c'))
    numsOut <- censorData(nums, varType='integer', rng=c(2.5, 7), rngFormat="vector")
    mOut <- censorData(m, varType='integer', rng=list(c(1,5),c(7,9)), rngFormat="list")
    
    expect_equal(residenceOut, factor(c(NA,"OR","OR","OR",NA,NA)))
    expect_equal(charsOut, factor(c('a','b','c','c',NA)))
    expect_equal(numsOut, c(2.5,2.5,3.0,4.0,5.0,6.0,7.0,7.0,7.0,7.0))
    expect_equal(mOut, matrix(c(1:5, 7,7,8,9,9), ncol=2))
})

test_that('censorData throws expected errors', {
    nums <- 1:10
    m <- matrix(11:20, ncol=2)
    
    expect_error(censorData(nums, 'integer'))
    expect_error(censorData(nums, 'integer', rng=c(1,10), rngFormat='list'))
    expect_error(censorData(nums, 'integer', rng=list(c(1,10),c(1,5)), rngFormat='list'))
    
    expect_error(censorData(m, 'integer', rng=c(1,10), rngFormat='vector'))
    expect_error(censorData(m, 'integer', rng=list(c(1,10),c(1,5)), rngFormat='vector'))
})

test_that('fillMissing', {
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