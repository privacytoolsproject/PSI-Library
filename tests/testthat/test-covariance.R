library(PSIlence)
context("covariance")

data(PUMS5extract10000)

######################### Input validation tests #######################################################

test_that('epsilon checks throw correct warning', {
    rangeIncome <- c(-10000, 713000)
    rangeEducation <- c(1, 16)
    range <- rbind(rangeIncome, rangeEducation)
    
    expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000, 
                                  epsilon = -0.1, columns = c("income", "education"), rng = range),
                 "Privacy parameter epsilon must be a value greater than zero.")
})

test_that('range checks throw correct warning', {
    rng <- c(100)
    rngStr <- paste('c(',toString(rng),')')
    errorStr <- paste('Error in range argument provided,', rngStr, ': requires upper and lower values as vector of length 2.')
    
    #fixed=TRUE forces error string to be interpreted as fixed string rather than regular expression, which is necessary due to parenthesis in error message.
    expect_error(dpCovariance$new(mechanism='mechanismLaplace', varType='numeric', n= 10000,
                                  epsilon = 0.1, columns = c("income", "education"),
                                  rng=rng),
                 errorStr, fixed=TRUE) 
    
    rng <- matrix(c(-10,0,100), nrow=1)
    rngStr <- paste('c(', toString(rng), ')')
    warningStr <- paste('Range argument of', rngStr, 'has more than two values.  Will proceed using min and max values as range.')
    
    #fixed=TRUE forces warning string to be interpreted as fixed string rather than regular expression, which is necessary due to parenthesis in error message.
    expect_warning(dpCovariance$new(mechanism='mechanismLaplace', varType='numeric', n=10000,
                                    epsilon=0.1, columns = c("income", "education"),
                                    rng=rng),
                   warningStr, fixed=TRUE)
})

# make sure error thrown when n not positive or a whole number
test_that('make sure error thrown when n not positive or a whole number',{
    expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = -1, 
                                  epsilon = 1, columns = c("income", "educ"), rng = range, formula='x~y'),
                 "n must be a positive whole number")
    expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 0.5, 
                                  epsilon = 1, columns = c("income", "educ"), rng = range, formula='x~y'),
                 "n must be a positive whole number")
})

######################### Workflow tests #######################################################

test_that('DP covariance workflow runs', {
    range.income <- range(PUMS5extract10000['income'])
    range.education <- range(PUMS5extract10000['educ'])
    range.age <- range(PUMS5extract10000['age'])
    range <- rbind(range.income, range.education, range.age)
    
    dpCov <- dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000,
                              epsilon = 1, columns = c("income", "educ", 'age'), rng = range, formula='income~educ')
    out <- dpCov$release(PUMS5extract10000)
    expect_equal(length(out),3)
})

test_that('coefficient release function operational in workflow', {
    range.income <- range(PUMS5extract10000['income'])
    range.education <- range(PUMS5extract10000['educ'])
    range.age <- range(PUMS5extract10000['age'])
    range <- rbind(range.income, range.education, range.age)
    
    #Next line expected to throw warning due to high epsilon val.
    dpCov <- dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000,
                              epsilon = 10000000000, columns = c("income", "educ", "age"), rng = range, formula='income~educ')
    out <- dpCov$release(PUMS5extract10000)
    coeffs <- coefficientRelease('income~age', out$release, n=10000)
    expect_equal(length(coeffs), 4)
    linreg <- lm(income~age, data=PUMS5extract10000)
    
    output <- as.numeric(coeffs$coefficients[[1]][1]) #extracts coefficient from output
    expectedOutput <- as.numeric(linreg[[1]][2])
    expect_equal(floor(output), floor(expectedOutput)) #check floor of values due to fact that there is some noise added here
})
