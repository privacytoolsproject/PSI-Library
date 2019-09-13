library(PSIlence)
context("covariance")

data(PUMS5extract10000)

######################### Input validation tests #######################################################

test_that('epsilon checks throw correct warning', {
    rangeIncome <- c(-10000, 713000)
    rangeEducation <- c(1, 16)
    range <- list(rangeIncome, rangeEducation)

    expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000,
                     epsilon = -0.1, columns = c("income", "education"), rng = range),
                 "Privacy parameter epsilon must be a value greater than zero.")
})

test_that('range checks throw correct error', {
    rng <- list(c(100))

    expect_error(dpCovariance$new(mechanism='mechanismLaplace', varType='numeric', n= 10000,
                                  epsilon = 0.1, columns = c("income", "education"),
                                  rng=rng),
                 "Input was expected to be of length  2  but is instead of length  1")

    rng <- matrix(c(-10,0,100), nrow=1)

    expect_error(dpCovariance$new(mechanism='mechanismLaplace', varType='numeric', n=10000,
                                    epsilon=0.1, columns = c("income", "education"),
                                    rng=rng),
                   "Input was expected to be of length  2  but is instead of length  1")
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
    range <- list(range.income, range.education, range.age)
    
    dpCov <- dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000,
                              epsilon = c(1,1,1), columns = c("income", "educ", 'age'), rng = range, formula='income~educ')
    out <- dpCov$release(PUMS5extract10000)
    expect_equal(length(out),3)
})

test_that('coefficient release function operational in workflow', {

    range.income <- range(PUMS5extract10000['income'])
    range.education <- range(PUMS5extract10000['educ'])
    range.age <- range(PUMS5extract10000['age'])
    range <- list(range.income, range.education, range.age)
    
    eps <- c(rep(10000000000,3))

    #Next line expected to throw warning due to high epsilon val.
    dpCov <- expect_warning(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000,
                              epsilon = eps, columns = c("income", "educ", "age"), rng = range, formula='income~educ'))
    out <- dpCov$release(PUMS5extract10000)
    coeffs <- coefficientRelease('income~age', out$release, n=10000)
    expect_equal(length(coeffs), 4)
    linreg <- lm(income~age, data=PUMS5extract10000)

    output <- as.numeric(coeffs$coefficients[[1]][1]) #extracts coefficient from output
    expectedOutput <- as.numeric(linreg[[1]][2])
    expect_equal(floor(output), floor(expectedOutput)) #check floor of values due to fact that there is some noise added here
})
