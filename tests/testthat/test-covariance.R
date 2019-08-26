library(PSIlence)
context("covariance")

data(PUMS5extract10000)

test_that('epsilon checks throw correct warning', {
  rangeIncome <- c(-10000, 713000)
  rangeEducation <- c(1, 16)
  range <- rbind(rangeIncome, rangeEducation)
  
  expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000, 
                                epsilon = -0.1, columns = c("income", "education"), rng = range),
               "Privacy parameter epsilon must be a value greater than zero.")
  
  expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000, 
                                epsilon = c(0.1, 0.5), columns = c("income", "education"), rng = range),
               "Privacy parameter epsilon must be a single value, but is currently a vector of length 2")
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

test_that('true covariance function is correct', {
  x1 <- c(0,5,3)
  x2 <- c(1,2,3)
  data <- data.frame(x1,x2)
  
  testCovar <- covar(data, intercept=FALSE)
  x<- covarianceFormatRelease(testCovar, c("x1", "x2"))
  y <- cov(data)
  
  expect_equal(as.matrix(x), y)
})

test_that('linear regression post-processing function is correct for 2x2 covariance matrix',{
  
  columns <- c('income', 'age')
  n <- 10000
  intercept <- FALSE
  formula <- 'income~age'
  data <- PUMS5extract10000[columns]
  
  covar <- covar(data, intercept=FALSE)

  formattedCovar <- covarianceFormatRelease(covar, columns)
  postLnReg <- covariancePostLinearRegression(formattedCovar, n, intercept, formula)
  output <- as.numeric(postLnReg[[1]][1]) #extracts coefficient from output
  trueLinearReg <- lm(formula, data=PUMS5extract10000)
  expectedOutput <- as.numeric(trueLinearReg[[1]][2]) #extracts coefficient from output
  expect_equal(floor(output), floor(expectedOutput)) #check floor of values due to floating point errors.
})

test_that('linear regression post-processing function is correct for 3x3 covariance matrix',{
  columns <- c('income', 'age', 'educ')
  n <- 10000
  intercept <- FALSE
  formula <- 'income~educ'
  data <- PUMS5extract10000[columns]
  
  covarMatrix <- covar(data, intercept=FALSE)
  
  formattedCovar <- covarianceFormatRelease(covarMatrix, columns)
  postLnReg <- covariancePostLinearRegression(formattedCovar, n, intercept, formula)
  output <- as.numeric(postLnReg[[1]][1]) #extracts coefficient from output
  
  trueLinearReg <- lm(formula, data=PUMS5extract10000)
  expectedOutput <- as.numeric(trueLinearReg[[1]][2]) #extracts coefficient from output
  expect_equal(floor(output), floor(expectedOutput)) #check floor of values due to floating point errors.
})

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

test_that('coefficient release function is correct', {
  range.income <- range(PUMS5extract10000['income'])
  range.education <- range(PUMS5extract10000['educ'])
  range.age <- range(PUMS5extract10000['age'])
  range <- rbind(range.income, range.education, range.age)
  
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

# make sure error thrown when n not positive or a whole number
test_that('make sure error thrown when n not positive or a whole number',{
  expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = -1, 
                                epsilon = 1, columns = c("income", "educ"), rng = range, formula='x~y'),
               "n must be a positive whole number")
  expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 0.5, 
                                epsilon = 1, columns = c("income", "educ"), rng = range, formula='x~y'),
               "n must be a positive whole number")
})
