library(PSIlence)
context("covariance")

data(PUMS5extract10000)

test_that('epsilon checks throw correct warning', {
  range.income <- c(-10000, 713000)
  range.education <- c(1, 16)
  range <- rbind(range.income, range.education)
  
  expect_error(dpCovariance$new(mechanism="mechanismLaplace",var.type = 'numeric', n = 10000, 
                                epsilon = -0.1, columns = c("income", "education"), rng = range),
               "Privacy parameter epsilon must be a value greater than zero.")
  
  expect_error(dpCovariance$new(mechanism="mechanismLaplace",var.type = 'numeric', n = 10000, 
                                epsilon = c(0.1, 0.5), columns = c("income", "education"), rng = range),
               "Privacy parameter epsilon must be a single value, but is currently a vector of length 2")
})

test_that('range checks throw correct warning', {
  expect_error(dpCovariance$new(mechanism='mechanismLaplace', var.type='numeric', n= 10000,
                                epsilon = 0.1, columns = c("income", "education"),
                                rng=c(100)),
               "range argument in error: requires upper and lower values as vector of length 2.")
  
  expect_warning(dpCovariance$new(mechanism='mechanismLaplace', var.type='numeric', n=10000,
                                  epsilon=0.1, columns = c("income", "education"),
                                  rng=c(-10,0,100)),
                 "range argument supplied has more than two values.  Will proceed using min and max values as range.")
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
  postLnReg <- covariance.postLinearRegression(formattedCovar, n, intercept, formula)
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
  postLnReg <- covariance.postLinearRegression(formattedCovar, n, intercept, formula)
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

  dpCov <- dpCovariance$new(mechanism="mechanismLaplace",var.type = 'numeric', n = 10000,
                            epsilon = 1, columns = c("income", "educ", 'age'), rng = range, formula='income~educ')
  out <- dpCov$release(PUMS5extract10000)
  expect_equal(length(out),3)
})

test_that('coefficient release function is correct', {
  range.income <- range(PUMS5extract10000['income'])
  range.education <- range(PUMS5extract10000['educ'])
  range.age <- range(PUMS5extract10000['age'])
  range <- rbind(range.income, range.education, range.age)
  
  dpCov <- dpCovariance$new(mechanism="mechanismLaplace",var.type = 'numeric', n = 10000,
                            epsilon = 10000000000, columns = c("income", "educ", "age"), rng = range, formula='income~educ')
  out <- dpCov$release(PUMS5extract10000)
  coeffs <- coefficient.release('income~age', out$release, n=10000)
  expect_equal(length(coeffs), 4)
  linreg <- lm(income~age, data=PUMS5extract10000)
  
  output <- as.numeric(coeffs$coefficients[[1]][1]) #extracts coefficient from output
  expectedOutput <- as.numeric(linreg[[1]][2])
  expect_equal(floor(output), floor(expectedOutput)) #check floor of values due to fact that there is some noise added here
})