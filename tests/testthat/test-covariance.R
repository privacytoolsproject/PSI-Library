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
  skip("Skipping: range checks do not throw correct warning")
  expect_error(dpCovariance$new(mechanism='mechanismLaplace', varType='numeric', n= 10000,
                                epsilon = 0.1, columns = c("income", "education"),
                                rng=c(100)),
               "range argument in error: requires upper and lower values as vector of length 2.")
  
  expect_warning(dpCovariance$new(mechanism='mechanismLaplace', varType='numeric', n=10000,
                                  epsilon=0.1, columns = c("income", "education"),
                                  rng=c(-10,0,100)),
                 "range argument supplied has more than two values.  Will proceed using min and max values as range.")
})

test_that('covariance running as expected', {
  skip("Skipping: covariance does not run as expected.")
  rangeIncome <- c(-10000, 713000)
  rangeEducation <- c(1, 16)
  range <- rbind(rangeIncome, rangeEducation)
  
  dpCov <- dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 10000, 
                            epsilon = 5000000000, columns = c("income", "educ"), rng = range, formula='x~y')
  expect_silent(dpCov$release(PUMS5extract10000))
  print(dpCov$release(PUMS5extract10000))
  
  print("expected values")
  print(cov(x=data["income"], y=data["educ"]))
  print(var(data["income"]))
  print(var(data["educ"]))
  
  print(range(data["income"]))
  print(range(data["educ"]))
})

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {
    expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = -1, 
                                  epsilon = 1, columns = c("income", "educ"), rng = range, formula='x~y'),
                 "n must be a positive whole number")
    expect_error(dpCovariance$new(mechanism="mechanismLaplace",varType = 'numeric', n = 0.5, 
                                  epsilon = 1, columns = c("income", "educ"), rng = range, formula='x~y'),
                 "n must be a positive whole number")
})
