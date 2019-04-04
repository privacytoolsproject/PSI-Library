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

# test_that('range checks throw correct warning', {
#   
#   expect_error(dpCovariance$new(mechanism='mechanismLaplace', var.type='numeric', n= 10000, 
#                                 epsilon = 0.1, columns = c("income", "education"),
#                                 rng=c(100)), 
#                "range argument in error: requires upper and lower values as vector of length 2.")
#   
#   expect_warning(dpCovariance$new(mechanism='mechanismLaplace', var.type='numeric', n=10000, 
#                                   epsilon=0.1, columns = c("income", "education"),
#                                   rng=c(-10,0,100)),
#                  "range argument supplied has more than two values.  Will proceed using min and max values as range.")
# })

test_that('covariance running as expected', {
  range.income <- c(-10000, 713000)
  range.education <- c(1, 16)
  range <- rbind(range.income, range.education)
  
  dpCov <- dpCovariance$new(mechanism="mechanismLaplace",var.type = 'numeric', n = 10000, 
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
