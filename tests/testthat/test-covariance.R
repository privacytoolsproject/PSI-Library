test_that('epsilon checks throw correct warning', {
  library(PSIlence)
  data(PUMS5extract10000, package = "PSIlence")
  range.income <- c(-10000, 713000)
  range.education <- c(1, 16)
  range <- rbind(range.income, range.education)
  expect_error(covariance.release(x = PUMS5extract10000, 
                                  var.type = 'numeric', n = 10000, 
                                  epsilon = -0.1, rng = range, 
                                  columns = c("income", "education"), 
                                  formulae = income ~ education), 
               "Privacy parameter epsilon must be a value greater than zero.")
  epsilon <- c(0.1, 0.5)
  expect_error(covariance.release(x = PUMS5extract10000, 
                                  var.type = 'numeric', n = 10000, 
                                  epsilon = epsilon, rng = range, 
                                  columns = c("income", "education"), 
                                  formulae = income ~ education), 
               paste("Privacy parameter epsilon must be a single value, but is currently a vector of length ", length(epsilon)))

})
