test_that('range check throws correct warning', {
  library(PSIlence)
  data(PUMS5extract10000, package = "PSIlence")
  expect_error(mean.release(PUMS5extract10000$income, 
                              var.type = 'numeric', 
                              epsilon = 0.1, n=10000,
                              rng = c(713000)), 
                 "range argument in error: requires upper and lower values as vector of length 2.")
  expect_warning(mean.release(PUMS5extract10000$income, 
                              var.type = 'numeric', 
                              epsilon = 0.1, n=10000,
                              rng = c(-10000, 13000, 713000)), 
                 "range argument supplied has more than two values.  Will proceed using min and max values as range.")
})

