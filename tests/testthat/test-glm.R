test_that('output is a list', {
  library(PSIlence)
  
  data(PUMS5extract10000, package = "PSIlence")
  ols.test <- glm.release(PUMS5extract10000, nrow(PUMS5extract10000), epsilon = 0.1,
                          formula = income ~ educ + age, objective = dp.ols)
  expect_output(str(ols.test), "List of 7")
})