library(PSIlence)
context("glm")

test_that('output is a list', {
  n <- 10000
  epsilon <- 0.5
  rng <- matrix(c(0, 200000, 18, 80, 0, 1), ncol=2, byrow=TRUE)
  form <- 'income ~ age + sex'
  
  model <- dpGLM$new(mechanism='mechanismObjective', varType='numeric', n=10000, epsilon=0.5, 
                     formula=form, rng=rng, objective='ols')
  out <- model$release(PUMS5extract10000)
  expect_output(str(out), "List of 4")
})

# make sure error thrown when n not positive or a whole number
test_that('error thrown when n not positive or whole number', {
    expect_error(dpGLM$new(mechanism='mechanismObjective', varType='numeric', n=-1, epsilon=0.5, 
                           formula=form, rng=rng, objective='ols'),
                 "n must be a positive whole number")
    expect_error(dpGLM$new(mechanism='mechanismObjective', varType='numeric', n=0.5, epsilon=0.5, 
                           formula=form, rng=rng, objective='ols'),
                 "n must be a positive whole number")
})
