library(PSIlence)
context("laplace")

identityFun = function(x) {
  return(x)
}

x <- c(0,1,2,3,4)
sens <- 0

post = function(out){
  return(out)
}

test_that('Laplace mechanism outputs unperturbed values when sensitivity artificially set to zero', {
  mech = mechanismLaplace$new(epsilon=1, varType='numeric', rng=c(0,5), rngFormat='vector')
  out = mech$evaluate(identityFun, x, sens, post)
  expect_equal(out$release, x)
})

test_that('Laplace mechanism truncates to inputted range', {
  trimmedMech = mechanismLaplace$new(epsilon=1, varType='numeric', rng=c(0,3), rngFormat='vector')
  out = trimmedMech$evaluate(identityFun, x, sens, post)
  expect_equal(out$release, c(0,1,2,3,3))
})

test_that('Laplace mechanism correctly imputes data', {
  x <- c(0,1,2,3,NA)
  imputedMech <- mechanismLaplace$new(epsilon=1, varType='integer', rng=c(0,5), rngFormat='vector')
  out <- imputedMech$evaluate(identityFun, x, sens, post)
  expect_length(out$release, length(x))
})

# TODO: the way the fillMissing function was implemented within the Laplace mechanism seemed to indicate that you could have, e.g., characters, in your input array.
# But you can't add Laplace noise to chars/factors/bools...
