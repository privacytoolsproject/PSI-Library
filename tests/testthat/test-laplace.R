library(PSIlence)
context("laplace")

identityFun = function(x) {
  return(x)
}

x = c(0,1,2,3,4)
sens = 0

post = function(out){
  return(out)
}

mech = mechanismLaplace$new(epsilon=1, var.type='numeric', rng=c(0,5))
out = mech$evaluate(identityFun, x, sens, post)
expect_equal(out$release, x)

trimmedMech = mechanismLaplace$new(epsilon=1, var.type='numeric', rng=c(0,3))
out = trimmedMech$evaluate(identityFun, x, sens, post)
expect_equal(out$release, c(0,1,2,3,3))