library(PSIlence)
context("Gaussian")

### Note: this is extremely cursory testing only to check that integrating new fillMissing function works properly 
#and should not be taken to be exhaustive testing of gaussian mechanism.

identityFun = function(x) {
  return(x)
}

x <- c(0,1,2,3,4)
sens <- 0

post = function(out){
  return(out)
}

mech <- mechanismGaussian$new(epsilon=1, delta=1, var.type='numeric', rng=c(0,5))
out <- mech$evaluate(identityFun, x, sens, post)
expect_equal(out$release, x)