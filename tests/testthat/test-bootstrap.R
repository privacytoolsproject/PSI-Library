library(PSIlence)
context("bootstrap")

test_that('bootstrap did not run, then thre NaNs, now produces result that is way too big', {
  data(PUMS5extract10000, package = "PSIlence")
  
  n.boot <- 25
  boot_mean <- dpMean$new(mechanism='mechanismBootstrap', var.type='numeric', 
                          variable='income', n=10000, epsilon=0.1, rng=c(0, 750000), 
                          n.boot=n.boot)
  boot_mean$release(PUMS5extract10000)
  print(boot_mean$result)
  
  print(mean(boot_mean$result$release))
})
