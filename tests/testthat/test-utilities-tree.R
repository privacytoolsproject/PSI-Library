context('tree utils')

test_that('adjacent element flip', {
  ls <- c(1,2,3,4)
  expect_equal(adjacentElements(ls), c(2,1,4,3))
})

test_that('estimation from below runs',{
  t <- list(c(10), c(6,4), c(3,3,1,3))
  
  w <- wBelow(t)
  expect_equal(w, c(4/7,2/3,1))
  
  c <- countBelow(t, w)
  expect_equal(t, c)
})

test_that('estimation from above runs', {
  t <- list(c(10), c(6,4), c(3,3,1,3))
  
  wB <- wBelow(t)
  wA <- wAbove(t, wB)
  expect_equal(wA, c(1, 5/8, 13/21))
  
  cB <- countBelow(t, wB)
  cA <- countAbove(t, cB, wA)
  expect_equal(t, cA)
})

test_that('optimal estimation runs', {
  t <- list(c(10), c(6,4), c(3,3,1,3))
  
  wB <- wBelow(t)
  wA <- wAbove(t, wB)
  cB <- countBelow(t, wB)
  cA <- countAbove(t, cB, wA)
  
  c <- optimalCount(t, wA, cA, cB)
  expect_equal(t,c)
})

test_that('optimal sigma est runs', {
  i <- inverseVariance(2,1)
  expect_equal(i, 1/2)
  
  t <- list(c(10), c(6,4), c(3,3,1,3))
  wB <- wBelow(t)
  wA <- wAbove(t, wB)
  s <- optimalSigma(wA, wB, 1)
  expect <- c(4/14, 1/3, 1/2) * sqrt(wA)
  expect_equal(expect, s)
})

test_that('optimal post-process script runs',{
  t <- list(c(10), c(6,4), c(3,3,1,3))
  out <- optimalPostProcess(t, 1)
  
  expect_equal(length(out), 4)
  expect_equal(out$optVariance, c(4/14, 1/3, 1/2) * sqrt(out$wAbove))
  expect_equal(out$optimalTree, t)
})

test_that('CDF created correctly', {
  # basic test case for public tree
  expect_equal(treePostCDF(list(c(10), c(6,4)),list(c(0,2,4)))$counts, c(0,6,10))
  expect_equal(treePostCDF(list(c(10), c(6,4)),list(c(0,2,4)))$proportions, c(0,6,10)/10)
  expect_equal(treePostCDF(list(c(10), c(6,4), c(3,3,1,3)),list(c(0,2,4,6,8)))$counts, c(0,3,6,7,10))
  
  # basic test cases for noisy tree to verify counts are pulled from correct spots
  expect_equal(treePostCDF(list(c(9), c(6,4)),list(c(0,2,4)))$counts, c(0,6,9))
  expect_equal(treePostCDF(list(c(10), c(6,3), c(1,4,2,5)), list(c(0,2,4,6,8)))$counts, c(0,1,6,8,10))
})

test_that('cdfMedian runs correctly', {
  t <- list(c(10), c(6,4), c(3,3,1,3))
  b <- list(c(0,2,4,6,8))
  cdf <- treePostCDF(t,b)
  med <- cdfMedian(cdf)
  expect_equal(med$val, 6)
  
  t <- list(c(10), c(5,5), c(2,3,3,2))
  b <- list(c(0,2,4,6,8))
  cdf <- treePostCDF(t,b)
  med <- cdfMedian(cdf)
  expect_equal(med$val, 4)
})

test_that('treeMean runs correctly', {
  t <- list(c(10), c(6,4), c(3,3,1,3))
  b <- list(c(0,2,4,6,8))
  expect_equal(treeMean(t, b), 3.8)
})