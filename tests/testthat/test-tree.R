context('tree statistic')

test_that('Tree binning runs',{
  out <- treeBins(c(0,10), 2, 10)
  expect_equal(length(out), 2)
  expect_equal(out[[1]], c(0,5,10))
})

test_that('Tree statistic initialization as expected', {
  x <- c(1:10)
  data.frame(x)
  stat <- dpTree$new('numeric', 'x', 10, 3, c(0,10), globalEps=1)
  expect_equal(stat$epsilon, 1/3)
})

test_that('Tree workflow runs', {
  x <- c(1:10)
  x <- data.frame(x)
  #will raise warning due to high epsilon
  stat <- expect_warning(dpTree$new('numeric', 'x', 10, 3, c(0,10), globalEps=10000))

  o <- stat$release(x)
  #print(o$optimalCounts)
  out <- o$release
  
  print(o)
  
  #expect_equal(length(out), 4)
  #expect_equal(out[[1]], 10)
  #expect_equal(length(out[[4]]), 8)
  #expect_equal(as.vector(round(out[[4]][8])), 2)
})