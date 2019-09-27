context('tree statistic')

test_that('Public tree statistic runs',{
  x <- c(1:10)
  tree(x, 3)
})