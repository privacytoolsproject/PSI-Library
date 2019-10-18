context('tree utils')

test_that('adjacent element flip', {
  ls <- c(1,2,3,4)
  expect_equal(adjacentElements(ls), c(2,1,4,3))
})

test_that('estimation from below runs',{
  t <- list(c(10), c(6,4), c(3,3,1,3))
  
})