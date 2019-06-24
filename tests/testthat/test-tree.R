library(PSIlence)
context('Tree statistic')

# 
# test_that('Nodes initialize',{
#   expect_equal(5,5)
#   x <- Node$new()
#   x$addRightChild()
# })
# 
test_that('Tree initializes properly', {
  expect_error(Tree$new(c(0,4), 4.3))
  expect_error(Tree$new(c(0,4), 17))
  expect_error(Tree$new(1,4))
  
  t <- Tree$new(c(0,4), 4)
  print(t)
})
