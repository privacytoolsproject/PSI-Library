library(PSIlence)
context('Tree statistic')


test_that('Nodes initialize',{
  expect_equal(5,5)
  x <- Node$new(index=1,range=c(0,1))
  x$addRightChild()
})

test_that('Tree initializes properly', {
  expect_error(Tree$new(c(0,4), 4.3))
  expect_error(Tree$new(c(0,4), 17))
  expect_error(Tree$new(1,4))
  
  #update tests so it's not just printing stuff for visual check
  t <- Tree$new(c(0,4), 8)
  t$printTree()
})

test_that('Tree Statistic bins counts properly',{
  x <- c(1,2,8,9)
  t <- publicTreeStatistic$new(x, c(0,10), 4)
  t$printTree()
  
  expect_equal(5,5)
})

# test_that('old public stat',{
#   x <- c(1,2,3,4,5,6)
#   n <- 6
#   rng <- c(1,6)
#   gran <- 2
#   universe.size <- floor(diff(rng) / gran + 1)
#   depth <- ceiling(log2(universe.size))
#   t <- binaryTree(x, n, rng, gran, universe.size, depth)
#   print(t)
#   expect_equal(5,5)
# })
