library(PSIlence)
context('Tree statistic')


test_that('Nodes initialize',{
  expect_equal(5,5)
  x <- Node$new(index=1,range=c(0,1))
  x$addChild()
})

test_that('Tree structure initializes properly', {
  expect_error(Tree$new(c(0,4), 4.3))
  expect_error(Tree$new(c(0,4), 17))
  expect_error(Tree$new(1,4))
  
  t <- Tree$new(c(0,4), 4)  
  leftChild = t$head$leftChild
  
  expect_equal(t$head$range,c(0,4))
  expect_equal(leftChild$range, c(0,2))
  expect_equal(leftChild$leftChild$range, c(0,1))
})

test_that('traverseLeft Tree method functions correctly',{
  t <- Tree$new(c(0,8),8)
  
  expect_true(t$traverseLeft(3,13))
  expect_true(t$traverseLeft(1,10))
  expect_true(t$traverseLeft(2,4))
  expect_true(t$traverseLeft(7,14))
  
  expect_false(t$traverseLeft(4,9))
  expect_false(t$traverseLeft(3,14))
  expect_false(t$traverseLeft(5,11))
  expect_false(t$traverseLeft(1,12))
}
)

test_that('Tree Statistic bins counts properly',{
  # general case
  x <- c(1,2,1,1)
  t <- publicTreeStatistic$new(x, c(0,4), 4)
  
  leftChild <- t$head$leftChild
  rightChild <- t$head$rightChild
  
  expect_equal(leftChild$weight,3)
  expect_equal(rightChild$weight,1)
  expect_equal(t$head$weight,length(x))
  
  #empty nodes
  x <- c(1)
  t <- publicTreeStatistic$new(x,c(0,4),4)
  
  expect_equal(t$head$rightChild$weight, 0)
  
})

test_that('traverseLeft for binning data points is correct',{
  t <- publicTreeStatistic$new(c(), c(0,4), 4)
  rightChild <- t$head$rightChild
  
  expect_true(t$traverseLeft(0.5, t$head))
  expect_false(t$traverseLeft(2.7, t$head))
  expect_false(t$traverseLeft(3.1, rightChild))
})

test_that('addValues functions correctly',{
  t <- Tree$new(c(0,4),2)
  toAdd <- c(1,3,7)
  t$addItems(toAdd)
  
  leftChild <- t$head$leftChild
  rightChild <- t$head$rightChild
  
  expect_equal(t$head$weight, toAdd[1])
  expect_equal(leftChild$weight, toAdd[2])
  expect_equal(rightChild$weight, toAdd[3])
  
  t$addItems(toAdd)
  
  expect_equal(t$head$weight, toAdd[1]*2)
  expect_equal(leftChild$weight, toAdd[2]*2)
  expect_equal(rightChild$weight, toAdd[3]*2)
})
