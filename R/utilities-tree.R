#' Calculates variance of Laplace noise. 
#' 
#' Recall that the variance of the Laplace distribution is 2b^2, where b is the classic
#' scaling parameter of the Laplace distribution (which is in fact the inverse of what we
#' in this library call the scaling parameter). Then, in our case, b = eps/sens.
#' 
#' @param sens Sensitivity of function that is to be perturbed
#' @param eps Privacy parameter used

inverseVariance <- function(sens, eps){
  b <- eps/sens
  return (2*b^2)
}

wBelow <- function(tree){
  weights <- numeric(length(tree)) # initialize with one weight per level of the tree
  i <- length(tree)
  while (i>0){
    if (i == length(tree)){
      weights[i] <- 1
    }
    else {
      prev <- weights[i+1]
      weights[i] <- (2*prev)/(2*weights+1)
    }
    i <- i - 1
  }
  return(weights)
}

countBelow <- function(tree, wBelows){
  counts <- vector('list', length(tree))
  i <- length(tree)
  while (i > 0){
    # taking advantage of R's vector operations to find all weighted counts per layer of tree at once
    if (i == length(tree)){
      counts[i] <- tree[i]
    }
    else{
      w <- wBelows[i]
      
      child <- tree[i+1]
      l <- length(child)
      childSum <- child(1:l-1) + child(2:l)              #sum all pairs of children counts (note this counts adjacent nodes that don't have same parent)
      childSum <- childSum[seq(1, l, by=2)]              # pick out the sums that are of children with same parents
      
      counts[i] <- w*tree[i]+(1-w)(childSum)
    }
    i <- i - 1
  }
  return(counts)
}

wAbove <- function(tree, wBelows){
  weights <- numeric(length(tree))
  i <- 1
  while (i <= length(tree)){
    if (i == 1){
      weights[i] <- 1
    }
    else {
      prevAbove <- weights[i-1]
      prevBelow <- wBelows[i]
      weights[i] <- 1/(1 + (prevAbove+prevBelow)^(-1))
    }
    i <- i + 1
  }
  return(weights)
}

countsAbove <- function(tree, countsBelow, wAboves){
  counts <- vector('list', length(tree))
  i <- 1
  while (i <= length(tree)){
    # taking advantage of R's vector operations to find all weighted counts per layer of tree at once
    if (i == 1){
      counts[i] <- tree[i]
    }
    else {
      w <- wAboves[i]
      
      parents <- rep(counts[i-1], each=2) #replicating each parent to make dimension correct
      adjacents <- adjacentElements(tree[i])
      counts[i] <- w*tree[i] + (1-w)(parents - adjacents)
    }
    i <- i + 1
  }
  return(i)
}

adjacentElements <- function(ls){
  adj <- vector("numeric", length=length(ls))
  i <- 1
  while (i <= length(ls)){
    if (i%%2 == 0){
      print('a')
      adj[i] = ls[i-1]
    }
    else {
      print('b')
      adj[i] = ls[i+1]
    }
    i <- i + 1
  }
  print(adj)
  return(adj)
}