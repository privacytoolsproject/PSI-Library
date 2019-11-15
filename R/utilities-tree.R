#### Optimal Post-Processed Counts ####

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

#' Creates array of adjacent elements in almost certainly a not optimally-efficient way
#'
#' @param ls 1D array of even length n
#'
#' @return 1D array of length n equal to ls with every other element swapped
#'
#' @examples
#' ls <- c(1,2,3,4)
#' adjacentElements(ls) #outputs c(2,1,3,4)
#' 
adjacentElements <- function(ls){
  adj <- vector("numeric", length=length(ls))
  i <- 1
  while (i <= length(ls)){
    if (i%%2 == 0){
      adj[i] = ls[i-1]
    }
    else {
      adj[i] = ls[i+1]
    }
    i <- i + 1
  }
  return(adj)
}

#' Recursive weight estimate from below
#'
#' See extra_docs/tree-post-processing for formula. This assumes that the variance is the same for every node at in the tree.
#' 
#' @param tree Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
#'
#' @return Single weight for each level of the tree.
#'
#' @examples
#' t <- list(c(10), c(6,4), c(3,3,1,3))
#' w <- wBelow(t) # should output c(4/7, 2/3, 1)
#' 
wBelow <- function(tree){
  weights <- numeric(length(tree)) # initialize with one weight per level of the tree
  i <- length(tree)
  while (i>0){
    if (i == length(tree)){
      weights[i] <- 1
    }
    else {
      prev <- weights[i+1]
      weights[i] <- (2*prev)/(2*prev+1)
    }
    i <- i - 1
  }
  return(weights)
}

#' Recursively compute counts from below
#' 
#' See extra_docs/tree-post-processing for formula. This assumes that the variance is the same for every node at in the tree.
#' 
#' @param tree Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
#' @param wBelows Array of weights of same length as tree, where the ith weight corresponds to the ith weight calculated from below.
#'
#' @return List of counts for each node of the tree.
#'
#' @examples
#' t <- list(c(10), c(6,4), c(3,3,1,3))
#' w <- wBelow(t)
#' c <- countBelow(t, w) # t and c should be equal, since counts in t have no added noise.
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
      child <- tree[[i+1]]
      l <- length(child)
      childSum <- child[1:l-1] + child[2:l]              #sum all pairs of children counts (note this counts adjacent nodes that don't have same parent)
      childSum <- childSum[seq(1, l, by=2)]              # pick out the sums that are of children with same parents
      counts[[i]] <- w*tree[[i]]+(1-w)*(childSum)
    }
    i <- i - 1
  }
  return(counts)
}

#' Recursive weight estimation from above
#'
#' See extra_docs/tree-post-processing for formula. This assumes that the variance is the same for every node at in the tree.
#' 
#' @param tree Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
#' @param wBelows Array of weights of same length as tree, where the ith weight corresponds to the ith weight calculated from below.
#'
#' @return Single weight for each level of the tree.
#'
#' @examples
#' t <- list(c(10), c(6,4), c(3,3,1,3))
#' wB <- wBelow(t)
#' wA <- wAbove(t, wB) # Should return c(1, 5/8, 13/21).
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

#' Recursively compute counts from above
#' 
#' See extra_docs/tree-post-processing for formula. This assumes that the variance is the same for every node at in the tree.
#' @param tree Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
#' @param countsBelow Array of counts of same length as tree, as calculated by countBelow function
#' @param wAboves Weights computed from above, assuming each node in tree has same variance.
#'
#' @return List of counts for each node of the tree.
#'
#' @examples
#' t <- list(c(10), c(6,4), c(3,3,1,3))
#' wB <- wBelow(t)
#' wA <- wAbove(t, wB)
#' cB <- countBelow(t, wB)
#' cA <- countAbove(t, cB, wA) #t and cA should be equal, since counts in t have no added noise.
countAbove <- function(tree, countsBelow, wAboves){
  counts <- vector('list', length(tree))
  i <- 1
  while (i <= length(tree)){
    # taking advantage of R's vector operations to find all weighted counts per layer of tree at once
    if (i == 1){
      counts[[i]] <- tree[[i]]
    }
    else {
      w <- wAboves[i]
      parents <- rep(counts[[i-1]], each=2) #replicating each parent to make dimension correct
      adjacents <- adjacentElements(tree[[i]])
      counts[[i]] <- w*tree[[i]] + (1-w)*(parents - adjacents)
    }
    i <- i + 1
  }
  return(counts)
}

#' Optimal counts for nodes of tree
#' 
#' See extra_docs/tree-post-processing for formula. This assumes that the variance is the same for every node at in the tree.
#'
#' @param tree Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
#' @param wA Array of weights calculated from above with wAbove
#' @param cA Array of counts calculated from below with countBelow 
#' @param wB Array of weights calculated from below with wBelow
#' @param cB Array of counts calculated from above with countAbove
#' 
#' @return List of counts for each node of the tree.
#'
#' @examples
#' t <- list(c(10), c(6,4), c(3,3,1,3))
#' wB <- wBelow(t)
#' wA <- wAbove(t, wB)
#' cB <- countBelow(t, wB)
#' cA <- countAbove(t, cB, wA)
#' c <- optimalCount(t, wA, cA, wB, cB) #will return t
#' 
optimalCount <- function(tree, wA, cA, wB, cB){
  wB <- wBelow(tree)
  cB <- countBelow(tree, wB)
  
  wA <- wAbove(tree, wB)
  cA <- countAbove(tree, cB, wA)
  
  counts <- vector('list', length(tree))
  i <- 1
  while (i <= length(tree)){
    if (i == 1){
      counts[[i]] <- tree[[i]]
    }
    else{
      w <- wA[i]
      parents <- rep(cA[[i-1]], each=2) #replicating each parent to make dimension correct
      adjacents <- adjacentElements(cB[[i]])
      counts[[i]] <- w*cB[[i]] + (1-w)*(parents - adjacents)
    }
    i <- i + 1
  }
  return(counts)
}

#' Optimal Sigma Estimation
#' 
#' See extra_docs/tree-post-processing for formula. This assumes that the variance is the same for every node at in the tree.
#'
#' @param wA Array of weights calculated from above
#' @param wB Array of weights calculated from below
#' @param epsilon Epsilon used for Laplace noise addition
#'
#' @return Value of variance for optimal estimate
optimalSigma <- function(wA, wB, epsilon){
  sB <- inverseVariance(2, epsilon)*wB #histogram of counts has sensitivity 2
  sOpt <- sB*sqrt(wA)
  return(sOpt)
}

optimalPostProcess <- function(tree, epsilon){
  wB <- wBelow(tree)
  wA <- wAbove(tree, wB)
  cB <- countBelow(tree, wB)
  cA <- countAbove(tree, cB, wA)
  
  out <- list('optimalTree'= optimalCount(tree, wA, cA, wB, cB))
  out$wBelow <- wB
  out$wAbove <- wA
  out$optVariance <- optimalSigma(wA, wB, epsilon)
  
  return(out)
}

### Post-Processed CDF ###

treePostCDF <- function(tree){
  print('a')
  # smallest granularity cdf possible uses leaf buckets
  vals <- tree[length(tree)][2]
  counts <- c()
  
  i <- 0
  print('b')
  while (i < length(vals)){
    print('c')
    # start out with min and max values at top of tree
    print(tree)
    m <- tree[1][2][1]
    M <- tree[2][2][2]
    print('d')
     # initialize count for i
    count <- 0
    # iterate through layers of tree
    index <- 1
    j <- 1
    print('e')
    while (j<=length(tree)){
      print('f')
      print(m)
      print(M)
      # determine if should traverse tree to left or right
      mid <- m + (M-m)/2
      print('g')
      # if looking at leftmost node in the tree, we know empirical cdf should evaluate to 0 not to bin size.
      if (i == 1){
        break
      }
      print('h')
      # if you don't need higher granularity, stop traversal
      if (vals[i] == M){
        count <- count + tree[j][1][index]
        break
      }
      print('i')
      # if at leaves of tree, record the count there
      if (j == len(tree)-1){
        count <- count + tree[j][1][index]
        }
      # print('j')
      # if traversing left
      else if (vals[i] <= mid){
        # reset max value to the mid, don't add to the count
        M <- mid
        # set next index of node to look at in next layer
        index <- index*2
        # if at end of tree, record count at that node ?
      }
      #print('k')
      #if traversing right, 
      else{
        # reset min value to the mid
        m <- mid
        # set to next index of node to look at in next layer
        index <- index*2 + 1
        count <- count + tree[j+1][0][index - 1] # add the node's left child to the count
      } 
      print('l')
      j <- j + 1
    }
    counts <- append(counts, count) # should change this to be more R friendly
    i <- i + 1
  }
  n <- tree[1][1][1] # pull public n from root of tree
  percents <- sapply(counts, (function(c) c/n))
  return (percents)
}

cdfMedian <- function(tree, cdf){
  if (0.5 %in% cdf){
    i <- which(2 == tree)
    val <- tree[length(tree)][2][i]
  }
  else{
    # otherwise, estimate cdf with closest value to 0.5th percentile
    distances <- sapply(cdf, (function(x) abs(x-0.5)))
    i <- which(min(distances) == distances)
    val <- tree[length(tree)][2][i]
  }
  return (val)
}