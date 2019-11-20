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

#' Optimal Post Processing
#' 
#' Wrapper function that generates optimal tree from noisy tree generated with the Laplace mechanism. The general idea is that
#' you can leverage the fact that the child counts at each node should sum to the parent count, and a node's count should be
#' equal to its parent count minus the count of the adjacent child node to generate less noisy counts at every node of the tree.
#' 
#' You can think of the leveraging of child node's information to get the parent node's counts as a recursive process that
#' begins at the leaves of a tree, which here we refer to with the tag "Below" in the helper function names. Similarly, leveraging
#' a parent node and the adjacent child nodes can be thought of as a recursive process that begins with the root node, which
#' is referred to here with the tage "Above" in the helper function names. A new count at every node in the tree can then be calculated
#' using the counts that are generated in this way, which each contribute an amount according to some weight, which are calculated
#' by the wBelow and wAbove functions respectively.
#' 
#' The theory behind this is explored in detail in the Honaker paper, whose implementation here is described in extra_docs/tree-post-processing. 
#' The implementation here assumes that the variance of the noise added is equal at every level of the tree, which also means that the weights
#' wBelow and wAbove are the same for nodes at the same level of the tree. Honaker provides a more general implementation in his work.
#' 
#' @references Honaker, James. "Efficient Use of Differentially Private Binary Trees." (2015).
#'
#' @param tree Differentially private tree generated from dpTree$release method
#' @param epsilon The epsilon value used for the noise addition at each node (note this is not the same as the global epsilon value.)
#'
#' @return A list that includes:
#'     optimalTree: a new tree with optimized counts at each level
#'     wBelow: the weights used to calculate the optimized counts from the leaf nodes to the top of the tree.
#'     wAbove: the weights used to calculate the optimized counts from the top of the tree to the leaf nodes.
#'     optVariance: the variance of the optimal tree counts.
#'
#' See extra_docs/tree-post-processing for more details on the meanings of wBelow, wAbove, and optVariance
optimalPostProcess <- function(tree, epsilon){
  wB <- wBelow(tree) # calculate the weights of counts from below
  wA <- wAbove(tree, wB) # calculate the weights of counts from above
  cB <- countBelow(tree, wB) # calculate counts from below
  cA <- countAbove(tree, cB, wA) # calculate counts from above
  
  out <- list('optimalTree'= optimalCount(tree, wA, cA, wB, cB)) # generate optimal count
  out$wBelow <- wB # save the weights used for counts from below
  out$wAbove <- wA # save the weights used for counts from above
  out$optVariance <- optimalSigma(wA, wB, epsilon) # calculate variance of optimal count
  
  return(out)
}

### Post-Processed CDF ###

#' Generates the least noisy CDF for the tree.
#' 
#' Note that for any numeric histogram, you can generate an empirical CDF for each of the maximal values of the bins by counting the number of items to the
#' left of that bin edge in the histogram. Here, we could do that by just using the histogram at the smallest level of the tree. This is not desirable
#' because you will have to sum many noisy things in order to get the counts. Instead, we can leverage the tree structure and minimize the number of counts
#' you need to sum for each count by traversing the tree.
#' 
#' For example, if a tree has bins with ranges [0,2), [2,4), [4,6), [6,8] at the leaf nodes, you could use the counts at the leaf nodes associated with [0,2) and [2,4)
#' to get an estimate for Pr(x < 4), or you could instead just use the count that is one level up on the tree which has range [0,4) to output this directly. Similarly,
#' if you want to get an estimate for Pr(x < 6), you could combine the count for [4,6) with the count for [0,4] to minimize the number of sums you must make. 
#' 
#' There are two edge cases here: the probability of a count that is less than the minimum of the tree is always 0 as we assume that the minimum and maximum values are always true, 
#' and the total sum is public knowledge.
#' 
#' Note that you could also minimize the number of counts you need to sum for each by calculating them as a difference from the total sum, and combine these two methods to get the most
#' optimal value for each of the counts. (E.g. in the previous example, you could get an estimate for Pr(x < 6) from the fact that you know Pr(x<8)=1 and the [6,8] count, which would
#' reduce the number of noisy counts used. We do not do this here. This code could certainly be optimized further in terms of runtime efficiency as well.
#'
#' @param tree (Optimized) differentially private tree of counts generated by dpTree$optimalPostProcess. Note that you could instead pass in the differentially private unoptimized tree
#' and the algorithm will still run correctly, but you will get more accurate results if you instead pass in the optimized version.
#' @param bins The bins of the dpTree object, which are stored as $binsByLevel.
#'
#' @return An empirical cdf that is as granular as the tree allows. The output consists of
#'    $bins: the values that the empirical cdf is evaluated at, which correspond to the ranges of the bins of the leaf nodes of the tree.
#'    $counts: the number of values that appear to the left of the corresponding value in $bins
#'    $proportions: the proportion of values that appear to the left of the corresponding value in $bins, i.e. proportions[i] is an approximation of Pr(x < bins[i])
treePostCDF <- function(tree, bins){
  counts <- c() # initializing counts
  vals <- bins[[length(bins)]] # vals are the values we can compute an empirical cdf for; i.e. the points that demarcate the bin ranges for the most granular bins in the tree
  i <- 1
  while (i <= length(vals)){ # generate the least noisy count for each of the values
    # start out with min and max values at top of tree
    m <- vals[1]
    M <- vals[length(vals)]

    #initialize count for i
    count <- 0
    
    # iterate through layers of tree
    index <- 1
    j <- 1
     while (j<=length(tree)){
        # determine if should traverse tree to left or right
        mid <- m + (M-m)/2
        if (i == 1){ # if looking at leftmost node in the tree, we know empirical cdf should evaluate to 0 not to bin size.
          count <- 0
          break
        }
        if (vals[i] == M){ # if you don't need higher granularity, stop traversal
          count <- count + as.numeric(tree[[j]][index]) # as.numeric is because R coerces the result of binary ops on named items to have the name of the first thing you are coercing to.
          break
        }
        if (j == length(tree)){ # if at leaves of tree, record the count there
          count <- count + as.numeric(tree[[j]][index])
          break
          }
        else if (vals[i] <= mid){ # if traversing left
          # reset max value to the mid, don't add to the count
          M <- mid
          # set next index of node to look at in next layer
          index <- index*2 - 1
        }
        else { #if traversing right
          # reset min value to the mid
          m <- mid
          # set to next index of node to look at in next layer
          index <- index*2
          count <- count + as.numeric(tree[[j+1]][index - 1]) # add the node's left child to the count
        }
        j <- j + 1
    }
    counts <- append(counts, count)
    i <- i + 1
  }
  out <- c()
  out$counts <- counts
  out$proportions <- counts/count[length(count)] # count[length(count)] is the total number of values in the data base (which is a publically known quantity)
  out$bins <- vals
  return (out)
}

#' (Estimated) Median of Input Data
#' 
#' Note that the median is just the cdf of a distribution evaluated at the 50th percentile. Given an empirical cdf genderated by the
#' treePostCDF function, we can report this value according to the noisy counts of the tree exactly if the 50th percentile is included.
#' If the granularity of the bins and counts of the data at the lowest level of the noisy tree are such that the 50th percentile is
#' not included in the output cdf, the function will instead output the closest count. 
#' 
#' E.g. if the input has
#' 
#' cdf$proportions <- c(0, 0.3. 0.4, 0.7, 1.0),
#' cdf$counts <- c(0,3,4,7,10)
#' 
#' then the function will output the cdf count associated with 0.4 as an estimate. The output is a tuple of both the estimated median, 
#' med$val, and the value of the proportion that was used to create the estimate. 
#'
#' @param cdf CDF generated with the treePostCDF function. 
#'
#' @return
#' @export
#'
#' @examples
cdfMedian <- function(cdf){
  med <- c()
  if (0.5 %in% cdf$proportions){
    med$val <- cdf$bins[cdf$proportions==0.5]
    med$prop <- 0.5
  }
  else{
    # otherwise, estimate cdf with closest value to 0.5th percentile
    distances <- sapply(cdf$proportions, (function(x) abs(x-0.5)))
    i <- which(min(distances) == distances)
    med$val <- cdf$bins[i+1]
    med$prop <- cdf$proportions[i]
  }
  return (med)
}

#' Estimated Mean of Input Data
#' 
#' Given a histogram of numeric values, you can estimate the mean of the underlying data by calculating the midpoint of each of the bins,
#' multiplying them by the bin's count, summing 
#' 
#' Let c_i be the count associated with bin i, and let mid_i be the midpoint of the data range that bin i represents.
#' Let n be the total number of data points. Then, 
#' \deqn{(1/n)\sum c_i * mid_i}  
#' is an estimate of the mean. 
#' 
#' Here, we calculate this estimate using the highest granularity bins possible, so the histogram created by the leaf nodes of the tree.
#' 
#' @param tree Differentially private tree of histogram counts generated with dpTree.
#' @param bins The bins by level of the tree, which corresponds to the binsByLevel attribute of dpTree.
#'
#' @return Estimate of the mean of the underlying data that created the dpTree.
treeMean <- function(tree, bins){
  #pull highest granularity bins from the list of all the tree's bins
  rng <- bins[[length(bins)]] 
  # find midpoints of each of the ranges of the bins
  mids <- (rng[2:length(rng)] - rng[1:length(rng)-1])/2 + rng[1:length(rng)-1]
  #estimate mean by weighting sum of bin midpoints by counts in lowest level of tree
  meanEst <- sum(mids * tree[[length(tree)]])/tree[[1]] 
  # convert to numeric to get rid of the bin labels that get coerced to the result due to R weirdness
  return(as.numeric(meanEst))
}