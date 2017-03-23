#' Binary tree in vector form

binary.tree <- function(x, var.type, n, rng, gran, universe.size, sensitivity, epsilon) {
    depth <- ceiling(log2(universe.size))
    tree <- rep(0, times=(2^depth + universe.size))
    for (i in 1:n) {
        idx <- ((x[i] - rng[1]) / gran) + 2^depth
        tree[idx] <- tree[idx] + 1
    }
    for (i in seq(2^depth, 2^depth - 1 + universe.size, 2)) {
        tree[i / 2] <- tree[i] + tree[i + 1]
    }
    depth.counter <- depth - 1
    while (depth.counter > 0) {
        for (i in seq(2^depth.counter, 2^(depth.counter + 1) - 1, 2)) {
            tree[i / 2] <- tree[i] + tree[i + 1]
        }
        depth.counter <- depth.counter - 1
    }
    out <- list('name' = 'quantile',
                'stat' = tree,
                'n' = n,
                'rng' = rng,
                'gran' = gran,
                'universe.size' = universe.size,
                'sensitivity' = sensitivity,
                'epsilon' = epsilon)
}


#' Release differentially private quantiles

quantile.release <- function(x, var.type, n, epsilon, rng, gran, cdf.step) {
    var.type <- check_variable_type(var.type, in_types=c('numeric', 'integer'))
    rng <- checkrange(rng)
    sensitivity <- (2 * log2(diff(rng) / gran + 1))
    universe.size <- floor(diff(rng) / gran + 1)
    release <- mechanism.laplace(
        fun=binary.tree,
        x=x,
        var.type=var.type,
        rng=rng,
        sensitivity=sensitivity,
        epsilon=epsilon,
        n=n,
        gran=gran,
        cdf.step=cdf.step,
        universe.size=universe.size
    )
    # post-process the noisy tree
    terminal <- tail(release$release, universe.size)
    terminal <- sapply(terminal, function(x) {round(max(0, x))})
    cdf <- rep(0, times=length(seq(rng[1] + cdf.step, rng[2], cdf.step)))
    for (step.idx in 1:length(cdf)) {
        leaf.idx <- (((rng[1] + cdf.step * (step.idx + 1)) - rng[1]) / gran) + 1
        cdf[leaf.idx] <- sum(head(terminal, leaf.idx))
    }
    cdf <- cdf / max(cdf)  # normalize
    release$release <- cdf[1:max(which(cdf != 0))]  # remove trailing zeros in vector
    return(release)
}


#' Release differentially private quantiles
#'
#' @param x A vector of the data
#' @param epsilon Epsilon value for differential privacy
#' @param n Number of observations of the variable
#' @param range The range of the universe as a vector (min, max)
#' @param cdfstep The step sized used in outputting the approximate CDF; the values output are [min, min + cdfstep], [min, min + 2 * cdfstep], etc.
#' @param gran The granularity of the universe between the min and max. It is assumed that the data input is rounded to the granularity
#' @return A vector whose values are the approximate counts of [min, min + cdfstep], [min, min + 2 * cdfstep], etc.
#' @author Nathan Manohar
#'
#' Releases an approximate CDF of a data set in a differentially private manner. Returns these values in a vector. 
#'   This is done by outputting an approximate CDF of the data which can be used to approximate any quantile. 
#'
#' @examples 
#' n <- 1000
#' x <- round(runif(n, min=0, max=100))/100
#' quantile.release(x=x, epsilon=0.1, n=n, range=c(0,1), cdfstep=.1, gran=.01)

#quantile.release <- function(x, epsilon, n, range, cdfstep, gran){
#  
#  range <- checkrange(range)
#  data <- censordata(x=x, range=range)
#  
#  #Create binary tree stored as an array
#  universe_size <- ((range[2] - range[1]) / gran) + 1
#  tree_depth <- ceiling(log2(universe_size))
#  tree <- rep(0, times = (2^tree_depth + universe_size))
#  
#  #Add the counts of the data to the binary tree
#  for(i in 1:n) {
#    index <- ((as.numeric(data[i]) - range[1]) / gran) + 2^tree_depth
#    tree[index] <- tree[index] + 1
#  }
#  
#  #Sum adjacent nodes to complete binary tree
#  for(i in seq(2^tree_depth, 2^tree_depth - 1 + universe_size, 2)) {
#    tree[i/2] <- tree[i] + tree[i+1]    
#  }
#  
#  tree_counter <- tree_depth - 1
#  
#  while(tree_counter > 0) {
#    for(i in seq(2^tree_counter, 2^(tree_counter + 1) - 1, 2)) {
#        tree[i/2] <- tree[i] + tree[i+1]
#    }
#    tree_counter <- tree_counter - 1
#  }
#  
#
### Need to think about how to incorporate ranges
#
#  #Add Laplace noise to the nodes of the tree
#  treeSensitivity <- 2*log2(universe_size) # This is value in original source.  Check this.
#  for(i in 1:(2^tree_depth - 1 + universe_size)) {
#    tree[i] <- tree[i] + rlaplace(n=1, sensitivity=treeSensitivity, epsilon=epsilon)      # DOESN'T THIS EPSILON NEED TO BE PARTITIONED?
#  }
#
#  # Could loop be replaced with:
#  n.tree <- length(tree)  #ought to be: 2^tree_depth - 1 + universe_size
#  tree <- tree + rlaplace(n=n.tree, sensitivity=treeSensitivity, epsilon=epsilon)
#
#
#  # Or one way to use mechanism code
#  myfunction=function(x){x}
#  for(i in 1:(2^tree_depth - 1 + universe_size)) {
#    tree[i] <- mechanism.laplace(fun=myfunction, x=tree[i], sensitivity=treeSensitivity, epsilon=epsilon) 
#  }
#  
#
#
#  returnValue <- rep(0, times = length(seq(range[1] + cdfstep, range[2], cdfstep)))
#  
#  #Call computeInterval function to obtain the counts in the desired intervals
#  for(i in 1:length(returnValue)) {
#    returnValue[i] <- quantile.computeInterval(tree, range[1] + (cdfstep * i), range, gran)
#  }
#  
#  return(returnValue)
#}

#' Get the accuracy of quantile statistic for a given value of epsilon
#'
#' @param epsilon Epsilon value for DP
#' @param beta The true value is within the accuracy range with
#'    probability 1-beta
#' @param range The range of the universe as a vector (min, max)
#' @param gran The granularity of the universe between the min and max
#' @param n Number of observations of the variable
#' @return The accuracy guaranteed by the given epsilon
#' @author Nathan Manohar
#'
#' The accuracy is interpreted as follows: The alpha value returned means that with 
#'   probability 1 - beta, simultaneously for all t with min <= t <= max, the algorithm's 
#'   estimate of the count in [min, t] is within alpha of the true value.
  
#quantile.getAccuracy <- function(epsilon, beta, range, gran, n) {
  #range <- checkrange(range)
  #universe_size <- ((range[2] - range[1]) / gran) + 1
  #accuracy <- ((4/epsilon) * log2(1/beta) * log2(universe_size)^(1.5))/(n*100) #/(n*100) added by JM on 8/5/14 for consistency among accuracies 
  #return(accuracy) 
#}

#' Get the epsilon value necessary to guarantee a desired level of accuracy of a quantile release
#'
#' @param alpha the accuracy parameter
#' @param beta the true value is within the accuracy range (alpha)
#    with probability 1-beta
#' @param range The range of the universe as a vector (min, max)
#' @param gran The granularity of the universe between the min and max
#' @param n Number of observations of the variable
#' @return The epsilon value necessary to gaurantee the given accuracy
#' @author Nathan Manohar

#quantile.getParameters <- function(alpha, beta, range, gran, n) {
  #range <- checkrange(range)
  #alpha <- alpha*100   #added by JM on 8/5/14 for consistency among accuracies
  #universe_size <- ((range[2] - range[1]) / gran) + 1
  #epsilon <- ((4/alpha) * log2(1/beta) * log2(universe_size)^(1.5))/n
  #return(epsilon)
#}
  

#' Utility function to compute the count in a range from a binary tree
#'
#' @param tree The binary tree stored as an array
#' @param val The value we want to compute the count in [min, val] for
#' @param range The range of the universe as a vector (min, max)
#' @param gran The granularity of the universe between the min and max
#' @return The count in [min, val] given by the binary tree.
#' @author Nathan Manohar

#quantile.computeInterval <- function(tree, val, range, gran) {
#  
#  range <- checkrange(range)
#  universe_size <- ((range[2] - range[1]) / gran) + 1   
#  tree_depth <- ceiling(log2(universe_size))
#  index <- ((val - range[1])/gran) + 1
#  
#  interval_count <- 0
#  done <- 0
#  current_index <- 1
#  nodes_below <- 2^tree_depth
#  
#  while(done == 0) {
#    if(index == nodes_below) {
#        interval_count <- interval_count + tree[current_index]
#        done <- 1
#    }
#    else if(index <= nodes_below/2) {
#        nodes_below <- nodes_below/2
#        current_index <- 2 * current_index
#    }
#    else {
#        nodes_below <- nodes_below/2
#        interval_count <- interval_count + tree[2 * current_index]
#        current_index <- (2 * current_index) + 1
#        index <- index - nodes_below
#    }
#  }
#  
#  return(interval_count)   
#}
