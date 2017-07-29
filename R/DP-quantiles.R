#' Binary tree
#' 
#' Binary tree in vector form.
#' 
#' @param x A numeric or integer vector.
#' @param var.type A character vector specifying variable type of \code{x}.
#'    Should be of length one.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param rng A numeric vector specifying an a priori estimate of the range
#'    of \code{x}. Should be of length two. 
#' @param gran 
#' @param universe.size 
#' @param sensitivity The difference of the range of \code{x} divided 
#'    by \code{n}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#'    
#' @return
#' @rdname binary.tree
#' @export
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
#' 
#' Function to calculate a differentially private release of quantiles.
#' 
#' @param x A numeric or integer vector.
#' @param var.type A character vector specifying variable type of \code{x}.
#'    Should be of length one.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param rng A numeric vector specifying an a priori estimate of the range
#'    of \code{x}. Should be of length two. 
#' @param gran
#' @param cdf.step
#' @param ... additional arguments passed to \code{quantile.release}.
#' 
#' @return 
#' @examples
#'  
#' @rdname quantile.release
#' @export
quantile.release <- function(x, var.type, n, epsilon, rng, gran, cdf.step, ...) {
    var.type <- check_variable_type(var.type, in_types=c('numeric', 'integer'))
    rng <- checkrange(rng)
    postlist <- list('accuracy' = 'getAccuracy',
                     'epsilon' = 'getParameters')
    sensitivity <- (2 * log2(diff(rng) / gran + 1))
    universe.size <- floor(diff(rng) / gran + 1)
    release <- mechanism.laplace(fun=binary.tree, x=x, var.type=var.type, rng=rng,
                                 sensitivity=sensitivity, epsilon=epsilon, n=n, gran=gran,
                                 cdf.step=cdf.step, universe.size=universe.size, postlist=postlist)
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


#' Quantiles accuracy
#' 
#' Get the accuracy of quantile statistic for a given value of epsilon.
#' 
#' The accuracy is interpreted as follows: The alpha value returned means that with 
#'   probability 1 - beta, simultaneously for all t with min <= t <= max, the algorithm's 
#'   estimate of the count in [min, t] is within alpha of the true value.
#'
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data.
#' @param universe.size Integer, the universe size
#' @param alpha A numeric vector specifying the statistical significance level.
#'    Default to 0.05.
#' 
#' @return The accuracy guaranteed by the given epsilon.
#' @author Nathan Manohar
#' @rdname quantile.getAccuracy
quantile.getAccuracy <- function(epsilon, n, universe.size, alpha=0.05) {
    accuracy <- ((4 / epsilon) * log(1 / alpha) * log2(universe.size)^1.5) / (n * 100)
    return(accuracy)
}


#' Quantiles epsilon
#' 
#' Get the epsilon value necessary to guarantee a desired level of accuracy of a quantile release.
#'
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data.
#' @param universe.size Integer, the universe size
#' @param alpha A numeric vector specifying the statistical significance level.
#'    Default to 0.05.
#'    
#' @return The epsilon value necessary to gaurantee the given accuracy.
#' @author Nathan Manohar
#' @rdname quantile.getParameters
quantile.getParameters <- function(accuracy, n, universe.size, alpha=0.05) {
    accuracy <- accuracy * 100  # added by JM on 8/5/14 for consistency among accuracies
    epsilon <- ((4 / accuracy) * log2(1 / alpha) * log2(universe.size)^1.5) / n
    return(epsilon)
}


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
