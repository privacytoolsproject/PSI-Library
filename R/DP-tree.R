#' Evaluate a differentially private binary tree
#'
#' @param x Vector of numeric observations
#' @param var.type Character vector specifying the variable type, should be
#'      one of \code{c('numeric', 'integer')}
#' @param n Number of observations
#' @param rng An a priori estimate of the range of \code{x}
#' @param epsilon Privacy parameter epsilon, should be between zero and one
#' @param sensitivity The sensitivity of the statistic
#' @param gran The granularity at which \code{x} is represented in the tree
#' @param variance The variance of the noise used to perturb tree nodes
#' @param percentiles Vector of percentiles used in the post-processing
#'      quantile function
#' @return List of values including the true value of the binary tree and the associated data
#'      needed for the noisy release and post-processing

dp.tree <- function(x, var.type, n, rng, epsilon, sensitivity, gran, variance, percentiles) {

    universe.size <- floor(diff(rng) / gran + 1)
    depth <- ceiling(log2(universe.size))
    tree <- binaryTree(x, n, rng, gran, universe.size, depth)

    out <- (list('name' = 'tree',
                 'stat' = tree$count,
                 'tree.data' = tree[, which(names(tree) != 'count')],
                 'n' = n,
                 'sensitivity' = sensitivity,
                 'gran' = gran,
                 'epsilon' = epsilon,
                 'rng' = rng,
                 'n.nodes' = nrow(tree),
                 'depth' = depth,
                 'terminal.index' = seq(2^(depth - 1), 2^depth - 1),
                 'sigma' = sqrt(variance),
                 'inv.sigma.sq' = 1 / variance,
                 'percentiles' = percentiles))
    return(out)
}


#' Function to release a differentially private tree and post-processing
#'
#' @param x Vector of numeric observations
#' @param var.type Character vector specifying the variable type, should be
#'      one of \code{c('numeric', 'integer')}
#' @param n Number of observations
#' @param epsilon Privacy parameter epsilon, should be between zero and one
#' @param rng An a priori estimate of the range of \code{x}
#' @param impute.rng Numeric range within which to impute missing values in \code{x}
#' @param gran The granularity at which \code{x} is represented in the tree
#' @param percentiles Vector of percentiles used in the post-processing
#'      quantile function. The default is \code{NULL}, in which case the
#'      quantile function is not executed.
#' @return List of values including the differentially private tree, the
#'      cumulative distribution function, the median, and optionally a
#'      vector of percentiles. Other attributes of the binary tree are
#'      also included.
#' @export
tree.release <- function(x, var.type, n, epsilon, rng, impute.rng, gran, percentiles=NULL) {
    var.type <- check_variable_type(var.type, in_types=c('numeric', 'integer'))
    rng <- checkrange(rng)
    impute.rng <- ifelse(is.null(impute.rng), rng, impute.rng)
    postlist <- list('release' = 'postFormatRelease',
                     'release' = 'postEfficientTree',
                     'cdf' = 'postCDF',
                     'median' = 'postMedian')
    if (!is.null(percentiles)) {
        postlist <- c(postlist, list('percentiles' = 'postPercentiles'))
    }
    sensitivity <- 2 * log2(diff(rng) / gran + 1)
    variance <- 2 * sensitivity / epsilon
    release <- mechanism.laplace(fun=dp.tree, x=x, var.type=var.type, n=n, rng=rng,
                                 sensitivity=sensitivity, epsilon=epsilon, gran=gran,
                                 variance=variance, percentiles=percentiles, postlist=postlist)
    return(release)
}


#' Function to truncate negative noisy node counts at zero
#'
#' @param release The differentially private noisy binary tree
#' @return Noisy binary tree truncated at zero

tree.postFormatRelease <- function(release) {
    release <- round(release)
    release[release < 0] <- 0
    return(release)
}


#' Function to efficiently estimate noisy node counts
#'
#' @param release The truncated differentially private noisy binary tree
#'      in vector form
#' @param tree.data Data frame with binary tree attributes, including depth
#'      and indicators of parent and adjacent nodes. Note that
#'      \code{nrow(tree.data) == length(release)}
#' @param n Number of observations
#' @param n.nodes Number of nodes in the binary tree, also \code{length(release)}
#' @param sigma The standard deviation of the noise used in perturbing nodes
#' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' @param terminal.index Vector of indices corresponding to the terminal
#'      leaf nodes of the binary tree
#' @return Efficient differentially private binary tree

tree.postEfficientTree <- function(release, tree.data, n, n.nodes, sigma, inv.sigma.sq, terminal.index) {
    tree <- cbind(tree.data, release)
    names(tree)[ncol(tree)] <- 'noisy' 
    tree <- estBottomUp(tree, min(terminal.index), n.nodes, sigma, inv.sigma.sq)
    tree <- estTopDown(tree, n, n.nodes, sigma, inv.sigma.sq)
    tree <- estEfficiently(tree, n, n.nodes, sigma, inv.sigma.sq)
    return(round(tree$est.efficient))
}


#' Function to derive CDF from efficient terminal node counts
#'
#' @param release Efficient differentially private binary tree
#' @param rng An a priori estimate of the range of the vector
#'      being represented as a binary tree
#' @param terminal.index Vector of indices corresponding to the terminal
#'      leaf nodes of the binary tree
#' @return Differentially private estimate of the empirical cumulative
#'      distribution function

tree.postCDF <- function(release, rng, terminal.index) {
    terminal <- release[terminal.index]
    step.size <- diff(rng) / length(terminal)
    cdf.steps <- seq(rng[1], rng[2], step.size)
    cdf <- c(0, cumsum(terminal) / sum(terminal))
    cdf <- data.frame(list('val' = cdf.steps, 'cdf' = cdf))
    return(cdf)
}


#' Function to evaluate the median using the DP CDF
#'
#' @param cdf Differentially private estimate of the empirical cumulative
#'      distribution function
#' @return Differentially private estimate of the median

tree.postMedian <- function(cdf) {
    out.median <- tree.postPercentiles(cdf, 0.5)$value
    return(out.median)
}


#' Quantile function using the DP CDF
#'
#' @param cdf Differentially private estimate of the empirical cumulative
#'      distribution function
#' @param percentiles Vector of probabilities given to the quantile function
#' @return Differnetially private estimate of the values corresponding to
#'      the provided probabilities

tree.postPercentiles <- function(cdf, percentiles) {
    absArgMin <- function(q, cdf) {
        target <- abs(q - cdf$cdf)
        out <- cdf$val[which(target == min(target))]
        return(c(q, mean(out)))
    }
    out.values <- lapply(percentiles, absArgMin, cdf)
    out.values <- data.frame(do.call(rbind, out.values))
    names(out.values) <- c('percentile', 'value')
    return(out.values)
}
