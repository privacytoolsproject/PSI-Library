#' Accuracy for a differentially private binary tree
#'
#' @param epsilon Numeric differential privacy parameter
#' @param rng Numeric a priori estimate of the variable range
#' @param gran Numeric granularity
#' @param alpha Numeric level of statistical significance, default 0.05
#' @return Accuracy guarantee for the tree given epsilon
#' @export tree.getAccuracy
#' @rdname tree.getAccuracy

tree.getAccuracy <- function(epsilon, rng, gran, alpha=0.05) {
    universe.size <- diff(rng) / gran + 1
    accuracy <- (2 * sqrt(2) / epsilon) * sqrt(log(2 / alpha)) * log2(universe.size)^(1.5)
    return(accuracy)
}


#' Epsilon for a differentially private binary tree
#'
#' @param accuracy Numeric accuracy needed
#' @param rng Numeric a priori estimate of the variable range
#' @param gran Numeric granularity
#' @param alpha Numeric level of statistical significance, default 0.05
#' @return Epsilon necessary to guarantee the given accuracy
#' @export tree.getParameters
#' @rdname tree.getParameters

tree.getParameters <- function(accuracy, rng, gran, alpha=0.05) {
    universe.size <- diff(rng) / gran + 1
    epsilon <- (2 * sqrt(2) / accuracy) * sqrt(log(2 / alpha)) * log2(universe.size)^(1.5)
    return(epsilon)
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


#' Function to efficiently estimate noisy node counts
#'
#' @param release The truncated differentially private noisy binary tree
#'      in vector form
#' @param tree.data Data frame with binary tree attributes, including depth
#'      and indicators of parent and adjacent nodes. Note that
#'      \code{nrow(tree.data) == length(release)}
#' @param n Number of observations
#' @param n.nodes Number of nodes in the binary tree, also \code{length(release)}
#' @param variance The variance of the noise used to perturb tree nodes
#' @param terminal.index Vector of indices corresponding to the terminal
#'      leaf nodes of the binary tree
#' @return Efficient differentially private binary tree

tree.postEfficient <- function(release, tree.data, n, variance, terminal.index) {
    n.nodes <- length(release)
    sigma <- sqrt(variance)
    inv.sigma.sq <- 1 / variance
    tree <- cbind(tree.data, release)
    names(tree)[ncol(tree)] <- 'noisy'
    tree <- estBottomUp(tree, min(terminal.index), n.nodes, sigma, inv.sigma.sq)
    tree <- estTopDown(tree, n, n.nodes, sigma, inv.sigma.sq)
    tree <- estEfficiently(tree, n, n.nodes, sigma, inv.sigma.sq)
    return(round(tree$est.efficient))
}


#' Differentially private binary tree
#'
#' @import methods
#' @export dpTree
#' @exportClass dpTree
#'
#' @include mechanisms.R

dpTree <- setRefClass(
    Class = 'dpTree',
    contains = 'mechanismLaplace'
)

dpTree$methods(
    initialize = function(mechanism, var.type, n, rng, gran, epsilon=NULL, accuracy=NULL, impute.rng=NULL, percentiles=NULL, alpha=0.05, ...) {
        .self$name <- 'Differentially private binary tree'
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$n <- n
        .self$rng <- rng
        .self$gran <- gran
        .self$alpha <- alpha
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- tree.getParameters(accuracy, rng, gran, alpha)
        } else {
            .self$epsilon <- epsilon
            .self$accuracy <- tree.getAccuracy(epsilon, rng, gran, alpha)
        }
        if (is.null(impute.rng)) {
            .self$impute.rng <- rng
        } else {
            .self$impute.rng <- impute.rng
        }
        .self$percentiles <- percentiles
})

dpTree$methods(
    release = function(x) {
        sens <- 2 * log2(diff(rng) / gran + 1)
        variance <- 2 * sens / epsilon
        universe.size <- floor(diff(rng) / gran + 1)
        depth <- ceiling(log2(universe.size))
        terminal.index <- seq(2^(depth - 1), 2^depth - 1)
        .self$result <- export(mechanism)$evaluate(.self$treeFun, x, sens, .self$postProcess, 
                                                   variance=variance, universe.size=universe.size, 
                                                   depth=depth, terminal.index=terminal.index, self=.self)
})

dpTree$methods(
    treeFun = function(x, universe.size, depth) {
        tree <- binaryTree(x, n, rng, gran, universe.size, depth)
        .self$tree.data <- tree[, which(names(tree) != 'count')]
        return(tree$count)
})

dpTree$methods(
    postProcess = function(out, ...) {
        out$release <- tree.postFormatRelease(out$release)
        ellipsis.vals <- getFuncArgs(list(...), tree.postEfficient)
        out$release <- do.call(tree.postEfficient, c(list(release=out$release, tree.data=tree.data, n=n), ellipsis.vals))
        ellipsis.vals <- getFuncArgs(list(...), tree.postCDF)
        out$cdf <- do.call(tree.postCDF, c(list(release=out$release, rng=rng), ellipsis.vals))
        out$median <- tree.postMedian(out$cdf)
        out$accuracy <- .self$accuracy
        out$epsilon <- .self$epsilon
        if (!is.null(percentiles)) {
            out$percentiles <- tree.postPercentiles(out$cdf, percentiles)
        }
        return(out)
})
