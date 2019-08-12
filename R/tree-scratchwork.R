#' #' Function to evaluate weights from the noise variance and standard errors in child nodes for the 
#' #'  node of a differentially private binary tree
#' #'
#' #' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param idx Index of the node for which the weight is evaluated
#' #' @return Weight
#' 
#' wBelow <- function(inv.sigma.sq, tree, idx) {
#'   left.idx <- 2 * idx
#'   right.idx <- left.idx + 1
#'   w <- inv.sigma.sq / (inv.sigma.sq + 1 / (tree$se.below[left.idx]^2 + tree$se.below[right.idx]^2))
#'   return(w)
#' }
#' 
#' 
#' #' Function to evaluate weights from the noise variance and standard errors in a parent and adjacent 
#' #'  nodes for the node of a differentially private binary tree
#' #'
#' #' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param parent Index of the parnet node
#' #' @param adjacent Index of the adjacent node
#' #' @return Weight
#' 
#' wAbove <- function(inv.sigma.sq, tree, parent, adjacent) {
#'   w <- inv.sigma.sq / (inv.sigma.sq + 1 / (tree$se.above[parent]^2 + tree$se.below[adjacent]^2))
#'   return(w)
#' }
#' 
#' 
#' #' Function to evaluate weights efficiently using the noise variance and standard errors in parent and adjacent 
#' #'  nodes as well child nodes for the node of a differentially private binary tree
#' #'
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param idx Index of the node for which the weight is evaluated
#' #' @param parent Index of the parnet node
#' #' @param adjacent Index of the adjacent node
#' #' @return Weight
#' 
#' wEfficient <- function(tree, idx, parent, adjacent) {
#'   w <- tree$se.below[idx]^(-2) / (tree$se.below[idx]^(-2) + (1 / (tree$se.above[parent]^2 + tree$se.below[adjacent]^2)))
#'   return(w)
#' }
#' 
#' 
#' #' Function to estimate the nodes of a tree using noisy child nodes
#' #'
#' #' @param w Weight used construct the estimate
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param idx Index of the node for which the estimate is evaluated
#' #' @return Noisy node estimate
#' 
#' estBelow <- function(w, tree, idx) {
#'   left.idx <- 2 * idx
#'   right.idx <- left.idx + 1
#'   est <- w * tree$noisy[idx] + (1 - w) * (tree$est.below[left.idx] + tree$est.below[right.idx])
#'   return(est)
#' }
#' 
#' 
#' #' Function to estimate the nodes of a tree using noisy parent and adjacent nodes
#' #'
#' #' @param w Weight used construct the estimate
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param idx Index of the node for which the estimate is evaluated
#' #' @param parent Index of the parnet node
#' #' @param adjacent Index of the adjacent node
#' #' @return Noisy node estimate
#' 
#' estAbove <- function(w, tree, idx, parent, adjacent) {
#'   est <- w * tree$noisy[idx] + (1 - w) * (tree$est.above[parent] - tree$est.below[adjacent])
#'   return(est)
#' }
#' 
#' 
#' #' Function to efficiently estimate the nodes of a tree using all available information in the tree
#' #'
#' #' @param w Weight used construct the estimate
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param idx Index of the node for which the estimate is evaluated
#' #' @param parent Index of the parnet node
#' #' @param adjacent Index of the adjacent node
#' #' @return Efficient noisy node estimate
#' 
#' estEfficient <- function(w, tree, idx, parent, adjacent) {
#'   est <- w * tree$est.below[idx] + (1 - w) * (tree$est.above[parent] - tree$est.below[adjacent])
#'   return(est)
#' }
#' 
#' 
#' #' Function to evaluate the standard error of a node estimate given a weight and the standard
#' #'  deviation of the noise used to perturb the nodes
#' #'
#' #' @param w Weight used construct the estimate
#' #' @param sigma Standard deviation of the noise used to perturb the estimates
#' #' @return Standard error of the node estimate
#' 
#' stErr <- function(w, sigma) {
#'   return(sigma * sqrt(w))
#' }
#' 
#' 
#' #' Function to estimate a noisy binary tree from the terminal nodes
#' #'
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param terminal.level.idx Index of the first terminal leaf node
#' #' @param n.nodes Number of nodes in the binary tree
#' #' @param sigma Standard deviation of the noise used to perturb the estimates
#' #' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' #' @return Bottom-up estimate of noisy binary tree in vector form
#' 
#' estBottomUp <- function(tree, terminal.level.idx, n.nodes, sigma, inv.sigma.sq) {
#'   tree$est.below <- c(rep(NA, (terminal.level.idx - 1)), tree$noisy[terminal.level.idx:nrow(tree)])
#'   tree$se.below <- c(rep(NA, (terminal.level.idx - 1)), rep(sigma, n.nodes - (terminal.level.idx - 1)))
#'   tree$w.below <- rep(NA, n.nodes)
#'   for (i in (terminal.level.idx - 1):2) {
#'     tree$w.below[i] <- wBelow(inv.sigma.sq, tree, i)
#'     tree$est.below[i] <- estBelow(tree$w.below[i], tree, i)
#'     tree$se.below[i] <- stErr(tree$w.below[i], sigma)
#'   }
#'   tree$est.below[tree$est.below < 0] <- 0
#'   return(tree)
#' }
#' 
#' 
#' #' Function to estimate a noisy binary tree from the top down
#' #'
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param n Number of observations in the vector represented by the binary tree
#' #' @param n.nodes Number of nodes in the binary tree
#' #' @param sigma Standard deviation of the noise used to perturb the estimates
#' #' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' #' @return Top-down estimate of noisy binary tree in vector form
#' 
#' estTopDown <- function(tree, n, n.nodes, sigma, inv.sigma.sq) {
#'   tree$est.above <- c(n, rep(NA, (n.nodes - 1)))
#'   tree$se.above <- c(0, rep(NA, (n.nodes - 1)))
#'   tree$w.above <- rep(NA, n.nodes)
#'   for (i in 2:n.nodes) {
#'     tree$w.above[i] <- wAbove(inv.sigma.sq, tree, tree$parent[i], tree$adjacent[i])
#'     tree$est.above[i] <- estAbove(tree$w.above[i], tree, i, tree$parent[i], tree$adjacent[i])
#'     tree$se.above[i] <- stErr(tree$w.above[i], sigma)
#'   }
#'   tree$est.above[tree$est.above < 0] <- 0
#'   return(tree)
#' }
#' 
#' 
#' #' Function to estimate a noisy binary tree efficiently using all available information in the tree
#' #'
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param n Number of observations in the vector represented by the binary tree
#' #' @param n.nodes Number of nodes in the binary tree
#' #' @param sigma Standard deviation of the noise used to perturb the estimates
#' #' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' #' @return Efficient estimate of noisy binary tree in vector form
#' 
#' estEfficiently <- function(tree, n, n.nodes, sigma, inv.sigma.sq) {
#'   tree$est.efficient <- c(n, rep(NA, (n.nodes - 1)))
#'   tree$se.efficient <- rep(NA, n.nodes)
#'   tree$w.efficient <- rep(NA, n.nodes)
#'   for (i in 2:n.nodes) {
#'     tree$w.efficient[i] <- wEfficient(tree, i, tree$parent[i], tree$adjacent[i])
#'     tree$est.efficient[i] <- estEfficient(tree$w.efficient[i], tree, i, tree$parent[i], tree$adjacent[i])
#'     tree$se.efficient[i] <- stErr(tree$w.efficient[i], sigma)
#'   }
#'   tree$est.efficient[tree$est.efficient < 0] <- 0
#'   return(tree)
#' }
#' 
#' 
#' #' Function to evaluate a binary tree
#' #'
#' #' @param x Numeric vector to be represented as a binary tree in vector form
#' #' @param n Number of observations in \code{x}
#' #' @param rng An a priori estimate of the range of \code{x}
#' #' @param gran The granularity at which \code{x} is represented in the tree
#' #' @param universe.size Difference in the range of \code{x} over the granularity, plus 1
#' #' @param depth The depth of the binary tree
#' #' @return A binary tree in vector form
#' 
#' binaryTree <- function(x, n, rng, gran, universe.size, depth) {
#'   tree <- rep(0, times = (2 ^ depth + universe.size))
#'   for (i in 1:n) {
#'     idx <- ((x[i] - rng[1]) / gran) + 2^depth
#'     tree[idx] <- tree[idx] + 1
#'   }
#'   d <- c()
#'   for (i in seq(2^depth, 2^depth - 1 + universe.size, 2)) {
#'     tree[i / 2] <- tree[i] + tree[i + 1]
#'     d <- c(d, depth)
#'   }
#'   depth.counter <- depth - 1
#'   while (depth.counter > 0) {
#'     for (i in seq(2^depth.counter, 2^(depth.counter + 1) - 1, 2)) {
#'       tree[i / 2] <- tree[i] + tree[i + 1]
#'       d <- c(d, depth.counter)
#'     }
#'     depth.counter <- depth.counter - 1
#'   } 
#'   tree <- data.frame(tree[1:(2^depth - 1)])
#'   names(tree) <- 'count'
#'   r <- c(0, rep(c(1, -1), nrow(tree) - 1))
#'   tree$depth <- 1
#'   tree$parent <- NA
#'   tree$adjacent <- NA
#'   for(i in 2:nrow(tree)) {
#'     tree$parent[i] <- trunc(i/2)
#'     tree$depth[i] <- trunc(log2(i)) + 1
#'     tree$adjacent[i] <- i + r[i]
#'   }
#'   return(tree)
#' }
#' 

#' #' Accuracy for a differentially private binary tree
#' #'
#' #' @param epsilon Numeric differential privacy parameter
#' #' @param rng Numeric a priori estimate of the variable range
#' #' @param gran Numeric granularity
#' #' @param alpha Numeric level of statistical significance, default 0.05
#' #' @return Accuracy guarantee for the tree given epsilon
#' #' @export tree.getAccuracy
#' #' @rdname tree.getAccuracy
#' 
#' tree.getAccuracy <- function(epsilon, rng, gran, alpha=0.05) {
#'     universe.size <- diff(rng) / gran + 1
#'     accuracy <- (2 * sqrt(2) / epsilon) * sqrt(log(2 / alpha)) * log2(universe.size)^(1.5)
#'     return(accuracy)
#' }
#' 
#' 
#' #' Epsilon for a differentially private binary tree
#' #'
#' #' @param accuracy Numeric accuracy needed
#' #' @param rng Numeric a priori estimate of the variable range
#' #' @param gran Numeric granularity
#' #' @param alpha Numeric level of statistical significance, default 0.05
#' #' @return Epsilon necessary to guarantee the given accuracy
#' #' @export tree.getParameters
#' #' @rdname tree.getParameters
#' 
#' tree.getParameters <- function(accuracy, rng, gran, alpha=0.05) {
#'     universe.size <- diff(rng) / gran + 1
#'     epsilon <- (2 * sqrt(2) / accuracy) * sqrt(log(2 / alpha)) * log2(universe.size)^(1.5)
#'     return(epsilon)
#' }
#' 
#' 
#' #' Function to truncate negative noisy node counts at zero
#' #'
#' #' @param release The differentially private noisy binary tree
#' #' @return Noisy binary tree truncated at zero
#' 
#' tree.postFormatRelease <- function(release) {
#'     release <- round(release)
#'     release[release < 0] <- 0
#'     return(release)
#' }
#' 
#' 
#' #' Function to derive CDF from efficient terminal node counts
#' #'
#' #' @param release Efficient differentially private binary tree
#' #' @param rng An a priori estimate of the range of the vector
#' #'      being represented as a binary tree
#' #' @param terminal.index Vector of indices corresponding to the terminal
#' #'      leaf nodes of the binary tree
#' #' @return Differentially private estimate of the empirical cumulative
#' #'      distribution function
#' 
#' tree.postCDF <- function(release, rng, terminal.index) {
#'     terminal <- release[terminal.index]
#'     step.size <- diff(rng) / length(terminal)
#'     cdf.steps <- seq(rng[1], rng[2], step.size)
#'     cdf <- c(0, cumsum(terminal) / sum(terminal))
#'     cdf <- data.frame(list('val' = cdf.steps, 'cdf' = cdf))
#'     return(cdf)
#' }
#' 
#' 
#' #' Function to evaluate the mean using the DP CDF
#' #'
#' #' @param cdf Differentially private estimate of the empirical cumulative
#' #'      distribution function
#' #' @param rng Numeric a priori estimate of the range
#' #' @param gran Granularity
#' #' @return Differentially private estimate of the mean
#' 
#' tree.postMean <- function(cdf, rng) {
#'     ecdf <- cdf$cdf
#'     pdf <- sapply(2:length(ecdf), function(i) ecdf[i] - ecdf[i - 1])
#'     p <- c(ecdf[1], pdf) * cdf$val
#'     return(sum(p))
#' }
#' 
#' 
#' #' Function to evaluate the median using the DP CDF
#' #'
#' #' @param cdf Differentially private estimate of the empirical cumulative
#' #'      distribution function
#' #' @return Differentially private estimate of the median
#' 
#' tree.postMedian <- function(cdf) {
#'     out.median <- tree.postPercentiles(cdf, 0.5)$value
#'     return(out.median)
#' }
#' 
#' 
#' #' Quantile function using the DP CDF
#' #'
#' #' @param cdf Differentially private estimate of the empirical cumulative
#' #'      distribution function
#' #' @param percentiles Vector of probabilities given to the quantile function
#' #' @return Differnetially private estimate of the values corresponding to
#' #'      the provided probabilities
#' 
#' tree.postPercentiles <- function(cdf, percentiles) {
#'     absArgMin <- function(q, cdf) {
#'         target <- abs(q - cdf$cdf)
#'         out <- cdf$val[which(target == min(target))]
#'         return(c(q, mean(out)))
#'     }
#'     out.values <- lapply(percentiles, absArgMin, cdf)
#'     out.values <- data.frame(do.call(rbind, out.values))
#'     names(out.values) <- c('percentile', 'value')
#'     return(out.values)
#' }
#' 
#' 
#' #' Function to efficiently estimate noisy node counts
#' #'
#' #' @param release The truncated differentially private noisy binary tree
#' #'      in vector form
#' #' @param tree.data Data frame with binary tree attributes, including depth
#' #'      and indicators of parent and adjacent nodes. Note that
#' #'      \code{nrow(tree.data) == length(release)}
#' #' @param n Number of observations
#' #' @param n.nodes Number of nodes in the binary tree, also \code{length(release)}
#' #' @param variance The variance of the noise used to perturb tree nodes
#' #' @param terminal.index Vector of indices corresponding to the terminal
#' #'      leaf nodes of the binary tree
#' #' @return Efficient differentially private binary tree
#' 
#' tree.postEfficient <- function(release, tree.data, n, variance, terminal.index) {
#'     n.nodes <- length(release)
#'     sigma <- sqrt(variance)
#'     inv.sigma.sq <- 1 / variance
#'     tree <- cbind(tree.data, release)
#'     names(tree)[ncol(tree)] <- 'noisy'
#'     tree <- estBottomUp(tree, min(terminal.index), n.nodes, sigma, inv.sigma.sq)
#'     tree <- estTopDown(tree, n, n.nodes, sigma, inv.sigma.sq)
#'     tree <- estEfficiently(tree, n, n.nodes, sigma, inv.sigma.sq)
#'     return(round(tree$est.efficient))
#' }
#' 
#' 
#' #' Differentially private binary tree
#' #'
#' #' @param mechanism Character, the privacy mechanism.
#' #' @param var.type Character, the R variable type. One of \code{'numeric'} or
#' #'   \code{'integer'}.
#' #' @param Variable Character, variable name.
#' #' @param n Integer, number of observations.
#' #' @param rng Numeric, a priori estimate of the range.
#' #' @param gran Numeric, the granularity of the variable.
#' #' @param epsilon Numeric, privacy cost parameter.
#' #' @param accuracy Numeric, accuracy guarantee given \code{epsilon}.
#' #' @param impute.rng Numeric, range within which missing values are imputed. If \code{NULL},
#' #'   the range provided in \code{rng} is used.
#' #' @param percentiles Numeric, the percentiles to evaluate in post-processing. Optional, 
#' #'    default \code{NULL}.
#' #' @param alpha Numeric, the level of statistical significance. Default 0.05.
#' #'
#' #' @import methods
#' #' @export dpTree
#' #' @exportClass dpTree
#' #'
#' #' @include mechanism.R
#' #' @include mechanism-laplace.R
#' 
#' dpTree <- setRefClass(
#'     Class = 'dpTree',
#'     contains = 'mechanismLaplace'
#' )
#' 
#' dpTree$methods(
#'     initialize = function(mechanism, var.type, variable, n, rng, gran, epsilon=NULL,
#'                           accuracy=NULL, impute.rng=NULL, percentiles=NULL, alpha=0.05, ...) {
#'         .self$name <- 'Differentially private binary tree'
#'         .self$mechanism <- mechanism
#'         .self$var.type <- var.type
#'         .self$variable <- variable
#'         .self$n <- n
#'         .self$rng <- rng
#'         .self$gran <- gran
#'         .self$alpha <- alpha
#'         if (is.null(epsilon)) {
#'             .self$accuracy <- accuracy
#'             .self$epsilon <- tree.getParameters(accuracy, rng, gran, alpha)
#'         } else {
#'             .self$epsilon <- epsilon
#'             .self$accuracy <- tree.getAccuracy(epsilon, rng, gran, alpha)
#'         }
#'         if (is.null(impute.rng)) {
#'             .self$impute.rng <- rng
#'         } else {
#'             .self$impute.rng <- impute.rng
#'         }
#'         .self$percentiles <- percentiles
#' })
#' 
#' dpTree$methods(
#'     release = function(data) {
#'         x <- data[, variable]
#'         sens <- 2 * log2(diff(rng) / gran + 1)
#'         variance <- 2 * sens / epsilon
#'         universe.size <- floor(diff(rng) / gran + 1)
#'         depth <- ceiling(log2(universe.size))
#'         terminal.index <- seq(2^(depth - 1), 2^depth - 1)
#'         .self$result <- export(mechanism)$evaluate(.self$treeFun, x, sens, .self$postProcess, 
#'                                                    variance=variance, universe.size=universe.size, 
#'                                                    depth=depth, terminal.index=terminal.index, self=.self)
#' })
#' 
#' dpTree$methods(
#'     treeFun = function(x, universe.size, depth) {
#'         tree <- binaryTree(x, n, rng, gran, universe.size, depth)
#'         .self$tree.data <- tree[, which(names(tree) != 'count')]
#'         return(tree$count)
#' })
#' 
#' dpTree$methods(
#'     postProcess = function(out, ...) {
#'         out$variable <- variable
#'         out$release <- tree.postFormatRelease(out$release)
#'         ellipsis.vals <- getFuncArgs(list(...), tree.postEfficient)
#'         out$release <- do.call(tree.postEfficient, c(list(release=out$release, tree.data=tree.data, n=n), ellipsis.vals))
#'         ellipsis.vals <- getFuncArgs(list(...), tree.postCDF)
#'         out$cdf <- do.call(tree.postCDF, c(list(release=out$release, rng=rng), ellipsis.vals))
#'         out$mean <- tree.postMean(out$cdf, rng)
#'         out$median <- tree.postMedian(out$cdf)
#'         out$accuracy <- .self$accuracy
#'         out$epsilon <- .self$epsilon
#'         if (!is.null(percentiles)) {
#'             out$percentiles <- tree.postPercentiles(out$cdf, percentiles)
#'         }
#'         return(out)
#' })











