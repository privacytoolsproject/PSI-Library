#' #' Function to efficiently estimate noisy node counts
#' #'
#' #' @param release The truncated differentially private noisy binary tree
#' #'      in vector form
#' #' @param treeData Data frame with binary tree attributes, including depth
#' #'      and indicators of parent and adjacent nodes. Note that
#' #'      \code{nrow(treeData) == length(release)}
#' #' @param n Number of observations
#' #' @param nNodes Number of nodes in the binary tree, also \code{length(release)}
#' #' @param variance The variance of the noise used to perturb tree nodes
#' #' @param terminalIndex Vector of indices corresponding to the terminal
#' #'      leaf nodes of the binary tree
#' #' @return Efficient differentially private binary tree
#' #' FIX THIS
#' treePostEfficient <- function(release, treeData, n, variance, terminalIndex) {
#'   nNodes <- length(release)
#'   sigma <- sqrt(variance)
#'   invSigmaSq <- 1 / variance
#'   tree <- cbind(treeData, release)
#'   names(tree)[ncol(tree)] <- 'noisy'
#'   tree <- estBottomUp(tree, min(terminalIndex), nNodes, sigma, invSigmaSq)
#'   tree <- estTopDown(tree, n, nNodes, sigma, invSigmaSq)
#'   tree <- estEfficiently(tree, n, nNodes, sigma, invSigmaSq)
#'   return(round(tree$est.efficient))
#' }
#' 
#' 
#' #' Function to truncate negative noisy node counts at zero
#' #'
#' #' @param release The differentially private noisy binary tree
#' #' @return Noisy binary tree truncated at zero
#' 
#' # treePostFormatRelease <- function(release) {
#' #     release <- round(release)
#' #     release[release < 0] <- 0
#' #     return(release)
#' # }
#' 
#' 
#' #' Function to derive CDF from efficient terminal node counts
#' #'
#' #' @param release Efficient differentially private binary tree
#' #' @param rng An a priori estimate of the range of the vector
#' #'      being represented as a binary tree
#' #' @param terminalIndex Vector of indices corresponding to the terminal
#' #'      leaf nodes of the binary tree
#' #' @return Differentially private estimate of the empirical cumulative
#' #'      distribution function
#' #' FIX THIS
#' treePostCDF <- function(release, rng, terminalIndex) {
#'   terminal <- release[terminalIndex]
#'   stepSize <- diff(rng) / length(terminal)
#'   cdfSteps <- seq(rng[1], rng[2], stepSize)
#'   cdf <- c(0, cumsum(terminal) / sum(terminal))
#'   cdf <- data.frame(list('val' = cdfSteps, 'cdf' = cdf))
#'   return(cdf)
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
#' #' FIX THIS
#' treePostMean <- function(cdf, rng) {
#'   ecdf <- cdf$cdf
#'   pdf <- sapply(2:length(ecdf), function(i) ecdf[i] - ecdf[i - 1])
#'   p <- c(ecdf[1], pdf) * cdf$val
#'   return(sum(p))
#' }
#' 
#' 
#' #' Function to evaluate the median using the DP CDF
#' #'
#' #' @param cdf Differentially private estimate of the empirical cumulative
#' #'      distribution function
#' #' @return Differentially private estimate of the median
#' #' FIX THIS
#' treePostMedian <- function(cdf) {
#'   outMedian <- treePostPercentiles(cdf, 0.5)$value
#'   return(outMedian)
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
#' #' FIX THIS
#' treePostPercentiles <- function(cdf, percentiles) {
#'   absArgMin <- function(q, cdf) {
#'     target <- abs(q - cdf$cdf)
#'     out <- cdf$val[which(target == min(target))]
#'     return(c(q, mean(out)))
#'   }
#'   outValues <- lapply(percentiles, absArgMin, cdf)
#'   outValues <- data.frame(do.call(rbind, outValues))
#'   names(outValues) <- c('percentile', 'value')
#'   return(outValues)
#' }
#' 
#' #' Function to evaluate weights from the noise variance and standard errors in child nodes for the 
#' #'  node of a differentially private binary tree
#' #'
#' #' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param idx Index of the node for which the weight is evaluated
#' #' @return Weight
#' 
#' wBelow <- function(invSigmaSq, tree, idx) {
#'     leftIndex <- 2 * idx
#'     rightIndex <- leftIndex + 1
#'     w <- invSigmaSq / (invSigmaSq + 1 / (tree$seBelow[leftIndex]^2 + tree$seBelow[rightIndex]^2))
#'     return(w)
#' }
#' 
#' 
#' #' Function to evaluate weights from the noise variance and standard errors in a parent and adjacent 
#' #'  nodes for the node of a differentially private binary tree
#' #'
#' #' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param parent Index of the parnet node
#' #' @param adjacent Index of the adjacent node
#' #' @return Weight
#' 
#' wAbove <- function(invSigmaSq, tree, parent, adjacent) {
#'     w <- invSigmaSq / (invSigmaSq + 1 / (tree$seAbove[parent]^2 + tree$seBelow[adjacent]^2))
#'     return(w)
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
#'     w <- tree$seBelow[idx]^(-2) / (tree$seBelow[idx]^(-2) + (1 / (tree$seAbove[parent]^2 + tree$seBelow[adjacent]^2)))
#'     return(w)
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
#'     leftIndex <- 2 * idx
#'     rightIndex <- leftIndex + 1
#'     est <- w * tree$noisy[idx] + (1 - w) * (tree$estBelow[leftIndex] + tree$estBelow[rightIndex])
#'     return(est)
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
#'     est <- w * tree$noisy[idx] + (1 - w) * (tree$estAbove[parent] - tree$estBelow[adjacent])
#'     return(est)
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
#'     est <- w * tree$estBelow[idx] + (1 - w) * (tree$estAbove[parent] - tree$estBelow[adjacent])
#'     return(est)
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
#'     return(sigma * sqrt(w))
#' }
#' 
#' 
#' #' Function to estimate a noisy binary tree from the terminal nodes
#' #'
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param terminalLevelIndex Index of the first terminal leaf node
#' #' @param nNodes Number of nodes in the binary tree
#' #' @param sigma Standard deviation of the noise used to perturb the estimates
#' #' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' #' @return Bottom-up estimate of noisy binary tree in vector form
#' 
#' estBottomUp <- function(tree, terminalLevelIndex, nNodes, sigma, invSigmaSq) {
#'     tree$estBelow <- c(rep(NA, (terminalLevelIndex - 1)), tree$noisy[terminalLevelIndex:nrow(tree)])
#'     tree$seBelow <- c(rep(NA, (terminalLevelIndex - 1)), rep(sigma, nNodes - (terminalLevelIndex - 1)))
#'     tree$wBelow <- rep(NA, nNodes)
#'     for (i in (terminalLevelIndex - 1):2) {
#'         tree$wBelow[i] <- wBelow(invSigmaSq, tree, i)
#'         tree$estBelow[i] <- estBelow(tree$wBelow[i], tree, i)
#'         tree$seBelow[i] <- stErr(tree$wBelow[i], sigma)
#'     }
#'     tree$estBelow[tree$estBelow < 0] <- 0
#'     return(tree)
#' }
#' 
#' 
#' #' Function to estimate a noisy binary tree from the top down
#' #'
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param n Number of observations in the vector represented by the binary tree
#' #' @param nNodes Number of nodes in the binary tree
#' #' @param sigma Standard deviation of the noise used to perturb the estimates
#' #' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' #' @return Top-down estimate of noisy binary tree in vector form
#' 
#' estTopDown <- function(tree, n, nNodes, sigma, invSigmaSq) {
#'     tree$estAbove <- c(n, rep(NA, (nNodes - 1)))
#'     tree$seAbove <- c(0, rep(NA, (nNodes - 1)))
#'     tree$wAbove <- rep(NA, nNodes)
#'     for (i in 2:nNodes) {
#'         tree$wAbove[i] <- wAbove(invSigmaSq, tree, tree$parent[i], tree$adjacent[i])
#'         tree$estAbove[i] <- estAbove(tree$wAbove[i], tree, i, tree$parent[i], tree$adjacent[i])
#'         tree$seAbove[i] <- stErr(tree$wAbove[i], sigma)
#'     }
#'     tree$estAbove[tree$estAbove < 0] <- 0
#'     return(tree)
#' }
#' 
#' 
#' #' Function to estimate a noisy binary tree efficiently using all available information in the tree
#' #'
#' #' @param tree Data frame with binary tree attributes and node values
#' #' @param n Number of observations in the vector represented by the binary tree
#' #' @param nNodes Number of nodes in the binary tree
#' #' @param sigma Standard deviation of the noise used to perturb the estimates
#' #' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' #' @return Efficient estimate of noisy binary tree in vector form
#' 
#' estEfficiently <- function(tree, n, nNodes, sigma, invSigmaSq) {
#'     tree$estEfficient <- c(n, rep(NA, (nNodes - 1)))
#'     tree$seEfficient <- rep(NA, nNodes)
#'     tree$wEfficient <- rep(NA, nNodes)
#'     for (i in 2:nNodes) {
#'         tree$wEfficient[i] <- wEfficient(tree, i, tree$parent[i], tree$adjacent[i])
#'         tree$estEfficient[i] <- estEfficient(tree$wEfficient[i], tree, i, tree$parent[i], tree$adjacent[i])
#'         tree$seEfficient[i] <- stErr(tree$wEfficient[i], sigma)
#'     }
#'     tree$estEfficient[tree$estEfficient < 0] <- 0
#'     return(tree)
#' }

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

parentNode <- function(tree, i, j){
  tree[[i-1]][floor(j/2)]
}

ChildNodeLeft <- function(tree, i, j){
  tree[[i+1]][j*2]
}

ChildNodeRight <- function(tree, i, j){
  tree[[i+1]][j*2+1]
}

wEfficient <- function(){
  
}
