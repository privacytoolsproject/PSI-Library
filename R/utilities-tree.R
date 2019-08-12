#' Function to evaluate weights from the noise variance and standard errors in child nodes for the 
#'  node of a differentially private binary tree
#'
#' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' @param tree Data frame with binary tree attributes and node values
#' @param idx Index of the node for which the weight is evaluated
#' @return Weight

wBelow <- function(invSigmaSq, tree, idx) {
    leftIndex <- 2 * idx
    rightIndex <- leftIndex + 1
    w <- invSigmaSq / (invSigmaSq + 1 / (tree$seBelow[leftIndex]^2 + tree$seBelow[rightIndex]^2))
    return(w)
}


#' Function to evaluate weights from the noise variance and standard errors in a parent and adjacent 
#'  nodes for the node of a differentially private binary tree
#'
#' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' @param tree Data frame with binary tree attributes and node values
#' @param parent Index of the parnet node
#' @param adjacent Index of the adjacent node
#' @return Weight

wAbove <- function(invSigmaSq, tree, parent, adjacent) {
    w <- invSigmaSq / (invSigmaSq + 1 / (tree$seAbove[parent]^2 + tree$seBelow[adjacent]^2))
    return(w)
}


#' Function to evaluate weights efficiently using the noise variance and standard errors in parent and adjacent 
#'  nodes as well child nodes for the node of a differentially private binary tree
#'
#' @param tree Data frame with binary tree attributes and node values
#' @param idx Index of the node for which the weight is evaluated
#' @param parent Index of the parnet node
#' @param adjacent Index of the adjacent node
#' @return Weight

wEfficient <- function(tree, idx, parent, adjacent) {
    w <- tree$seBelow[idx]^(-2) / (tree$seBelow[idx]^(-2) + (1 / (tree$seAbove[parent]^2 + tree$seBelow[adjacent]^2)))
    return(w)
}


#' Function to estimate the nodes of a tree using noisy child nodes
#'
#' @param w Weight used construct the estimate
#' @param tree Data frame with binary tree attributes and node values
#' @param idx Index of the node for which the estimate is evaluated
#' @return Noisy node estimate

estBelow <- function(w, tree, idx) {
    leftIndex <- 2 * idx
    rightIndex <- leftIndex + 1
    est <- w * tree$noisy[idx] + (1 - w) * (tree$estBelow[leftIndex] + tree$estBelow[rightIndex])
    return(est)
}


#' Function to estimate the nodes of a tree using noisy parent and adjacent nodes
#'
#' @param w Weight used construct the estimate
#' @param tree Data frame with binary tree attributes and node values
#' @param idx Index of the node for which the estimate is evaluated
#' @param parent Index of the parnet node
#' @param adjacent Index of the adjacent node
#' @return Noisy node estimate

estAbove <- function(w, tree, idx, parent, adjacent) {
    est <- w * tree$noisy[idx] + (1 - w) * (tree$estAbove[parent] - tree$estBelow[adjacent])
    return(est)
}


#' Function to efficiently estimate the nodes of a tree using all available information in the tree
#'
#' @param w Weight used construct the estimate
#' @param tree Data frame with binary tree attributes and node values
#' @param idx Index of the node for which the estimate is evaluated
#' @param parent Index of the parnet node
#' @param adjacent Index of the adjacent node
#' @return Efficient noisy node estimate

estEfficient <- function(w, tree, idx, parent, adjacent) {
    est <- w * tree$estBelow[idx] + (1 - w) * (tree$estAbove[parent] - tree$estBelow[adjacent])
    return(est)
}


#' Function to evaluate the standard error of a node estimate given a weight and the standard
#'  deviation of the noise used to perturb the nodes
#'
#' @param w Weight used construct the estimate
#' @param sigma Standard deviation of the noise used to perturb the estimates
#' @return Standard error of the node estimate

stErr <- function(w, sigma) {
    return(sigma * sqrt(w))
}


#' Function to estimate a noisy binary tree from the terminal nodes
#'
#' @param tree Data frame with binary tree attributes and node values
#' @param terminalLevelIndex Index of the first terminal leaf node
#' @param nNodes Number of nodes in the binary tree
#' @param sigma Standard deviation of the noise used to perturb the estimates
#' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' @return Bottom-up estimate of noisy binary tree in vector form

estBottomUp <- function(tree, terminalLevelIndex, nNodes, sigma, invSigmaSq) {
    tree$estBelow <- c(rep(NA, (terminalLevelIndex - 1)), tree$noisy[terminalLevelIndex:nrow(tree)])
    tree$seBelow <- c(rep(NA, (terminalLevelIndex - 1)), rep(sigma, nNodes - (terminalLevelIndex - 1)))
    tree$wBelow <- rep(NA, nNodes)
    for (i in (terminalLevelIndex - 1):2) {
        tree$wBelow[i] <- wBelow(invSigmaSq, tree, i)
        tree$estBelow[i] <- estBelow(tree$wBelow[i], tree, i)
        tree$seBelow[i] <- stErr(tree$wBelow[i], sigma)
    }
    tree$estBelow[tree$estBelow < 0] <- 0
    return(tree)
}


#' Function to estimate a noisy binary tree from the top down
#'
#' @param tree Data frame with binary tree attributes and node values
#' @param n Number of observations in the vector represented by the binary tree
#' @param nNodes Number of nodes in the binary tree
#' @param sigma Standard deviation of the noise used to perturb the estimates
#' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' @return Top-down estimate of noisy binary tree in vector form

estTopDown <- function(tree, n, nNodes, sigma, invSigmaSq) {
    tree$estAbove <- c(n, rep(NA, (nNodes - 1)))
    tree$seAbove <- c(0, rep(NA, (nNodes - 1)))
    tree$wAbove <- rep(NA, nNodes)
    for (i in 2:nNodes) {
        tree$wAbove[i] <- wAbove(invSigmaSq, tree, tree$parent[i], tree$adjacent[i])
        tree$estAbove[i] <- estAbove(tree$wAbove[i], tree, i, tree$parent[i], tree$adjacent[i])
        tree$seAbove[i] <- stErr(tree$wAbove[i], sigma)
    }
    tree$estAbove[tree$estAbove < 0] <- 0
    return(tree)
}


#' Function to estimate a noisy binary tree efficiently using all available information in the tree
#'
#' @param tree Data frame with binary tree attributes and node values
#' @param n Number of observations in the vector represented by the binary tree
#' @param nNodes Number of nodes in the binary tree
#' @param sigma Standard deviation of the noise used to perturb the estimates
#' @param invSigmaSq Inverse variance of the noise used in perturbing nodes
#' @return Efficient estimate of noisy binary tree in vector form

estEfficiently <- function(tree, n, nNodes, sigma, invSigmaSq) {
    tree$estEfficient <- c(n, rep(NA, (nNodes - 1)))
    tree$seEfficient <- rep(NA, nNodes)
    tree$wEfficient <- rep(NA, nNodes)
    for (i in 2:nNodes) {
        tree$wEfficient[i] <- wEfficient(tree, i, tree$parent[i], tree$adjacent[i])
        tree$estEfficient[i] <- estEfficient(tree$wEfficient[i], tree, i, tree$parent[i], tree$adjacent[i])
        tree$seEfficient[i] <- stErr(tree$wEfficient[i], sigma)
    }
    tree$estEfficient[tree$estEfficient < 0] <- 0
    return(tree)
}


#' Function to evaluate a binary tree
#'
#' @param x Numeric vector to be represented as a binary tree in vector form
#' @param n Number of observations in \code{x}
#' @param rng An a priori estimate of the range of \code{x}
#' @param gran The granularity at which \code{x} is represented in the tree
#' @param universeSize Difference in the range of \code{x} over the granularity, plus 1
#' @param depth The depth of the binary tree
#' @return A binary tree in vector form

binaryTree <- function(x, n, rng, gran, universeSize, depth) {
    tree <- rep(0, times=(2^depth + universeSize))
    for (i in 1:n) {
        idx <- ((x[i] - rng[1]) / gran) + 2^depth
        tree[idx] <- tree[idx] + 1
    }
    d <- c()
    for (i in seq(2^depth, 2^depth - 1 + universeSize, 2)) {
        tree[i / 2] <- tree[i] + tree[i + 1]
        d <- c(d, depth)
    }
    depthCounter <- depth - 1
    while (depthCounter > 0) {
        for (i in seq(2^depthCounter, 2^(depthCounter + 1) - 1, 2)) {
            tree[i / 2] <- tree[i] + tree[i + 1]
            d <- c(d, depthCounter)
        }
        depthCounter <- depthCounter - 1
    } 
    tree <- data.frame(tree[1:(2^depth - 1)])
    names(tree) <- 'count'
    r <- c(0, rep(c(1, -1), nrow(tree) - 1))
    tree$depth <- 1
    tree$parent <- NA
    tree$adjacent <- NA
    for(i in 2:nrow(tree)) {
        tree$parent[i] <- trunc(i/2)
        tree$depth[i] <- trunc(log2(i)) + 1
        tree$adjacent[i] <- i + r[i]
    }
    return(tree)
}
