#' Function to evaluate weights from the noise variance and standard errors in child nodes for the 
#'  node of a differentially private binary tree
#'
#' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' @param tree Data frame with binary tree attributes and node values
#' @param idx Index of the node for which the weight is evaluated
#' @return Weight

wBelow <- function(inv.sigma.sq, tree, idx) {
  left.idx <- 2 * idx
  right.idx <- left.idx + 1
  w <- inv.sigma.sq / (inv.sigma.sq + 1 / (tree$se.below[left.idx]^2 + tree$se.below[right.idx]^2))
  return(w)
}


#' Function to evaluate weights from the noise variance and standard errors in a parent and adjacent 
#'  nodes for the node of a differentially private binary tree
#'
#' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' @param tree Data frame with binary tree attributes and node values
#' @param parent Index of the parnet node
#' @param adjacent Index of the adjacent node
#' @return Weight

wAbove <- function(inv.sigma.sq, tree, parent, adjacent) {
  w <- inv.sigma.sq / (inv.sigma.sq + 1 / (tree$se.above[parent]^2 + tree$se.below[adjacent]^2))
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
  w <- tree$se.below[idx]^(-2) / (tree$se.below[idx]^(-2) + (1 / (tree$se.above[parent]^2 + tree$se.below[adjacent]^2)))
  return(w)
}


#' Function to estimate the nodes of a tree using noisy child nodes
#'
#' @param w Weight used construct the estimate
#' @param tree Data frame with binary tree attributes and node values
#' @param idx Index of the node for which the estimate is evaluated
#' @return Noisy node estimate

estBelow <- function(w, tree, idx) {
  left.idx <- 2 * idx
  right.idx <- left.idx + 1
  est <- w * tree$noisy[idx] + (1 - w) * (tree$est.below[left.idx] + tree$est.below[right.idx])
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
  est <- w * tree$noisy[idx] + (1 - w) * (tree$est.above[parent] - tree$est.below[adjacent])
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
  est <- w * tree$est.below[idx] + (1 - w) * (tree$est.above[parent] - tree$est.below[adjacent])
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
#' @param terminal.level.idx Index of the first terminal leaf node
#' @param n.nodes Number of nodes in the binary tree
#' @param sigma Standard deviation of the noise used to perturb the estimates
#' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' @return Bottom-up estimate of noisy binary tree in vector form

estBottomUp <- function(tree, terminal.level.idx, n.nodes, sigma, inv.sigma.sq) {
  tree$est.below <- c(rep(NA, (terminal.level.idx - 1)), tree$noisy[terminal.level.idx:nrow(tree)])
  tree$se.below <- c(rep(NA, (terminal.level.idx - 1)), rep(sigma, n.nodes - (terminal.level.idx - 1)))
  tree$w.below <- rep(NA, n.nodes)
  for (i in (terminal.level.idx - 1):2) {
    tree$w.below[i] <- wBelow(inv.sigma.sq, tree, i)
    tree$est.below[i] <- estBelow(tree$w.below[i], tree, i)
    tree$se.below[i] <- stErr(tree$w.below[i], sigma)
  }
  tree$est.below[tree$est.below < 0] <- 0
  return(tree)
}


#' Function to estimate a noisy binary tree from the top down
#'
#' @param tree Data frame with binary tree attributes and node values
#' @param n Number of observations in the vector represented by the binary tree
#' @param n.nodes Number of nodes in the binary tree
#' @param sigma Standard deviation of the noise used to perturb the estimates
#' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' @return Top-down estimate of noisy binary tree in vector form

estTopDown <- function(tree, n, n.nodes, sigma, inv.sigma.sq) {
  tree$est.above <- c(n, rep(NA, (n.nodes - 1)))
  tree$se.above <- c(0, rep(NA, (n.nodes - 1)))
  tree$w.above <- rep(NA, n.nodes)
  for (i in 2:n.nodes) {
    tree$w.above[i] <- wAbove(inv.sigma.sq, tree, tree$parent[i], tree$adjacent[i])
    tree$est.above[i] <- estAbove(tree$w.above[i], tree, i, tree$parent[i], tree$adjacent[i])
    tree$se.above[i] <- stErr(tree$w.above[i], sigma)
  }
  tree$est.above[tree$est.above < 0] <- 0
  return(tree)
}


#' Function to estimate a noisy binary tree efficiently using all available information in the tree
#'
#' @param tree Data frame with binary tree attributes and node values
#' @param n Number of observations in the vector represented by the binary tree
#' @param n.nodes Number of nodes in the binary tree
#' @param sigma Standard deviation of the noise used to perturb the estimates
#' @param inv.sigma.sq Inverse variance of the noise used in perturbing nodes
#' @return Efficient estimate of noisy binary tree in vector form

estEfficiently <- function(tree, n, n.nodes, sigma, inv.sigma.sq) {
  tree$est.efficient <- c(n, rep(NA, (n.nodes - 1)))
  tree$se.efficient <- rep(NA, n.nodes)
  tree$w.efficient <- rep(NA, n.nodes)
  for (i in 2:n.nodes) {
    tree$w.efficient[i] <- wEfficient(tree, i, tree$parent[i], tree$adjacent[i])
    tree$est.efficient[i] <- estEfficient(tree$w.efficient[i], tree, i, tree$parent[i], tree$adjacent[i])
    tree$se.efficient[i] <- stErr(tree$w.efficient[i], sigma)
  }
  tree$est.efficient[tree$est.efficient < 0] <- 0
  return(tree)
}


#' Function to evaluate a binary tree
#'
#' @param x Numeric vector to be represented as a binary tree in vector form
#' @param n Number of observations in \code{x}
#' @param rng An a priori estimate of the range of \code{x}
#' @param gran The granularity at which \code{x} is represented in the tree
#' @param universe.size Difference in the range of \code{x} over the granularity, plus 1
#' @param depth The depth of the binary tree
#' @return A binary tree in vector form

binaryTree <- function(x, n, rng, gran, universe.size, depth) {
  tree <- rep(0, times=(2^depth + universe.size))
  for (i in 1:n) {
    idx <- ((x[i] - rng[1]) / gran) + 2^depth
    tree[idx] <- tree[idx] + 1
  }
  d <- c()
  for (i in seq(2^depth, 2^depth - 1 + universe.size, 2)) {
    tree[i / 2] <- tree[i] + tree[i + 1]
    d <- c(d, depth)
  }
  depth.counter <- depth - 1
  while (depth.counter > 0) {
    for (i in seq(2^depth.counter, 2^(depth.counter + 1) - 1, 2)) {
      tree[i / 2] <- tree[i] + tree[i + 1]
      d <- c(d, depth.counter)
    }
    depth.counter <- depth.counter - 1
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