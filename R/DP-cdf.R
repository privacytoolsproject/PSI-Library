#' Function to evaluate the CDF with a binary tree

dp.tree <- function(x, var.type, n, rng, epsilon, sensitivity, gran, variance) {

    universe.size <- floor(diff(rng) / gran + 1)
    depth <- ceiling(log2(universe.size))
    tree <- binaryTree(x, n, universe.size, depth)

    out <- (list('name' = 'tree',
                 'stat' = tree$count,
                 'tree.data' = tree[, which(names(tree) != 'count')],
                 'sensitivity' = sensitivity,
                 'gran' = gran,
                 'epsilon' = epsilon,
                 'rng' = rng,
                 'n.nodes' = nrow(tree),
                 'depth' = depth,
                 'terminal.index' = seq(2^(depth - 1), 2^depth - 1),
                 'sigma' = sqrt(variance),
                 'inv.sigma.sq' = 1 / variance))
    return(out)
}


#' Function to release a differentially private tree
#'
#' The perturbation takes place on the binary tree once constructed, then efficient estimation
#' takes place. This is why the efficient estimation is considered a post-processing step.

tree.release <- function(x, var.type, n, epsilon, rng, gran) {
    var.type <- check_variable_type(var.type, in_types=c('numeric', 'integer'))
    postlist <- list('efficient' = 'postEfficientTree')
    sensitivity <- 2 * log2(diff(rng) / gran + 1)
    variance = 2 * sensitivity / epsilon,
    tree <- mechanism.laplace(fun=dp.tree, x=x, var.type=var.type, rng=rng,
                              sensitivity=sensitivity, epsilon=epsilon, gran=gran,
                              cdf.step=cdf.step, variance=variance, postlist=postlist)
    return(release)
}

## to do
## [1] make sure no negative counts in noisy tree from mechanism output

tree.postEfficientTree <- function(release, tree.data, n.nodes, sigma, inv.sigma.sq, terminal.index) {
    tree <- cbind(tree.data, release)
    names(tree)[ncol(tree)] <- 'noisy' 
    tree <- estBottomUp(tree, min(terminal.index), n.nodes, sigma, inv.sigma.sq)
    tree <- estTopDown(tree, n.nodes, sigma, inv.sigma.sq)
    tree <- estEfficiently(tree, n.nodes, sigma, inv.sigma.sq)
    return(tree)
}
