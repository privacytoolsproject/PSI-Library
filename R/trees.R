wBelow <- function(inv.sigma.sq, tree, idx) {
    left.idx <- 2 * idx
    right.idx <- left.idx + 1
    w <- inv.sigma.sq / (inv.sigma.sq + 1 / (tree$se.below[left.idx]^2 + tree$se.below[right.idx]^2))
    return(w)
}

wAbove <- function(inv.sigma.sq, tree, parent, adjacent) {
    w <- inv.sigma.sq / (inv.sigma.sq + 1 / (tree$se.above[parent]^2 + tree$se.below[adjacent]^2))
    return(w)
}

wEfficient <- function(tree, idx, parent, adjacent) {
    w <- tree$se.below[idx]^(-2) / (tree$se.below[idx]^(-2) + (1 / (tree$se.above[parent]^2 + tree$se.below[adjacent]^2)))
    return(w)
}

estBelow <- function(w, tree, idx) {
    left.idx <- 2 * idx
    right.idx <- left.idx + 1
    est <- w * tree$noisy[idx] + (1 - w) * (tree$est.below[left.idx] + tree$est.below[right.idx])
    return(est)
}

estAbove <- function(w, tree, idx, parent, adjacent) {
    est <- w * tree$noisy[idx] + (1 - w) * (tree$est.above[parent] - tree$est.below[adjacent])
    return(est)
}

estEfficient <- function(w, tree, idx, parent, adjacent) {
    est <- w * tree$est.below[idx] + (1 - w) * (tree$est.above[parent] - tree$est.below[adjacent])
    return(est)
}

stErr <- function(w, sigma) {
    return(sigma * sqrt(w))
}

estBottomUp <- function(tree, terminal.level.idx, n.nodes, st.dev, inv.sigma.sq) {
    tree$est.below <- c(rep(NA, (terminal.level.idx - 1)), tree$noisy[terminal.level.idx:nrow(tree)])
    tree$se.below <- c(rep(NA, (terminal.level.idx - 1)), rep(st.dev, n.nodes - (terminal.level.idx - 1)))
    tree$w.below <- rep(NA, n.nodes)
    for (i in (terminal.level.idx - 1):2) {
        tree$w.below[i] <- wBelow(inv.sigma.sq, tree, i)
        tree$est.below[i] <- estBelow(tree$w.below[i], tree, i)
        tree$se.below[i] <- stErr(tree$w.below[i], st.dev)
    }
    tree$est.below[tree$est.below < 0] <- 0
    return(tree)
}

estTopDown <- function(tree, n, n.nodes, st.dev, inv.sigma.sq) {
    tree$est.above <- c(n, rep(NA, (n.nodes - 1)))
    tree$se.above <- c(0, rep(NA, (n.nodes - 1)))
    tree$w.above <- rep(NA, n.nodes)
    for (i in 2:n.nodes) {
        tree$w.above[i] <- wAbove(inv.sigma.sq, tree, tree$parent[i], tree$adjacent[i])
        tree$est.above[i] <- estAbove(tree$w.above[i], tree, i, tree$parent[i], tree$adjacent[i])
        tree$se.above[i] <- stErr(tree$w.above[i], st.dev)
    }
    tree$est.above[tree$est.above < 0] <- 0
    return(tree)
}

estEfficiently <- function(tree, n, n.nodes, st.dev, inv.sigma.sq) {
    tree$est.efficient <- c(n, rep(NA, (n.nodes - 1)))
    tree$se.efficient <- rep(NA, n.nodes)
    tree$w.efficient <- rep(NA, n.nodes)
    for (i in 2:n.nodes) {
        tree$w.efficient[i] <- wEfficient(tree, i, tree$parent[i], tree$adjacent[i])
        tree$est.efficient[i] <- estEfficient(tree$w.efficient[i], tree, i, tree$parent[i], tree$adjacent[i])
        tree$se.efficient[i] <- stErr(tree$w.efficient[i], st.dev)
    }
    tree$est.efficient[tree$est.efficient < 0] <- 0
    return(tree)
}

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
