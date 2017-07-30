#' Draw cryptographically secure random variates from uniform distribution
#'
#' @param n Integer giving number of variates needed
#' @param seed Integer indicating a seed for R's PNRG, defaults to \code{NULL}
#'
#' Draws secure random variates from the uniform distribution through \code{openssl}.
#' If a seed is provided, the \code{runif} function is used to draw the random variates.
#'
#' @examples
#' uniform_secure <- dpUnif(n=1000)
#' uniform_repeatable <- dpUnif(n=1, seed=75436)

dpUnif <- function(n, seed=NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
        return(runif(n))
    }
    return(openssl::rand_num(n))
}


#' Draw cryptographically secure random variates
#'
#' @param n Integer giving number of variates needed
#' @param scale Numeric scale for the distribution
#' @param dist Character specifying the distribution from which to draw the noise
#' @param seed Integer indicating a seed for R's PNRG, defaults to \code{NULL}
#'
#' @examples
#' laplace_noise <- dpNoise(n=1000, scale=1, dist='laplace')
#' gaussian_noise <- dpNoise(n=1000, scale=1, dist='gaussian')
#' laplace_noise_repeatable <- dpNoise(n=1, scale=1, dist='Laplace', seed=96845)

dpNoise <- function(n, scale, dist, seed=NULL) {
    u <- dpUnif(n, seed)
    if (dist == 'laplace') {
        return(qlap(u, b=scale))
    } else if (dist == 'gaussian') {
        return(qnorm(u, sd=scale))
    } else {
        stop(sprintf('Distribution "%s" not understood', dist))
    }
}



#' Random draw from Laplace distribution
#'
#' @param sensitivity numeric
#' @param epsilon numeric
#' @param n integer, number of draws
#' @return Random draws from Laplace distribution
#' @examples
#' rlaplace(sensitivity=1, epsilon=0.1)

rlaplace = function(n=1, sensitivity, epsilon) {
    flip <- sample(c(-1, 1), size=n, replace=TRUE)
    expon <- rexp(n=n, rate=(epsilon / sensitivity))
    return(flip * expon)
}


#' Random draw from Laplace distribution
#'
#' @param mu numeric, center of the distribution
#' @param b numeric, spread
#' @param n integer, number of draws
#' @return Random draws from Laplace distribution
#' @examples
#' rlap(size=1000)

rlap = function(mu=0, b=1, size=1) {
    p <- runif(size) - 0.5
    draws <- mu - b * sgn(p) * log(1 - 2 * abs(p))
    return(draws)
}


#' Probability density for Laplace distribution
#'
#' @param x numeric, value
#' @param mu numeric, center of the distribution
#' @param b numeric, spread
#' @return Density for elements of x
#' @examples
#' x <- seq(-3, 3, length.out=61)
#' dlap(x)

dlap <- function(x, mu=0, b=1) {
    dens <- 0.5 * b * exp(-1 * abs(x - mu) / b)
    return(dens)
}


#' Cumulative distribution function for Laplace distribution
#'
#' @param x numeric, value
#' @param mu numeric, center of the distribution
#' @param b numeric, spread
#' @return Probability less than or equal to x
#' @examples
#' x <- 0
#' plap(x)

plap <- function(x, mu=0, b=1) {
    cdf <- 0.5 + 0.5 * sgn(x - mu) * (1 - exp(-1 * (abs(x - mu) / b)))
    return(cdf)
}


#' Quantile function for Laplace distribution
#'
#' @param p Numeric, vector of probabilities
#' @param mu numeric, center of the distribution
#' @param b numeric, spread
#' @return Quantile function
#' @examples
#' probs <- c(0.05, 0.50, 0.95)
#' qlap(probs)

qlap <- function(p, mu=0, b=1) {
    q <- ifelse(p < 0.5, mu + b * log(2 * p), mu - b * log(2 - 2 * p))
    return(q)
}


#' Sign function
#'
#' @param x numeric, value or vector or values
#' @return The sign of passed values
#' @examples
#' sgn(rnrom(10))

sgn <- function(x) {
    return(ifelse(x < 0, -1, 1))
}


#' Utility function for checking that range is ordered pair
#'
#' @param range A vector, that ought to be an ordered pair
#' @return An ordered pair
#'
#' Checks if a supplied range is an ordered pair.  Coerces any vector of length two
#'   or greater into an ordered pair, and issues an error for shorter vectors.
#'
#' @examples
#' checkrange(1:3)
#' \dontrun{checkrange(1)}

checkrange <- function(rng) {
    if (NCOL(rng) > 1) {
        for (i in 1:nrow(rng)) {
            rng[i, ] <- sort(rng[i, ])
        }
    } else {
        if (length(rng) < 2) {
            stop("range argument in error: requires upper and lower values as vector of length 2.")
        } else if (length(rng) > 2) {
            warning("range argument supplied has more than two values.  Will proceed using min and max values as range.")
            rng <- c(min(rng), max(rng))
        } else {
            rng <- sort(rng)
        }
    }
	return(rng)
}


#' Utility function for checking that epsilon is acceptably defined
#'
#' @param epsilon A vector, that ought to be positive and length of 1
#' @return The supplied epsilon if acceptable, otherwise an error message interupts
#'
#' @examples
#' checkepsilon(0.1)
#' \dontrun{checkepsilon(-2)}
#' \dontrun{checkepsilon(c(0.1,0.5))}

checkepsilon = function(epsilon) {
	if (epsilon <= 0) {
		stop("Privacy parameter epsilon must be a value greater than zero.")
	}
	if (length(epsilon) > 1) {
		stop(paste("Privacy parameter epsilon must be a single value, but is currently a vector of length ", length(epsilon)))
	}
	return(epsilon)
}


#' Utility function for censoring data
#'
#' @param x A vector of numeric or categorial values to censor
#' @param var_type Character indicating the variable type
#' @param rng For numeric vectors, a vector (min, max) of the bounds of the range. For numeric matrices with nrow N and ncol P, a Px2 matrix of (min, max) bounds.
#' @param levels For categorical types, a vector containing the levels to be returned
#' @return Original vector with values outside the bounds censored to the bounds
#'
#' For numeric types, checks if x is in range = (min, max) and censors values to either min
#' or max if it is out of the range. For categorical types, values not in `levels` are coded NA.
#'
#' @examples
#' censordata(x=1:10, var_type='integer', range=c(2.5, 7))
#' censordata(x=c('a', 'b', 'c', 'd'), var_type='character', levels=c('a', 'b', 'c'))

censordata = function(x, var_type, rng=NULL, levels=NULL) {
    if (var_type %in% c('character', 'factor')) {
        if (is.null(levels)) {
            stop('`levels` are required for categorical types')
        }
        x <- factor(x, levels=levels, exclude=NULL)
    } else {
        if (is.null(rng)) {
            stop('range `rng` is required for numeric types')
        }
        if (NCOL(x) > 1) {
            for (j in 1:ncol(x)) {
                rng[j, ] <- checkrange(rng[j, ])
                x[, j][x[, j] < rng[j, 1]] <- rng[j, 1]
                x[, j][x[, j] > rng[j, 2]] <- rng[j, 2]
            }
        } else {
            rng <- checkrange(rng)
            x[x < rng[1]] <- rng[1]
            x[x > rng[2]] <- rng[2]
        }
    }
    return(x)
}


#' Utility function to check type of variable is within set of acceptable types 
#'
#' @param type Character specifying the type of the variable
#' @param in_types Vector of acceptable types 
#' @return The original character string indicating the variable type
#' 
#' Verifies that the variable is an element in the set of acceptable types
#' 
#' @examples 
#' check_variable_type(type='Numeric', in_types=c('Numeric', 'Factor'))

check_variable_type = function(type, in_types) { 
    if (!(type %in% in_types)) {
        stop(paste('Variable type', type, 'should be one of', paste(in_types, collapse = ', ')))
    } 
    return(type)
} 


#' Utility function to verify that a variable is dichotomous
#'
#' @param x Vector of values
#' @return Logical vector coded 0-1
#'
#' This function effectively allows the user to ask for any variable containing
#' at most two unique values to treat the variable as logical. If the variable
#' contains numeric values, the highest value is recoded 1 and the the lower
#' value is recoded 0. If the variable is categorical and contains only two unique
#' values, the least frequently observed is recoded 1.
#'
#' @examples
#' make_logical(sample(c('cat', 'dog'), size=8, replace=TRUE))
#' make_logical(sample(c(0, 1), size=8, replace=TRUE))
#' make_logical(sample(c(-6.87, 3.23), size=8, replace=TRUE)

make_logical <- function(x) {
    if (!length(unique(x)) <= 2) { # how to handle if contains 1 value only?
        stop('Variable has more than two values')
    }
    if (class(x) %in% c('character', 'factor')) {
        tab <- data.frame(table(x))
        min_level <- tab[which.min(tab[, 2]), 1]
        x <- ifelse(x == min_level, 1, 0)
    } else {
        x <- ifelse(x == max(x), 1, 0)
    }
    return(x)
}


#' Utility function to verify the type of histogram mechanism
#'
#' @param mechanism Character string specifying the mechanism
#' 
#' Verifies that the mechanism is one of `noisy`, `stability`, or `random` and returns 
#' the mechanism if so, else throws an error 
#' 
#' @examples 
#' check_histogram_mechanism('stability')

check_histogram_mechanism <- function(mechanism) { 
    if (!(is.null(mechanism)) && !(mechanism %in% c('noisy', 'stability', 'random'))) { 
        stop('`mechanism` must be one of `noisy`, `stability`, `random`')
    } 
    return(mechanism)
} 


#' Utility function to include NA level for categorical types when vector of bins
#' does not include all observed levels in the data vector.
#'
#' @param x Vector, categorical type
#' @param bins Vector, depositor-provided list of levels for which to count values

check_histogram_categorical <- function(x, bins) {
    x <- factor(x, levels=bins, exclude=NULL)
    return(x)
}


#' Utility function to check bins argument to histogram
#' 
#' @param n_bins Number of cells in which to tabulate values
#' @param n Number of observations
#' @return Number of bins
#' 
#' If number of bins is not provided, use the Sturges method

check_histogram_bins <- function(n_bins, n) {
    if (is.null(n_bins)) {
        n_bins <- ceiling(log2(n)) + 1
    } else if (n_bins < 2) {
        stop('`n_bins` must be at least 2')
    } else if (as.logical(n_bins %% 1)) {
        warning('non-integer value `n_bins` converted to next highest integer value')
        n_bins <- ceiling(n_bins)
    }
    return(n_bins)
}


#' Utility function to check sufficient n 
#' 

check_histogram_n <- function(accuracy, n, n_bins, epsilon, delta, alpha) { 
    cond1 <- (8 / accuracy) * (0.5 - log(delta) / epsilon)
    cond2 <- 4 * log(min(n_bins, (4 / accuracy)) / alpha) / (accuracy * epsilon)
    if (n < max(cond1, cond2, na.rm=TRUE)) { 
        return(FALSE)
        #stop('number of rows insufficient to provide privacy or accuracy with given parameters')
    } 
    return(TRUE)
} 


#' Utility function to match arguments of a function with list output of another function
#'
#' @param output List with output of a function
#' @param target.func Character name of the function with arguments that need to be filled by output
#' @return List of arguments and values needed for specification of \code{target.func}

getFuncArgs <- function(output, target.func) {
    spec <- list()
    for (element in names(output)) {
        if (element %in% names(formals(target.func))) {
            spec[[element]] <- output[[element]]
        }
    }
    return(spec)
}


#' Sweep operator (from Amelia)

amsweep <- function(g, m, reverse=FALSE) {
    if (identical(m, vector(mode='logical', length=length(m)))) {
        return(g)
    } else {
        p <- nrow(g)
        rowsm <- sum(m)
        if (rowsm == p) {
            h <- solve(g)
            h <- -h
        } else {
            kseq <- 1:p
            k <- kseq[m]
            kcompl <- kseq[-k]
            g11 <- g[k, k, drop=FALSE]
            g12 <- g[k, kcompl, drop=FALSE]
            g21 <- t(g12)
            g22 <- g[kcompl, kcompl, drop=FALSE]
            h11a <- try(solve(g11), silent=TRUE)
            if (inherits(h11a, "try-error")) {
                h11a <- mpinv(g11)
            }
            h11 <- as.matrix((-h11a))
            if (reverse) {sgn2 <- -1} else {sgn2 <- 1}
            h12 <- as.matrix(sgn2 * (h11a %*% g12))
            h21 <- as.matrix(t(h12))
            h22 <- g22 - (g21 %*% h11a %*% g12)
            hwo <- rbind(cbind(h11, h12), cbind(h21, h22))
            xordering <- c(k, kcompl)
            h <- matrix(0, p, p)
            h[xordering, xordering] <- hwo
        }
        return(h)
    }
}

#' Function to obtain indices in data frame for dependent & independent variables from a formula
#'
#' @param formula Formula
#' @param data Data frame, in this case being the data frame of a private covariance matrix
#' @param intercept Logical indicating whether the intercept is included
#' @return Named list with names corresponding to labels and locations (i.e., columns) for variables
#'  in the specification.
#'
#' @examples
#'
#' y <- rnorm(100) * 2
#' x <- (y + rnorm(100)) > 0
#' data <- data.frame(cbind(y, x))
#' f <- as.formula('y ~ x')
#' extract.indices(f, data, FALSE)

extract.indices <- function(formula, data, intercept) {
    t <- terms(formula, data=data)
    y.loc <- attr(t, 'response')
    x.loc <- which(names(data) %in% attr(t, 'term.labels'))
    x.label <- names(data)[x.loc]
    if (intercept) {
        intercept.loc <- which(names(data) == 'intercept')
        x.loc <- c(intercept.loc, x.loc)
        x.label <- append(x.label, 'Intercept', after=(intercept.loc - 1))
        if (intercept.loc <= y.loc) { y.loc <- y.loc + 1 }
    }
    return(list('y.loc' = y.loc,
                'x.loc' = x.loc,
                'x.label' = x.label))
}


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
#' @param terminal.level.index Index of the first terminal leaf node
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
