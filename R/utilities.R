#' Differentially Private Uniform Draw 
#' 
#' Draw cryptographically secure random variates from a uniform distribution.
#'
#' @param n An integer giving number of variates needed.
#' @param seed An integer indicating a seed for R's PNRG, defaults to \code{NULL}.
#' 
#' @return Random numeric vector of length \code{n} containing values between
#'    zero and one.
#'
#' Draws secure random variates from the uniform distribution through \code{openssl}.
#' If a seed is provided, the \code{runif} function is used to draw the random variates.
#' @examples
#' 
#' uniform_secure <- dpUnif(n=1000)
#' uniform_repeatable <- dpUnif(n=1, seed=75436)
#' @seealso \code{\link{dpNoise}}
#' @rdname dpUnif
#' @export
dpUnif <- function(n, seed=NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
        return(runif(n))
    }
    return(openssl::rand_num(n))
}


#' Differentially Private Noise Generator
#' 
#' Compile noise from a cryptographically secure random variates to achieve
#'    differentially private statistics.
#'
#' @param n An integer giving number of variates needed.
#' @param scale Numeric, the scale for the distribution.
#' @param dist A character specifying the distribution from which to draw the 
#'    noise.
#' @param shape An integer giving the shape parameter for the gamma
#'    distribution. Default to \code{NULL}.
#' @param seed An integer indicating a seed for R's PNRG, defaults 
#'    to \code{NULL}.
#'
#' @return Cryptographically secure noise vector or matrix.
#' @examples
#'
#' laplace_noise <- dpNoise(n=1000, scale=1, dist='laplace')
#' gaussian_noise <- dpNoise(n=1000, scale=1, dist='gaussian')
#' laplace_noise_repeatable <- dpNoise(n=1, scale=1, dist='laplace', seed=96845)
#' @seealso \code{\link{dpUnif}}
#' @rdname dpNoise
#' @export
dpNoise <- function(n, scale, dist, shape=NULL, seed=NULL) {
    u <- dpUnif(n, seed)
    if (dist == 'laplace') {
        return(qlap(u, b=scale))
    } else if (dist == 'gaussian') {
        return(qnorm(u, sd=scale))
    } else if (dist == 'gamma') {
        return(qgamma(u, scale=scale, shape=shape))
    } else {
        stop(sprintf('Distribution "%s" not understood', dist))
    }
}


#' Fill missing values
#'
#' Impute uniformly in the range of the provided variable
#'
#' @param x Vector with missing values to impute
#' @param var.type Character specifying the variable type
#' @param lower Numeric lower bound, default NULL
#' @param upper Numeric upper bound, default NULL
#' @param categories Set of possible categories from which to choose,
#'      default NULL
#' @param seed Integer indicating a seed for R's PRNG, default NULL
#' @return Vector \code{x} with missing values imputed
#' @examples
#'
#' # numeric example
#' y <- rnorm(100)
#' y[sample(1:100, size=10)] <- NA
#' y_imputed <- fillMissing(x=y, var.type='numeric', lower=-1, upper=1)
#'
#' # categorical example
#' cats <- as.factor(c('dog', 'zebra', 'bird', 'hippo'))
#' s <- sample(cats, size=100, replace=TRUE)
#' s[sample(1:100, size=10)] <- NA
#' s_imputed <- fillMissing(x=s, var.type='factor', categories=cats)
#'
#' @seealso \code{\link{dpUnif}}
#' @rdname dpNoise
#' @export
fillMissing <- function(x, var.type, lower=NULL, upper=NULL, categories=NULL, seed=NULL) {
    miss_idx <- is.na(x)
    if (sum(miss_idx) == 0) { return(x) }
    u <- dpUnif(sum(miss_idx), seed)
    if (var.type %in% c('character', 'factor')) {
        lower <- 1
        upper <- length(categories)
    }
    out <- u * (upper - lower) + lower
    if (var.type %in% c('logical', 'integer')) {
        out <- as.integer(out)
    } else if (var.type %in% c('character', 'factor')) {
        out <- categories[as.integer(out)]
    }
    x[miss_idx] <- out
    return(x)
}


#' Random draw from Laplace distribution
#'
#' @param mu numeric, center of the distribution
#' @param b numeric, spread
#' @param size integer, number of draws
#' 
#' @return Random draws from Laplace distribution
#' @examples
#' 
#' rlap(size=1000)
#' @export
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
#' 
#' @return Density for elements of x
#' @examples
#' 
#' x <- seq(-3, 3, length.out=61)
#' dlap(x)
#' @export
dlap <- function(x, mu=0, b=1) {
    dens <- 0.5 * b * exp(-1 * abs(x - mu) / b)
    return(dens)
}


#' LaPlace Cumulative Distribution Function
#' 
#' Determines the probability a draw from a LaPlace distribution is less than 
#'    or equal to the specified value.
#'
#' @param x Numeric, the value(s) at which the user wants to know the CDF height.
#' @param mu Numeric, the center of the LaPlace distribution, defaults to 0.
#' @param b Numeric, the spread of the LaPlace distribution, defaults to 1.
#' 
#' @return Probability the LaPlace draw is less than or equal to \code{x}.
#' @examples
#' 
#' x <- 0
#' plap(x)
#' @rdname plap
#' @export
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
#' @export
qlap <- function(p, mu=0, b=1) {
    q <- ifelse(p < 0.5, mu + b * log(2 * p), mu - b * log(2 - 2 * p))
    return(q)
}


#' Sign function
#' 
#' Function to determine what the sign of the passed values should be.
#'
#' @param x numeric, value or vector or values
#' @return The sign of passed values
#' @examples
#' sgn(rnorm(10))
#' @export
sgn <- function(x) {
    return(ifelse(x < 0, -1, 1))
}


#' Range Parameter Check
#' 
#' Checks if a supplied range is an ordered pair. Coerces any vector of length 
#'    two or greater into an ordered pair, and issues an error for
#'    shorter vectors.
#'    
#' @param rng A numeric vector of length two, that ought to be an 
#'    ordered pair.
#' 
#' @return An ordered pair.
#' @examples
#'
#' checkrange(c(1,3))
#' checkrange(1:3)
#' \dontrun{checkrange(1)}
#' @rdname checkrange
#' @export
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


#' Epsilon Parameter Check
#' 
#' Utility function for checking that epsilon is acceptably defined.
#'
#' @param epsilon A vector, that ought to be positive and of length one.
#' 
#' @return The supplied epsilon if acceptable, otherwise an error 
#'    message interupts.
#'
#' @examples
#' 
#' checkepsilon(0.1)
#' \dontrun{checkepsilon(-2)}
#' \dontrun{checkepsilon(c(0.1,0.5))}
#' @rdname checkepsilon
#' @export
checkepsilon <- function(epsilon) {
	if (epsilon <= 0) {
		stop("Privacy parameter epsilon must be a value greater than zero.")
	}
	if (length(epsilon) > 1) {
		stop(paste("Privacy parameter epsilon must be a single value, but is currently a vector of length ", length(epsilon)))
	}
	return(epsilon)
}


#' Censoring data
#' 
#' For numeric types, checks if x is in rng = (min, max) and censors values to 
#'    either min or max if it is out of the range. For categorical types, 
#'    values not in `levels` are coded NA.
#'
#' @param x A vector of numeric or categorial values to censor.
#' @param var_type Character indicating the variable type of \code{x}.
#' @param rng For numeric vectors, a vector (min, max) of the bounds of the 
#'    range. For numeric matrices with nrow N and ncol P, a Px2 matrix of 
#'    (min, max) bounds.
#' @param levels For categorical types, a vector containing the levels to 
#'    be returned.
#' 
#' @return Original vector with values outside the bounds censored to the bounds.
#' @examples
#' 
#' censordata(x=1:10, var_type='integer', rng=c(2.5, 7))
#' censordata(x=c('a', 'b', 'c', 'd'), var_type='character', levels=c('a', 'b', 'c'))
#' @rdname censordata
#' @export
censordata <- function(x, var_type, rng=NULL, levels=NULL) {
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


#' Checking variable types
#' 
#' Verifies that the variable is an element in the set of acceptable types.
#' 
#' @param type A character specifying the type of the variable.
#' @param in_types A vector of acceptable types of variables.
#' 
#' @return The original character string indicating the variable type.
#' @examples 
#' 
#' check_variable_type(type='Numeric', in_types=c('Numeric', 'Factor'))
#' @rdname check_variable_type
#' @export
check_variable_type <- function(type, in_types) { 
    if (!(type %in% in_types)) {
        stop(paste('Variable type', type, 'should be one of', paste(in_types, collapse = ', ')))
    } 
    return(type)
} 


#' Logical variable check
#' 
#' Utility function to verify that a variable is dichotomous.
#'
#' This function effectively allows the user to ask for any variable containing
#' at most two unique values to treat the variable as logical. If the variable
#' contains numeric values, the highest value is recoded 1 and the the lower
#' value is recoded 0. If the variable is categorical and contains only two unique
#' values, the least frequently observed is recoded 1.
#'
#' @param x Vector containing two unique types of values.
#' 
#' @return Logical form of \code{x} coded 0-1.
#' @examples
#' 
#' make_logical(sample(c('cat', 'dog'), size=8, replace=TRUE))
#' make_logical(sample(c(0, 1), size=8, replace=TRUE))
#' make_logical(sample(c(-6.87, 3.23), size=8, replace=TRUE))
#' @rdname make_logical
#' @export
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
#' 
#' check_histogram_mechanism('stability')
#' @export
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


#' Check histogram bins argument
#' 
#' Utility function to check bins argument to histogram. If number of bins 
#'    is not provided, the Sturges method is used.
#' 
#' @param n_bins The number of cells in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#'
#' @return Number of bins
#' @rdname check_histogram_bins
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


#' Histogram N check
#' 
#' Utility function to check sufficient N in data.
#' 
#' @param accuracy A numeric vector of length one representing the accuracy of 
#'    the noisy estimate
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#' @param n_bins A numeric vector of length one specifying the number of cells 
#'    in which to tabulate values.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#' @param alpha A numeric vector of length one specifying the numeric 
#'    statistical significance level.
#'    
#' @return A logical vector indicating whether the number of observations is 
#'    sufficient to provide desired privacy and accuracy with the given
#'    parameters.
#' @rdname check_histogram_n
check_histogram_n <- function(accuracy, n, n_bins, epsilon, delta, alpha) { 
    cond1 <- (8 / accuracy) * (0.5 - log(delta) / epsilon)
    cond2 <- 4 * log(min(n_bins, (4 / accuracy)) / alpha) / (accuracy * epsilon)
    if (n < max(cond1, cond2, na.rm=TRUE)) { 
        return(FALSE)
        #stop('number of rows insufficient to provide privacy or accuracy with given parameters')
    } 
    return(TRUE)
}


#' Extract function arguments
#' 
#' Utility function to match arguments of a function with list output of another function.
#'
#' @param output A list with the output of a function.
#' @param target.func A character vector containing the name of the function 
#'    with arguments that need to be filled by output.
#'    
#' @return List of arguments and values needed for specification of \code{target.func}
#' @rdname getFuncArgs
getFuncArgs <- function(output, target.func) {
    spec <- list()
    for (element in names(output)) {
        if (element %in% names(formals(target.func))) {
            spec[[element]] <- output[[element]]
        }
    }
    return(spec)
}


#' Linear regression on covariance matrix
#' 
#' Function to extract regression coefficients from the covariance matrix via 
#'    the sweep operator.
#'
#' @param formula An object of class 'formula' containing the desired 
#'    regression formula.
#' @param release A numeric private release of the covariance matrix.
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#' @param intercept A logical vector indicating whether the intercept is 
#'    included in \code{formula}.
#'    
#' @return A numeric vector of regression coefficients corresponding 
#'    to \code{formula}.
#' @rdname linear.reg
linear.reg <- function(formula, release, n, intercept) {
  if (!all(eigen(release)$values > 0)) {  # could do is.positive.definite() but that requires matrixcalc package
    coefs <- "The input matrix is not invertible"
    return(coefs)
  } else {
    xy.locs <- extract.indices(formula, release, intercept)
    x.loc <- xy.locs$x.loc
    y.loc <- xy.locs$y.loc
    loc.vec <- rep(TRUE, (length(x.loc) + 1))
    loc.vec[y.loc] <- FALSE
    sweep <- amsweep((as.matrix(release) / n), loc.vec)
    coefs <- sweep[y.loc, x.loc]
    se <- sqrt(sweep[y.loc, y.loc] * diag(solve(release[x.loc, x.loc])))
    coefs <- data.frame(cbind(coefs, se))
    coefs <- format(round(coefs, 5), nsmall=5)
    rownames(coefs) <- xy.locs$x.label
    names(coefs) <- c('Estimate', 'Std. Error')
    return(coefs)
  }
}


#' Moore Penrose Inverse Function ###Must assign authorship to this###
#' 
#' Generate the Moore-Penrose pseudoinverse matrix of \code{X}.
#' 
#' @param X A numeric, symmetric covariance matrix.
#' @param tol Convergence requirement. 

mpinv <- function(X, tol = sqrt(.Machine$double.eps)) {
  ## Moore-Penrose Inverse function (aka Generalized Inverse)
  ##   X:    symmetric matrix
  ##   tol:  convergence requirement
  s <- svd(X)
  e <- s$d
  e[e > tol] <- 1/e[e > tol]
  s$v %*% diag(e,nrow=length(e)) %*% t(s$u)
}


#' Sweep operator ###Check if we need to assign authorship to this###
#' 
#' Sweeps a covariance matrix to extract regression coefficients.
#' 
#' @param g Each unit of a numeric, symmetric covariance matrix divided by
#'    the number of observations in the data the covariance matrix was 
#'    derived from.
#' @param m A logical vector of length equal to the number of rows in \code{g}
#'    in which the \code{TRUE} values correspond to the x values in the matrix
#'    and the \code{FALSE} values correspond to the y value in the matrix.
#' @param reverse Logical vector specifying whether the sign of the matrix 
#'    should be flipped. Default to \code{FALSE}.
#' 
#' @return The coefficients from \code{g}.
#' @rdname amsweep
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


#' Extract regression indices
#' 
#' Function to obtain indices in data frame for dependent & independent variables from a formula.
#'
#' @param formula An object of class 'formula' containing the desired 
#'    regression formula.
#' @param data A numeric data frame with at least two columns.
#' @param intercept A logical vector indicating whether the intercept is 
#'    included in \code{formula}.
#' 
#' @return A named list with names corresponding to labels and locations 
#'    (i.e., columns) for variables in the specification.
#' @examples
#'
#' y <- rnorm(100) * 2
#' x <- (y + rnorm(100)) > 0
#' data <- data.frame(cbind(y, x))
#' f <- as.formula('y ~ x')
#' extract.indices(f, data, FALSE)
#' @rdname extract.indices
#' @export
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

#' Function to convert factor variables to binary indicators
#'
#' @param df Data frame
#' @return List with data frame with factor columns converted to dummy indicators and the names of
#'      the columns of the transformed data frame.
#'
#' For each factor variable in the data frame, a binary indicator is generated for (k - 1) of its
#' k levels. The first level is dropped. The original factor variable is dropped. The names of the
#' binary indicators are the result of combining the name of the original factor variable and the
#' level represented by the indicator.

makeDummies <- function(df) {
    factors <- sapply(df, class) == 'factor'
    factors <- names(df)[factors]
    for (col in factors) {
        col.levels <- levels(df[, col])
        dummy.list <- lapply(col.levels, function(x) as.numeric(df[, col] == x))
        dummy.df <- data.frame(dummy.list)
        names(dummy.df) <- sapply(col.levels, function(x) paste(col, x, sep='_'))
        df <- cbind(df, dummy.df[, 2:length(col.levels)])  # drop first level
    }
    df <- df[, !(names(df) %in% factors)]
    return(list('data' = df, 'names' = names(df)))
}


#' Function to evaluate the p-norm of vectors in a matrix
#'
#' @param X matrix of numeric values
#' @param p The p in p-norm
#' @param margin The subscripts over which the norm is applied, where 1 gives row norms
#'      and 2 gives the column norms
#' @return A vector of norms

vectorNorm <- function(X, p=1, margin=1) {
    fun <- function(vec, p) {
        vec.norm <- sum(abs(vec)^p)^(1 / p)
        return(vec.norm)
    }
    if (!is.null(dim(X))) {
        p.norm <- apply(X, margin, fun, p)
    } else {
        p.norm <- fun(X, p)
    }
    return(p.norm)
}


#' Function to map rows of a numeric matrix to the unit ball
#'
#' This function will ensure that the row having the largest p-norm will be located on the unit ball
#'
#' @param X Numeric matrix
#' @param p The p in p-norm
#' @return List that includes the transformed matrix and the value corresponding to the maximum 
#'      observed p-norm

mapMatrixUnit <- function(X, p=1) {
    max.norm <- max(vectorNorm(X, p=p))
    normed.matrix <- X / max.norm
    return(list('matrix' = normed.matrix, 'max.norm' = max.norm))
}


#' Function to trim lower and upper regions of a vector of values
#'
#' @param vec Numeric vector
#' @param alpha Numeric proportion of vector to be trimmed, specifically the 
#'      least and greatest \code{alpha / 2} are trimmed
#' @return Trimmed vector

trimVector <- function(vec, alpha) {
    alpha <- alpha / 2
    lower <- quantile(vec, probs=alpha)
    upper <- quantile(vec, probs=(1 - alpha))
    trimmed <- vec[vec >= lower & vec <= upper]
    return(trimmed)
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
