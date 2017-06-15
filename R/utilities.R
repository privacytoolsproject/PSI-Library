#' Random draw from Laplace distribution
#'
#' @param n integer, number of draws
#' @param sensitivity Numeric, the sensitivity of the estimate
#' @param epsilon Numeric, epsilon parameter for differential privacy
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
#' @param size number of observations
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
#' @param rng A vector, that ought to be an ordered pair
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
#' @param epsilon Numeric, epsilon parameter for differential privacy
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
#' censordata(x=1:10, var_type='integer', rng=c(2.5, 7))
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
#' @return something here
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
#' @return something here

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
#' @param accuracy The accuracy we need to guarantee
#' @param n Number of observations
#' @param n_bins Number of cells in which to tabulate values
#' @param epsilon Numeric, epsilon parameter for differential privacy
#' @param delta something here
#' @param alpha something here
#' @return something here

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
#' 
#' @param g something here
#' @param m something here
#' @param reverse something here
#' @return something here

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
                h11a <- MASS::ginv(g11)
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
