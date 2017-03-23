#' Function to evaluate a histogram and specify arguments
#'
#' @param x Vector of categorical or numeric values
#' @param var.type Character string indicating the variable type
#' @param stability Logical indicating if the stability mechanism is to be used
#' @param bins Vector of bins
#' @param n.bins Integer indicating the number of bins
#' @param n Integer indicating the number of observations in \code{x}
#' @return List with the true value of the statistic and arguments to be passed to other functions

dp.histogram <- function(x, var.type, stability, bins, n.bins, n, sensitivity, epsilon) {
    if (var.type %in% c('numeric', 'integer')) {
        values <- table(cut(x, breaks=bins, include.lowest=TRUE, right=TRUE))
    } else {
        values <- table(x, useNA='ifany')
    }
    out <- list('name' = 'histogram',
                'stat' = values,
                'stability' = stability,
                'n.bins' = n.bins,
                'n' = n,
                'bins' = bins,
                'sensitivity' = sensitivity,
                'epsilon' = epsilon)
    return(out)
}


#' Release differentially private histogram
#'
#' @param x Vector, numeric or categorical
#' @param var.type String, specifies the type of x
#' @param n Integer, number of observations in x
#' @param epsilon Float, Epsilon value for differential privacy
#' @param rng Tuple, range of x, required for numeric types
#' @param bins Vector of bins for which values are counted, required for categorical types
#' @param n.bins Integer, Number of cells in which to tabulate values in x, ignored if \code{var.type \%in\% c('factor', 'categorical')}
#'
#' If the variable is categorical, bins are assumed to be provided by the depositor, and these bin values
#' used to construct the table. The vector is pre-processed so that observed levels not specified in these
#' bins are recoded to `NA`. Thus, any observed levels not specified in the `bins` argument show up as `NA`
#' in the output table.
#'
#' If the variable is numeric, the number of bins `n_bins` is set by the user optionally, else the Sturges
#' method is used to select the number of bins given the number of observations `n`. The bins are then
#' constructed to be equal intervals between the provided range.
#'
#' Uses the Laplace mechanism. If the stability mechanism improves accuracy, its value is used.
#'
#' @examples
#' # numeric types
#' x_num <- rnorm(100)
#' x_num_na <- x_num
#' x_num_na[sample(1:length(x_num_na), size=10, replace=FALSE)] <- NA
#' x_int <- as.integer(round(x_num * 20))
#' r_num <- histogram.release(x_num, var.type='numeric', range=c(-2, 2), n=100, epsilon=0.1)
#' r_num_na <- histogram.release(x_num_na, var.type='numeric', range=c(-2, 2), n=100, epsilon=0.1)
#' r_int <- histogram.release(x_int, var.type='integer', range=c(-40, 40), n=100, epsilon=0.1)
#' r_num_random <- histogram.release(x_num, var.type='numeric', range=c(-2, 2), n=100, epsilon=0.1, mechanism='random')
#' # accuracy is returning inf, which filters the entire release for stability histogram
#' r_num_stability <- histogram.release(x_num, var.type='numeric', range=c(-2, 2), n=100, epsilon=0.1, mechanism='stability')
#' r_num_noisy <- histogram.release(x_num, var.type='numeric', range=c(-2, 2), n=100, epsilon=0.1, mechanism='noisy')
#'
#' # categorical types
#' x_char <- c(rep('a', 40), rep('b', 25), rep('c', 15), rep('d', 12), rep('e', 5), rep('f', 2), rep('g', 1))
#' x_fac <- factor(x_char)
#' bins <- c('a', 'b', 'c', 'd', 'e')
#' r_char <- histogram.release(x_char, var.type='character', n=100, epsilon=0.1, bins=bins)
#' r_fac <- histogram.release(x_fac, var.type='factor', n=100, epsilon=0.1, bins=bins)

histogram.release <- function(x, var.type, n, epsilon, rng=NULL, bins=NULL, n.bins=NULL) {
    var.type <- check_variable_type(var.type, in_types=c('numeric', 'integer', 'factor', 'character'))
    if (var.type %in% c('numeric', 'integer')) {
        n.bins <- check_histogram_bins(n.bins, n)
        bins <- seq(rng[1], rng[2], length.out=(n.bins + 1))
    } else {
        n.bins <- length(bins)
    }
    release.noisy <- mechanism.laplace(
        fun=dp.histogram,
        x=x, var.type=var.type, rng=rng,
        sensitivity=1, epsilon=epsilon,
        stability=FALSE, bins=bins, n.bins=n.bins, n=n)
    release.stability <- mechanism.laplace(
        fun=dp.histogram,
        x=x, var.type=var.type, rng=rng,
        sensitivity=1, epsilon=epsilon,
        stability=TRUE, bins=bins, n.bins=n.bins, n=n)
    stability.accurate <- release.stability$accuracy < release.noisy$accuracy
    stability.check <- check_histogram_n(release.stability$accuracy, n, n.bins, epsilon, delta=2^-30, alpha=0.05)
    if (stability.accurate && stability.check) {
        release <- release.stability
        a <- release$accuracy * n / 2
        release$release <- release$release[release$release >= a]
    } else {
        release <- release.noisy
        release$release <- ifelse(release$release < 0, 0, round(release$release))
    }
    return(release)
}


#' Accuracy of release
#'
#' @param n.bins Integer indicating number of cells in which to tabulate values
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon value for differential privacy
#' @param stability Logical indicating whether stability mechanism is used
#' @param delta Numeric delta value for differential privacy, fixed
#' @param alpha Numeric statistical significance level, fixed
#' @param error Numeric, fixed
#' @return Accuracy

histogram.getAccuracy <- function(n.bins, n, epsilon, stability, delta=2^-30, alpha=0.05, error=1e-9) {
    if (stability) {
        lo <- 0
        hi <- 1
        eval <- alpha + error
        while ((eval <= alpha - error) || (eval > alpha)) {
            acc <- (hi + lo) / 2
            eval <- min((4 / acc), n.bins) * exp(-acc * n * epsilon / 4)
            ifelse(eval < alpha, (hi <- acc), (lo <- acc))
            if (hi - lo <= 0) { return(Inf) }
        }
        acc <- max(acc, (8 / n) * (0.5 - log(delta) / epsilon))
    } else {
        acc <- -2 * log(1 - (1 - alpha)^(1 / n.bins)) / (n * epsilon)
    }
    return(acc)
}


#' Privacy parameters
#'
#' @param n.bins Integer indicating number of cells in which to tabulate values
#' @param n Integer indicating number of observations
#' @param accuracy Numeric
#' @param delta Numeric delta value for differential privacy, fixed
#' @param alpha Numeric statistical significance level, fixed
#' @param error Numeric, fixed
#' @return Differential privacy parameter epsilon

histogram.getParameters <- function(n.bins, n, accuracy, stability, delta=2^-30, alpha=0.05, error=1e-9) {
    if (stability) {
        lo <- 0
        hi <- 1
        eval <- n + 1
        while ((eval <= n * (1 - error)) || (eval > n)) {
            eps <- (hi + lo) / 2
            eval <- 8 / accuracy * (0.5 - log(delta) / eps)
            ifelse(eval < n, (hi <- eps), (lo <- eps))
            if (hi - lo <= 0) { return(Inf) }
        }
        eps <- max(eps, 4 * log(min(n.bins, 4 / accuracy) / alpha) / (accuracy * n))
    } else {
        eps <- -2 * log(1 - (1 - alpha)^(1 / n.bins)) / (n * accuracy)
    }
    return(eps)
}


#' Confidence interval
#'

histogram.getCI <- function(release, n.bins, n, accuracy) {
    accxn <- accuracy * n
    out <- list()
    for (k in 1:n.bins) {
        bin.count <- release[k]
        if (bin.count == 0) {
            out[[k]] <- c(0, accxn)
        } else {
            out[[k]] <- c(max(0, bin.count - accxn), accxn + bin.count)
        }
    }
    out <- data.frame(do.call(rbind, out))
    names(out) <- c('lower', 'upper')
    return(out)
}


#' JSON doc for histogram
#'
#' @return JSON for histogram function

histogram.getJSON <- function(output.json=TRUE) {
    out <- list()
    out$statistic <- 'Histogram'
    out$description <- 'Differentially Private Histogram'
    out$mechanisms <- c('noisy', 'random', 'stability')
    out$variableTypes <- list('numeric' = list(), 'categorical' = list())
    out$variableTypes$numeric$rTypes <- c('numeric', 'integer')
    out$variableTypes$numeric$fields <- list(
            'n' = 'Number of observations',
            'range' = 'Ordered pair indicating effective lower and upper bounds',
            'n_bins' = 'Number of cells in output (optional, default Sturges method)'
        )
    out$variableTypes$categorical$rTypes <- c('character', 'factor')
    out$variableTypes$categorical$fields <- list(
            'n' = 'Number of observations',
            'bins' = 'Vector indicating levels for which to produce frequencies'
        )
    if (output.json) {
        out <- jsonlite::toJSON(out, pretty=TRUE)
    }
    return(out)
}
