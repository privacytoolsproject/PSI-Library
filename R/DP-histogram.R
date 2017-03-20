#' Function to evaluate a histogram for a numeric variable
#'
#' @param x Vector of numeric values
#' @param var.levels Vector specifying the bins
#' @return List with fields `name` specifying the statistic and `stat` with the value of the statistic

dp.histogram.numeric <- function(x, var.levels, stability, n.bins, n) {
    values <- table(cut(x, breaks=var_levels, include.lowest=TRUE, right=TRUE))
    out <- list('name' = 'histogram',
                'stat' = values,
                'stability' = stability,
                'n.bins' = n.bins,
                'n' = n)
    return(out)
}


#' Function to evaluate a histogram for a categorical variable
#'
#' @param x Vector of categorical values
#' @param var.levels Vector specifying the bins
#' @return List with fields `name` specifying the statistic and `stat` with the value of the statistic

dp.histogram.categorical <- function(x, var.levels, stability, n.bins, n) {
    k <- length(var.levels)  # avoid unused argument error in mechanism
    values <- table(x, useNA='ifany')
    out <- list('name' = 'histogram',
                'stat' = values,
                'stability' = stability,
                'n.bins' = n.bins,
                'n' = n)
    return(out)
}


#' Function to evaluate a histogram and specify arguments
#'
#' @param x Vector of categorical or numeric values
#' @param var.type Character string indicating the variable type
#' @param stability Logical indicating if the stability mechanism is to be used
#' @param bins Vector of bins
#' @param n.bins Integer indicating the number of bins
#' @param n Integer indicating the number of observations in \code{x}
#' @return List with the true value of the statistic and arguments to be passed to other functions

dp.histogram <- function(x, var.type, stability, bins, n.bins, n) {
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
                'bins' = bins)
    return(out)
}


#' Release differentially private histogram
#'
#' @param x Vector, numeric or categorical
#' @param var_type String, specifies the type of x
#' @param range Tuple, range of x
#' @param n Integer, number of observations in x
#' @param epsilon Float, Epsilon value for differential privacy
#' @param delta Float, Delta value for differential privacy
#' @param beta Float, level of significance
#' @param bins Vector of bins for which values are counted
#' @param n_bins Integer, Number of cells in which to tabulate values in x, ignored if `var_type` \in {`factor`, `categorical`}
#' @param mechanism String, Differentially private mechanism
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
#' If the mechanism is not explicitly provided, use the mechanism with highest accuracy given the epsilon.
#'
#' @examples
#' # numeric types
#' x_num <- rnorm(100)
#' x_num_na <- x_num
#' x_num_na[sample(1:length(x_num_na), size=10, replace=FALSE)] <- NA
#' x_int <- as.integer(round(x_num * 20))
#' r_num <- histogram.release(x_num, var_type='numeric', range=c(-2, 2), n=100, epsilon=0.1)
#' r_num_na <- histogram.release(x_num_na, var_type='numeric', range=c(-2, 2), n=100, epsilon=0.1)
#' r_int <- histogram.release(x_int, var_type='integer', range=c(-40, 40), n=100, epsilon=0.1)
#' r_num_random <- histogram.release(x_num, var_type='numeric', range=c(-2, 2), n=100, epsilon=0.1, mechanism='random')
#' # accuracy is returning inf, which filters the entire release for stability histogram
#' r_num_stability <- histogram.release(x_num, var_type='numeric', range=c(-2, 2), n=100, epsilon=0.1, mechanism='stability')
#' r_num_noisy <- histogram.release(x_num, var_type='numeric', range=c(-2, 2), n=100, epsilon=0.1, mechanism='noisy')
#'
#' # categorical types
#' x_char <- c(rep('a', 40), rep('b', 25), rep('c', 15), rep('d', 12), rep('e', 5), rep('f', 2), rep('g', 1))
#' x_fac <- factor(x_char)
#' bins <- c('a', 'b', 'c', 'd', 'e')
#' r_char <- histogram.release(x_char, var_type='character', n=100, epsilon=0.1, bins=bins)
#' r_fac <- histogram.release(x_fac, var_type='factor', n=100, epsilon=0.1, bins=bins)

#histogram.release <- function(x, var_type, n, epsilon, delta=2^-30, beta=0.05, range=NULL, bins=NULL, n_bins=NULL, mechanism=NULL) {
#    var_type <- check_variable_type(var_type, in_types=c('numeric', 'integer', 'factor', 'character'))
#    if (is.null(mechanism)) {
#        mechs <- c('noisy', 'stability', 'random')
#        accuracies <- c(
#            histogram.getAccuracy('noisy', n_bins, n, epsilon),
#            histogram.getAccuracy('stability', n_bins, n, epsilon),
#            histogram.getAccuracy('random', n_bins, n, epsilon)
#        )
#        mechanism <- mechs[which.min(accuracies)]
#    }
#    mechanism <- check_histogram_mechanism(mechanism)
#    if (mechanism == 'random') {
#        release <- mechanism.histogram.random(x, var_type, epsilon, levels, n_bins)
#    } else {
#        if (var_type %in% c('factor', 'character')) {
#            release <- mechanism.laplace(dp.histogram.categorical, x, var_type, range, sensitivity=1, epsilon=epsilon, var.levels=bins)
#        } else {
#            n.bins <- check_histogram_bins(n_bins=n_bins, n=n)
#            var.levels <- seq(range[1], range[2], length.out=(n.bins + 1))
#            release <- mechanism.laplace(dp.histogram.numeric, x, var_type, range, sensitivity=1, epsilon=epsilon, var.levels=var.levels)
#        }
#        if (mechanism == 'noisy') {
#            release <- ifelse(release < 0, 0, round(release))
#        } else {
#            accuracy <- histogram.getAccuracy(mechanism, n_bins, n, epsilon)
#            if (check_histogram_n(accuracy, n, n_bins, epsilon, delta, beta)) {
#                a <- accuracy * n / 2
#                release <- release[release >= a]
#            }
#        }
#    }
#    return(release)
#}


# new release function
histogram.release <- function(x, var.type, n, epsilon, rng, bins=NULL, n.bins=NULL) {
    var.type <- check_variable_type(var.type, in_types=c('numeric', 'integer', 'factor', 'character'))
    if (var.type %in% c('numeric', 'integer')) {
        n.bins <- check_histogram_bins(n.bins, n)
        bins <- seq(rng[1], rng[2], length.out=(n.bins + 1))
    } else {
        n.bins <- length(bins)
    }
    release.noisy <- mechanism.laplace(
        fun=dp.histogram,
        x=x, var_type=var.type, rng=rng,
        sensitivity=1, epsilon=epsilon,
        stability=FALSE, bins=bins, n.bins=n.bins, n=n)
    release.stability <- mechanism.laplace(
        fun=dp.histogram,
        x=x, var_type=var.type, rng=rng,
        sensitivity=1, epsilon=epsilon,
        stability=TRUE, bins=bins, n.bins=n.bins, n=n)
    if (release.stability$accuracy < release.noisy$accuracy) {
        release <- release.stability
        if (check_histogram_n(release$accuracy, n, n.bins, epsilon, delta=2^-30, alpha=0.05)) {
            a <- release$accuracy * n / 2
            release$release <- release$release[release$release >= a]
        } else {
            release <- release.noisy
            release$release <- ifelse(release$release < 0, 0, round(release$release))
        }
    } else {
        release <- release.noisy
        release$release <- ifelse(release$release < 0, 0, round(release$release))
    }
    return(release)
}


# new accuracy function
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


#' Accuracy of release 
#' 
#' @param mechanism Specify the histogram mechanism 
#' @param n_bins Number of cells in which to tabulate values 
#' @param epsilon Epsilon value for differential privacy
#' @param delta Delta value for differential privacy
#' @param beta

#histogram.getAccuracy <- function(mechanism, n_bins, n, epsilon, delta=2^-30, beta=0.05, error=1e-9) { 
#    mechanism <- check_histogram_mechanism(mechanism)
#
#    if (mechanism == 'noisy') { 
#        acc <- -2 * log(1 - (1 - beta)^(1 / n_bins)) / (n * epsilon)
#    } else if (mechanism == 'stability') { 
#        lo <- 0
#        hi <- 1
#        eval <- beta + error
#        while ((eval <= beta - error) || (eval > beta)) { 
#            acc <- (hi + lo) / 2
#            eval <- min((4 / acc), n_bins) * exp(-acc * n * epsilon / 4)
#            if (eval < beta) { 
#                hi <- acc
#            } else { 
#                lo <- acc
#            } 
#            if (hi - lo <= 0) { 
#                return(Inf)
#            } 
#        } 
#        acc <- max(acc, (8 / n) * (0.5 - log(delta) / epsilon))
#    } else { 
#        acc <- Inf
#    } 
#    return(acc) 
#} 


# new get parameters fn
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



#' Privacy parameters
#'
#' @param mechanism Specify the histogram mechanism
#' @param n_bins Number of cells in which to tabulate values
#' @param accuracy
#' @param delta Delta value for differential privacy
#' @param beta
#' @return epsilon Differential privacy parameter

#histogram.getParameters <- function(mechanism, n_bins, n, accuracy, delta, beta=0.05, error=1e-9) {
#    mechanism <- check_histogram_mechanism(mechanism)
#    if (mechanism == 'noisy') {
#        eps <- -2 * log(1 - (1 - beta)^(1 / n_bins)) / (n * accuracy)
#    } else if (mechanism == 'stability') {
#        lo <- 0
#        hi <- 1
#        eval <- n + 1
#        while ((eval <= n * (1 - error)) || (eval > n)) {
#            eps <- (hi + lo) / 2
#            eval <- 8 / accuracy * (0.5 - log(delta) / eps)
#            if (eval < n) {
#                hi <- eps
#            } else {
#                lo <- eps
#            }
#            if (hi - lo <= 0) {
#                return(Inf)
#            }
#        }
#        eps <- max(eps, 4 * log(min(n_bins, 4 / accuracy) / beta) / (accuracy * n))
#    } else {
#        eps <- log((n_bins + 1) / (n_bins - 1))
#    }
#    return(eps)
#}


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
