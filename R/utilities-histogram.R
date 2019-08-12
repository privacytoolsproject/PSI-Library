#' Get accuracy for a stable statistic (only used for histogram statistic)
#' 
#' Function to find the accuracy guarantee of a stable statistic release at a given epsilon 
#' value.
#' 
#' @param sensitivity the sensitivity of the statistic
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#'    Default to 2^-30.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @return Accuracy guarantee for statistic release given epsilon.

stability.getAccuracy <- function(sensitivity, epsilon, delta = 2^-30, alpha=0.05) {
    accuracy <- (sensitivity/epsilon) * log(2/(alpha*delta)) + 1
    return(accuracy)
}


#' Get epsilon for a stable statistic (only used for histogram statistic)
#' 
#' Function to find the epsilon value necessary to meet a desired level of 
#' accuracy for a stable statistic release.
#' 
#' @param sensitivity the sensitivity of the statistic
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#'    Default to 2^-30.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @return The scalar epsilon necessary to guarantee the needed accuracy.

stability.getEpsilon <- function(sensitivity, accuracy, delta = 2^-30, alpha=0.05) {
    epsilon <- sensitivity * log(2/(alpha*delta)) / (accuracy - 1)
    return(epsilon)
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
