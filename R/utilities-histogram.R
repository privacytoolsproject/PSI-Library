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

stabilityGetAccuracy <- function(sensitivity, epsilon, delta = 2^-30, alpha=0.05) {
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

stabilityGetEpsilon <- function(sensitivity, accuracy, delta = 2^-30, alpha=0.05) {
    epsilon <- sensitivity * log(2/(alpha*delta)) / (accuracy - 1)
    return(epsilon)
}


#' Utility function to include NA level for categorical types when vector of bins
#' does not include all observed levels in the data vector.
#'
#' @param x Vector, categorical type
#' @param bins Vector, depositor-provided list of levels for which to count values

histogramCategoricalBins <- function(x, bins) {
    x <- factor(x, levels=bins, exclude=NULL)
    return(x)
}


#' Check histogram bins argument
#' 
#' Utility function to check bins argument to histogram. If number of bins 
#'    is not provided, the Sturges method is used.
#' 
#' @param nBins The number of cells in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#'
#' @return Number of bins
#' @rdname checkHistogramNBins
checkHistogramNBins <- function(nBins, n) {
    if (!is.null(nBins)) { # nBins may be null, in which case we do not want to change it
        if (nBins < 2) {
            stop('number of bins must be at least 2')
        } else if (as.logical(nBins %% 1)) {
            warning('non-integer value for number of bins converted to next highest integer value')
            return(ceiling(nBins))
        } else {
            return(nBins) # if no issues with input value, return input value
        }
    } else {
        return(nBins) # return the input value if the input nBins is null
    }
}
