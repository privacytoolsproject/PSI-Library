#' Heavyhitters accuracy
#' 
#' Get the accuracy of heavyhitters statistic for a given value of epsilon.
#' The accuracy is interpreted as follows: The alpha value returned means that with
#' probability 1 - beta, simultaneously for all t with min <= t <= max, the algorithm's
#' estimate of the count in [min, t] is within alpha of the true value.
#'
#' @param gap The difference between the \code{k} category and the 
#'    \code{k+1} category.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param delta A numeric vector representing the probability of an arbitrary
#'    leakage of information from the data. Should be of length one 
#'    and should be a very small value.
#'    
#' @return Accuracy guarantee for heavyhitters release given epsilon.
#' @author Victor Balcer
#' @rdname heavyhittersGetAccuracy

heavyhittersGetAccuracy <- function(gap, epsilon, delta) {
    accuracy <- exp(-epsilon * gap / 2) / delta
    return(accuracy)
}


#' Heavyhitters epsilon
#' 
#' Get the epsilon value necessary to guarantee a desired level of accuracy of a heavyhitters release
#'
#' @param gap The difference between the \code{k} category and the 
#'    \code{k+1} category.
#' @param delta A numeric vector representing the probability of an arbitrary
#'    leakage of information from the data. Should be of length one 
#'    and should be a very small value.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @return The epsilon value necessary to gaurantee the given accuracy.
#' @author Victor Balcer
#' @rdname heavyhittersGetParameters

heavyhittersGetParameters <- function(gap, delta, alpha=0.05) {
  epsilon <- -2 * log(alpha * delta) / gap
  return(epsilon)
}


#' Heavyhitters function
#' 
#' Function to evaluate most common categorical values
#' 
#' @param x Vector of categorical values
#' @return Sorted histogram with counts for each level

funHeavy <- function(x) {
    histogram <- table(x, useNA='ifany')
    histogram <- sort(-histogram) * -1
    return(histogram)
}


#' Differentially private heavy hitters
#'
#' @param mechanism Character, the mechanism used to perturb the output.
#' @param varType Character, the variable type, one of \code{'character'} or
#'    \code{'factor'}.
#' @param variable Character, the variable name in the data frame.
#' @param n Integer, the number of observations.
#' @param epsilon Numeric, the privacy loss parameter.
#' @param k Integer, the number of bins to release.
#' @param bins Character, the available bins, or the levels of the categorical variable.
#' @param alpha Numeric, level of statistical significance, default 0.05.
#' @param delta Numeric, probability of privacy loss beyond \code{epsilon}.
#'
#' @import methods
#' @export dpHeavyHitters
#' @exportClass dpHeavyHitters
#'
#' @include mechanism.R
#' @include mechanism-exponential.R

dpHeavyHitters <- setRefClass(
    Class = 'dpHeavyHitters',
    contains = 'mechanismExponential'
)

dpHeavyHitters$methods(
    initialize = function(mechanism, varType, variable, n, epsilon, 
                          k, bins, alpha=0.05, delta=NULL) {
        .self$name <- 'Differentially private heavy hitters'
        .self$mechanism <- mechanism
        .self$variable <- variable
        .self$varType <- checkVariableType(varType, inTypes=c('character', 'factor'))
        .self$n <- checkN(n)
        .self$epsilon <- epsilon
        .self$k <- k
        .self$bins <- bins
        .self$alpha <- alpha
        .self$delta <- checkDelta(.self$mechanism, delta)
})

dpHeavyHitters$methods(
    release = function(data) {
        x <- data[, variable]
        .self$result <- export(mechanism)$evaluate(funHeavy, x, 2)
s})

dpHeavyHitters$methods(
    postProcess = function(out, gap) { #is gap allowed to be public?? it isn't released from the exponential mechanism.
        out$variable <- variable
        out$accuracy <- heavyhittersGetAccuracy(gap, epsilon, delta)
        out$epsilon <- heavyhittersGetParameters(gap, delta, alpha)
        out$delta <- .self$delta
})
