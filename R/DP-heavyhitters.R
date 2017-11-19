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
#' @rdname heavyhitters.getAccuracy
heavyhitters.getAccuracy <- function(gap, epsilon, delta) {
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
#' @rdname heavyhitters.getParameters
heavyhitters.getParameters <- function(gap, delta, alpha=0.05) {
  epsilon <- -2 * log(alpha * delta) / gap
  return(epsilon)
}


#' Heavyhitters function
#' 
#' Function to evaluate most common categorical values
#' 
#' @param x Vector of categorical values
#' @return Sorted histogram with counts for each level
fun.heavy <- function(x) {
    histogram <- table(x, useNA='ifany')
    histogram <- sort(-histogram) * -1
    return(histogram)
}

dpHeavyHitters <- setRefClass(
    Class = 'dpHeavyHitters',
    contains = 'mechanismExponential'
)

dpHeavyHitters$methods(
    initialize = function(mechanism, var.type, n, epsilon, k, bins, alpha=0.05, delta=1e-5) {
        .self$name <- 'Differentially private heavy hitters'
        .self$mechanism <- mechanism
        .self$var.type <- check_variable_type(var.type, in_types=c('character', 'factor'))
        .self$n <- n
        .self$epsilon <- epsilon
        .self$k <- k
        .self$bins <- bins
        .self$alpha <- alpha
        .self$delta <- delta
})

dpHeavyHitters$methods(
    release = function(x) {
        .self$result <- export(mechanism)$evaluate(fun.heavy, x, 2, .self$postProcess)
})

dpHeavyHitters$methods(
    postProcess = function(out, gap) {
        out$accuracy <- heavyhitters.getAccuracy(gap, epsilon, delta)
        out$epsilon <- heavyhitters.getParameters(gap, delta, alpha)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
