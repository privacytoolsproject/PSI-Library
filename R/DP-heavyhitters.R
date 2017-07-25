#' DP Heavyhitters
#' 
#' Function to evaluate most common values and specify arguments to post-processing
#'
#' @param x A vector of categorical data (character, factor).
#' @param var.type A character vector specifying variable type of \code{x}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param delta A numeric vector representing the probability of an arbitrary
#'    leakage of information from \code{x}. Should be of length one 
#'    and should be a very small value.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param k An integer querying the most common \code{k} categories.
#' @param sensitivity The difference of the range of the \code{x} divided 
#'    by \code{n}.
#'    
#' @return List with true value of statistic and parameters to be passed to post-processing
#' @rdname dp.heavyhitters
dp.heavyhitters <- function(x, var.type, epsilon, delta, n, k, sensitivity) {
    hist <- table(x, useNA='ifany')
    if (k > length(hist) - 1) { stop('failure: k too large') }
    hist <- sort(-hist) * -1
    gap <- as.numeric(hist[k] - hist[k+1])
    failed <- gap < -2 / epsilon * log(delta)
    out <- list('name' = 'heavyhitters',
                'stat' = hist,
                'var.type' = var.type,
                'k' = k,
                'epsilon' = epsilon,
                'delta' = delta,
                'gap' = gap,
                'failed' = failed)
    return(out)
}


#' Release differentially private heavyhitters
#' 
#' A function for privately releasing the top-k heavy hitters, that is,
#'  the k most common unique values.  For example, k=1 releases the mode.
#'
#' @param x A vector of categorical data (character, factor).
#' @param var.type A character vector specifying variable type of \code{x}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param k An integer querying the most common \code{k} categories.
#' @param bins A vector of categories from which the top \code{k} categories 
#'    are evaluated.
#' @param delta A numeric vector representing the probability of an arbitrary
#'    leakage of information from \code{x}. Should be of length one 
#'    and should be a very small value. Default to 10^-6.
#'    
#' @return A vector whose values are \code{k} most common categories, or a 
#'    failure message.
#' @examples
#' 
#' N <- 1000
#' epsilon <- 0.5
#' delta <- 1e-7
#' observed.levels <- bins <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
#' probs <- c(0.40, 0.25, 0.15, 0.10, 0.04, 0.03, 0.02, 0.01)
#' y <- sample(observed.levels, size=N, prob=probs, replace=TRUE)
#' release <- heavyhitters.release(y, 'character', epsilon, delta, N, k=3, bins=bins)
#' @rdname heavyhitters.release
#' @export
heavyhitters.release <- function(x, var.type, epsilon, n, k, bins, delta=0.000001) {
    var.type <- check_variable_type(var.type, in_types=c('character', 'factor'))
    postlist <- list('accuracy' = 'getAccuracy',
                     'epsilon' = 'getParameters',
                     'failed' = 'postNoteFailure')
    release <- mechanism.exponential(fun=dp.heavyhitters, x=x, var.type=var.type,
                                     epsilon=epsilon, delta=delta, sensitivity=2, n=n, k=k,
                                     bins=bins, postlist=postlist)
    return(release)
}


#' Heavyhitters accuracy
#' 
#' Get the accuracy of heavyhitters statistic for a given value of epsilon.
#' The accuracy is interpreted as follows: The alpha value returned means that with
#' probability 1 - beta, simultaneously for all t with min <= t <= max, the algorithm's
#' estimate of the count in [min, t] is within alpha of the true value.
#'
#' @param gap 
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
#' @param gap 
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


#' Heavyhitters failure
#' 
#' Function to see if the heavyhitters mechanism failed.
#' 
#' @param failed A logical vector.
#' 
#' @return \code{failed}. If \code{TRUE}, the mechanism failed. If \code{FALSE},
#'    the function succeeded.
#' @rdname heavyhitters.postNoteFailure
heavyhitters.postNoteFailure <- function(failed) {
    return(failed)
}

# --------------------------------------------------------- #
# --------------------------------------------------------- #
# dp heavyhitters class

fun.heavy <- function(x, var.type) {
    hist <- table(x, useNA='ifany')
    hist <- sort(-hist) * -1
    return(hist)
}

dpHeavyHitters <- setRefClass(
    Class = 'dpHeavyHitters',
    contains = 'mechanismExponential'
)

dpHeavyHitters$methods(
    initialize = function(mechanism, var.type, n, epsilon, k, bins, alpha=0.05, delta=1e-7) {
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
        .self$result <- export(mechanism)$evaluate(fun.hist, x, 2, .self$postProcess)
})

dpHeavyHitters$methods(
    postProcess = function(out) {
        gap <- as.numeric(out$release[k] - out$release[k + 1])
        if (gap < -2 / epsilon * log(delta)) {
            out$release <- NULL
            return(out)
        }
        out$accuracy <- heavyhitters.getAccuracy(gap, epsilon, delta)
        out$epsilon <- heavyhitters.getParameters(gap, delta, alpha)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
