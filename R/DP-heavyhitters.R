#' Function to evaluate most common values and specify arguments to post-processing
#'
#' @param x Vector of categorical data (character, factor)
#' @param var.type Character string indicating data type
#' @param epsilon Epsilon value for differential privacy
#' @param delta something here
#' @param n Integer indicating number of observations
#' @param k Integer querying the most common \code{k} categories
#' @param sensitivity Numeric the sensitivity of the statistic
#' @return List with true value of statistic and parameters to be passed to post-processing

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


#' Release differentially private most common values
#'
#' @param x Vector of categorical data (character, factor)
#' @param var.type Character string indicating data type
#' @param epsilon Epsilon value for differential privacy
#' @param delta something here
#' @param n Integer indicating number of observations
#' @param k Integer querying the most common \code{k} categories
#' @param bins Vector of categories from which the top \code{k} categories are evaluated
#' @return A vector whose values are \code{k} most common categories, or a failure message
#'
#' Contains functions for privately releasing the top-k heavy hitters, that is,
#'  the k most common unique values.  For example, k=1 releases the mode.
#'
#' @examples
#' N <- 1000
#' epsilon <- 0.5
#' delta <- 1e-7
#' observed.levels <- bins <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
#' probs <- c(0.40, 0.25, 0.15, 0.10, 0.04, 0.03, 0.02, 0.01)
#' y <- sample(observed.levels, size=N, prob=probs, replace=TRUE)
#' release <- heavyhitters.release(y, 'character', epsilon, delta, N, k=3, bins=bins)

heavyhitters.release <- function(x, var.type, epsilon, delta, n, k, bins) {
    var.type <- check_variable_type(var.type, in_types=c('character', 'factor'))
    postlist <- list('accuracy' = 'getAccuracy',
                     'epsilon' = 'getParameters',
                     'failed' = 'postNoteFailure')
    release <- mechanism.exponential(fun=dp.heavyhitters, x=x, var.type=var.type,
                                     epsilon=epsilon, delta=delta, sensitivity=2, n=n, k=k,
                                     bins=bins, postlist=postlist)
    return(release)
}


#' Get the accuracy of heavyhitters statistic for a given value of epsilon
#'
#' @param gap Gap
#' @param epsilon Epsilon value for differential privacy
#' @param delta Delta value for differential privacy
#' @return The accuracy guaranteed by the given epsilon
#' @author Victor Balcer
#'
#' The accuracy is interpreted as follows: The alpha value returned means that with
#'   probability 1 - beta, simultaneously for all t with min <= t <= max, the algorithm's
#'   estimate of the count in [min, t] is within alpha of the true value.

heavyhitters.getAccuracy <- function(gap, epsilon, delta) {
    accuracy <- exp(-epsilon * gap / 2) / delta
    return(accuracy)
}


#' Get the epsilon value necessary to guarantee a desired level of accuracy of a heavyhitters release
#'
#' @param gap Gap
#' @param delta Delta value for differential privacy
#' @param alpha The true value is within the accuracy range with probability 1 - \code{alpha}
#' @return The epsilon value necessary to gaurantee the given accuracy
#' @author Victor Balcer

heavyhitters.getParameters <- function(gap, delta, alpha=0.05) {
  epsilon <- -2 * log(alpha * delta) / gap
  return(epsilon)
}


#' Function to see failed
#' @param failed something here
#' @return something here

heavyhitters.postNoteFailure <- function(failed) {
    return(failed)
}
