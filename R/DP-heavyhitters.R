#' Function to evaluate most common values and specify arguments to post-processing
#'
#' @param x Vector of categorical data (character, factor)
#' @param var.type Character string indicating data type
#' @param epsilon Epsilon value for differential privacy
#' @param n Integer indicating number of observations
#' @param k Integer querying the most common \code{k} categories
#' @param sensitivity Numeric the sensitivity of the statistic
#' @return List with true value of statistic and parameters to be passed to post-processing

dp.heavyhitters <- function(x, var.type, epsilon, n, k, sensitivity) {
    hist <- table(x, useNA='ifany')
    if (k > length(hist) - 1) { stop('failure: k too large') }
    idx <- 1:(k+1)
    hist <- sort(-hist, partial=idx)[idx] * -1
    gap <- as.numeric(hist[k] - hist[k+1])
    out <- list('name' = 'heavyhitters',
                'stat' = gap,
                'var.type' = var.type,
                'k' = k,
                'epsilon' = epsilon,
                'heavyhitters' = names(hist)[1:k])
    return(out)
}


#' Release differentially private most common values
#'
#' @param x Vector of categorical data (character, factor)
#' @param var.type Character string indicating data type
#' @param epsilon Epsilon value for differential privacy
#' @param n Integer indicating number of observations
#' @param k Integer querying the most common \code{k} categories
#' @param bins Vector of categories from which the top \code{k} categories are evaluated
#' @return A vector whose values are \code{k} most common categories, or a failure message
#' @author Victor Balcer
#'
#' Contains functions for privately releasing the top-k heavy hitters, that is,
#'  the k most common unique values.  For example, k=1 releases the mode.
#'
#' @examples
#' N <- 1000
#' epsilon <- 0.5
#' observed.levels <- bins <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
#' probs <- c(0.40, 0.25, 0.15, 0.10, 0.04, 0.03, 0.02, 0.01)
#' y <- sample(observed.levels, size=N, prob=probs, replace=TRUE)
#' release <- heavyhitters.release(y, 'character', epsilon, N, k=3, bins=bins)

heavyhitters.release <- function(x, var.type, epsilon, n, k, bins) {
    var.type <- check_variable_type(var.type, in_types=c('character', 'factor'))
    postlist <- list('accuracy' = 'getAccuracy',
                     'epsilon' = 'getParameters')
    release <- mechanism.laplace(fun=dp.heavyhitters, x=x, var.type=var.type, rng=NULL,
                                 epsilon=epsilon, sensitivity=1, n=n, k=k, bins=bins,
                                 postlist=postlist)
    if (release$release < -2 / epsilon * log(1e-7)) {
        release$release <- 'failure: gap too small'  # or throw error here?
    } else {
        release$release <- release$heavyhitters
    }
    release <- release[names(release) != 'heavyhitters']  # remove redundant element
    return(release)
}


#' Get the accuracy of heavyhitters statistic for a given value of epsilon
#'
#' @param release The noisy estimate
#' @param epsilon Epsilon value for differential privacy
#' @param delta Delta value for differential privacy
#' @return The accuracy guaranteed by the given epsilon
#' @author Victor Balcer
#'
#' The accuracy is interpreted as follows: The alpha value returned means that with
#'   probability 1 - beta, simultaneously for all t with min <= t <= max, the algorithm's
#'   estimate of the count in [min, t] is within alpha of the true value.

heavyhitters.getAccuracy <- function(release, epsilon, delta=1e-7) {
    accuracy <- exp(-epsilon * release / 2) / delta
    return(accuracy)
}


#' Get the epsilon value necessary to guarantee a desired level of accuracy of a heavyhitters release
#'
#' @param release The noisy estimate
#' @param delta Delta value for differential privacy
#' @param alpha The true value is within the accuracy range with probability 1 - \code{alpha}
#' @return The epsilon value necessary to gaurantee the given accuracy
#' @author Victor Balcer

heavyhitters.getParameters <- function(release, delta=1e-7, alpha=0.05) {
  epsilon <- -2 * log(alpha * delta) / release
  return(epsilon)
}
