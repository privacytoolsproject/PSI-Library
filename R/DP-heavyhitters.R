# might require having the mechanism use and amend the list returned by the dp.x functions, 
# instead of assigning a new list that is a subset. this would allow the 'heavyhitters' element
# to be viewed by the post-processing calls in the release fn. probably would not need a separate
# and redundant postprocessing function then. 
#
# just need to be sure that attributes in the release aren't duplicated (e.g., epsilon)

#' Function to evaluate most common values and specify arguments to post-processing
#'
#' @param x

dp.heavyhitters <- function(x, var.type, n, epsilon, sensitivity, k) {
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
#' @param x A vector of the data
#' @param epsilon Epsilon value for differential privacy
#' @param delta Delta value for differential privacy
#' @param k top-k elements to return
#' @param symb Symbol to return on failure (default : NA)
#' @return A vector whose values are k most common values or a failure return symbol
#' @author Victor Balcer
#'
#' Contains functions for privately releasing the top-k heavy hitters, that is,
#'  the k most common unique values.  For example, k=1 releases the mode.
#'
#' @examples
#' n <- 10000
#' range <- c(0,20)
#' x <- rbinom(n, size=max(range), prob=0.7)
#' heavyhitters.release(x=x, epsilon=.1, delta=.0000001)

heavyhitters.release <- function(x, var.type, epsilon, n, k, bins) {
    var.type <- check_variable_type(var.type, in_types=c('character', 'factor'))
    postlist <- list('accuracy' = 'getAccuracy',
                     'epsilon' = 'getParameters')
    release <- mechanism.laplace(fun=dp.heavyhitters, x=x, var.type=var.type, rng=NULL,
                                 epsilon=epsilon, sensitivity=1, n=n, k=k, bins=bins,
                                 postlist=postlist)
    if (release$release < -2 / epsilon * log(1e-7)) {
        release$release <- 'failure: gap too small'
    } else {
        release[['release']] <- release$heavyhitters
    }
    return(release)
}


#heavyhitters.release <- function(x, epsilon, delta, k=1, symb=NA){
#
#  hist <- table(factor(x, exclude=c()), useNA="ifany")
#
#  # TODO: possible violation of privacy
#  if(k > length(hist) - 1)
#  {
#    return(symb)
#  }
#
#  # TODO: should be replaced with a partial sorting algorithm
#  hist <- sort(hist, decreasing=TRUE)[1:(k+1)]
#
#  gap <- hist[k] - hist[k+1] + rlaplace(n=1, sensitivity=2, epsilon=epsilon)
#
#  testfailure <- gap < -2 / epsilon * log(delta)
#
#  if(testfailure){
#    return(symb)
#  }else{
#    return(names(hist[1:k]))
#  }
#
#}


#' Get the accuracy of heavyhitters statistic for a given value of epsilon
#'
#' @param epsilon Epsilon value for differential privacy
#' @param delta Delta value for differential privacy
#' @param beta The true value is within the accuracy range with probability 1-beta
#' @param gap
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
#' @param delta Delta value for DP
#' @param beta the true value is within the accuracy range (alpha) with probability 1-beta
#' @param gap something here
#' @return The epsilon value necessary to gaurantee the given accuracy
#' @author Victor Balcer

heavyhitters.getParameters <- function(release, delta=1e-7, alpha=0.05) {
  epsilon <- -2 * log(alpha * delta) / release
  return(epsilon)
}
