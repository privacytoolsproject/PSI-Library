#' Function to evaluate true variance & specify post-processing parameters
#'
#' @param x Numeric vector
#' @param var.trpe Character string indicating variable type
#' @param n Integer indicating number of observations in \code{x}
#' @param sensitivity Numeric, the sensitivity of the estimate
#' @param epsilon Numeric, epsilon parameter for differential privacy
#' @return List with fields providing values for the true statistic and post-processing

dp.variance <- function(x, var.type, n, sensitivity, epsilon) {
    out <- list(
        'name' = 'variance',
        'stat' = var(x),
        'var.type' = var.type,
        'n' = n,
        'sensitivity' = sensitivity,
        'epsilon' = epsilon
        )
    return(out)
}


#' Function for differentially private release of variance
#'
#' @param x Numeric vector
#' @param var.trpe Character string indicating variable type
#' @param n Integer indicating number of observations in \code{x}
#' @param epsilon Numeric, epsilon parameter for differential privacy
#' @param rng Numeric 2-tuple indicating range of numeric variable
#' @return List output containing noisy estimate and post-processing values from the Laplace mechanism

variance.release <- function(x, var.type, n, epsilon, rng) {
    rng <- checkrange(rng)
    sensitivity <- (n - 1) / n^2 * diff(rng)^2
    postlist <- c('post.std')
    release <- mechanism.laplace(
        fun=dp.variance,
        x=x,
        var.type=var.type,
        rng=rng,
        sensitivity=sensitivity,
        epsilon=epsilon,
        n=n,
        postlist=postlist
    )
    return(release)
}


#' Function to extract standard deviation from noisy estimate of variance
#'
#' @param release Numeric noisy estimate of variance
#' @return Noisy estimate of standard deviation

post.std <- function(release) {
    std <- sqrt(release)
    return(std)
}
