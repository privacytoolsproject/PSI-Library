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
    postlist <- list('std' = 'post.std')
    release <- mechanism.laplace(fun=dp.variance, x=x, var.type=var.type, rng=rng,
                                 sensitivity=sensitivity, epsilon=epsilon, n=n, postlist=postlist)
    return(release)
}


#' Function to extract standard deviation from noisy estimate of variance
#'
#' @param release Numeric noisy estimate of variance
#' @return Noisy estimate of standard deviation

variance.post.std <- function(release) {
    std <- sqrt(release)
    return(std)
}


# --------------------------------------------------------- #
# --------------------------------------------------------- #
# Reference class implementation of variance

dpVariance <- setRefClass(
    Class = 'dpVariance',
    contains = c('mechanismLaplace')
)

dpVariance$methods(
    initialize = function(mechanism, var.type, n, epsilon, rng, alpha=0.05) {
        .self$name <- 'Differentially private variance'
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$n <- n
        .self$epsilon <- epsilon
        .self$rng <- rng
        .self$alpha <- alpha
})

dpVariance$methods(
    release = function(x) {
        sens <- (n - 1) / n^2 * diff(rng)^2
        .self$result <- export(mechanism)$evaluate(var, x, sens, .self$postProcess)
})

dpVariance$methods(
    postProcess = function(out) {
        out$std <- sqrt(out$release)
        return(out)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
