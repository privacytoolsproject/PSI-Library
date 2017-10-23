#' DP Covariance
#' 
#' Function to evaluate true variance & specify post-processing parameters.
#'
#' @param x A numeric, logical, or integer vector.
#' @param var.type A character vector specifying variable type of \code{x}. 
#'    Should be of length one and should contain either 'numeric', 
#'    'logical', or 'integer'.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param sensitivity The difference of the range of \code{x} divided 
#'    by \code{n}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#'    
#' @return A list with fields providing values for the true statistic and post-processing.
#' @rdname dp.variance
#' @export
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


#' Release differentially private variance
#' 
#' Function to calculate a differentially private release of variance.
#'
#' @param x A numeric, logical, or integer vector.
#' @param var.type A character vector specifying variable type of \code{x}. 
#'    Should be of length one and should contain either 'numeric', 
#'    'logical', or 'integer'.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param rng A numeric vector specifying an a priori estimate of the range
#'    of \code{x}. Should be of length two.
#' @param impute.rng Numeric range within which missing values of \code{x} are 
#'    imputed. Defaults to \code{rng} if \code{NULL}
#'    
#' @return A list output containing noisy estimate and post-processing values from the 
#'    Laplace mechanism.
#' @examples
#' 
#' n <- 1000
#' x_num <- runif(n)
#' x_bool <- x_num >= 0.5
#' r_num <- variance.release(x=x_num, var.type = 'numeric', n=n, epsilon = 0.1, rng = c(0,1))
#' r_bool <- variance.release(x=x_bool, var.type='logical', epsilon=0.5, n=n, rng=c(0, 1))
#' @rdname variance.release
#' @export
variance.release <- function(x, var.type, n, epsilon, rng, impute.rng=NULL) {
    rng <- checkrange(rng)
    if (is.null(impute.rng)) { impute.rng <- rng }
    sensitivity <- (n - 1) / n^2 * diff(rng)^2
    postlist <- list('std' = 'postStandardDeviation')
    release <- mechanism.laplace(fun=dp.variance, x=x, var.type=var.type, rng=rng, impute.rng=impute.rng,
                                 sensitivity=sensitivity, epsilon=epsilon, n=n, postlist=postlist)
    return(release)
}


#' Postprocessed variance standard deviation
#' 
#' Function to extract standard deviation from noisy estimate of variance.
#'
#' @param release Differentially private release of variance.
#' 
#' @return Noisy estimate of the standard deviation of \code{release}.
#' @rdname variance.postStandardDeviation
variance.postStandardDeviation <- function(release) {
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
