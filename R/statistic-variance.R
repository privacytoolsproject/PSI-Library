#' Postprocessed variance standard deviation
#' 
#' Function to extract standard deviation from noisy estimate of variance.
#'
#' @param release Differentially private release of variance.
#' 
#' @return Noisy estimate of the standard deviation of \code{release}.
#' @rdname variancePostStandardDeviation

variancePostStandardDeviation <- function(release) {
    std <- sqrt(release)
    return(std)
}


#' Differentially private variance
#'
#' @param mechanism Character, the privacy mechanism.
#' @param varType Character, the R variable type. One of \code{c('numeric',
#'   'integer', 'logical')}.
#' @param Variable Character, variable name.
#' @param n Integer, number of observations
#' @param rng Numeric, a priori estimate of the range
#' @param epsilon Numeric, privacy cost parameter
#' @param accuracy Numeric, accuracy guarantee given \code{epsilon}
#' @param imputeRng Numeric, range within which missing values are imputed. If \code{NULL},
#'   the range provided in \code{rng} is used.
#' @param alpha Numeric, the level of statistical significance. Default 0.05.
#'
#' @import methods
#' @export dpVariance
#' @exportClass dpVariance
#'
#' @include mechanism.R
#' @include mechanism-laplace.R

dpVariance <- setRefClass(
    Class = 'dpVariance',
    contains = c('mechanismLaplace')
)

dpVariance$methods(
    initialize = function(mechanism, varType, variable, n, rng=NULL, epsilon=NULL, accuracy=NULL,
                          imputeRng=NULL, alpha=0.05) {
        .self$name <- 'Differentially private variance'
        .self$mechanism <- mechanism
        .self$varType <- varType
        .self$variable <- variable
        .self$n <- checkN(n)
        .self$rng <- checkRange(rng, varType)
        .self$sens <- (n - 1) / n^2 * diff(.self$rng)^2
        
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- laplaceGetEpsilon(.self$sens, .self$accuracy, alpha)
        } else {
            checkEpsilon(epsilon)
            .self$epsilon <- epsilon
            .self$accuracy <- laplaceGetAccuracy(.self$sens, .self$epsilon, alpha)
        }
        
        if (is.null(imputeRng)) {
            .self$imputeRng <- rng
        } else {
            .self$imputeRng <- checkImputationRange(imputeRng, .self$rng, .self$varType)
        }
        
        .self$alpha <- alpha
})

dpVariance$methods(
    release = function(data) {
        x <- data[, variable]
        .self$result <- export(mechanism)$evaluate(var, x, sens, .self$postProcess)
})

dpVariance$methods(
    postProcess = function(out) {
        out$variable <- variable
        out$std <- sqrt(out$release)
        out$epsilon <- .self$epsilon
        out$accuracy <- .self$accuracy
        return(out)
})
