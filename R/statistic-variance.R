#' Sensitivity of population variance
#' 
#' For a detailed derivation of the sensitivity, see /extra_docs/variance_sensitivity.pdf.
#' @param n Numeric vector of length one; the number of datapoints in the database.
#' @param rng Numeric vector of length two; first entry is minimal bound on the database entries, second is maximal bound on the database entries.
#' @return Numeric vector of length one; a maximal bound on the sensitivity of the population variance.
#'
#' @examples
#' varianceSensitivity(2,c(0,10)) #should return 50
varianceSensitivity <- function(n, rng){
  return(diff(rng)^2/n)
}

#' Postprocessed variance standard deviation
#' 
#' Function to extract standard deviation from noisy estimate of variance.
#'
#' @param release Differentially private release of variance.
#' 
#' @return Noisy estimate of the standard deviation of \code{release}.
#' @rdname postStandDev

postStandDev <- function(release) {
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
    initialize = function(mechanism, varType, variable, n, rng=NULL, epsilon,
                          imputeRng=NULL, alpha=0.05) {
        .self$name <- 'Differentially private variance'
        .self$mechanism <- mechanism
        .self$varType <- varType
        .self$variable <- variable
        .self$n <- checkNValidity(n)
        .self$rng <- checkRange(rng, .self$varType)
        .self$sens <- varianceSensitivity(.self$n, .self$rng)
        
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
        out$std <- postStandDev(out$release)
        return(out)
})
