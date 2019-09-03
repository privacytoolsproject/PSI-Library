#' Differentially private mean
#'
#' @param mechanism Character, the privacy mechanism. For \code{dpMean}, one
#'   of \code{c('mechanismLaplace', 'mechanismBootstrap')}.
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
#' @param nBoot Integer, the number of bootstrap replications if using the bootstrap
#'   bootstrap mechanism, ignored otherwise. Default 20.
#'
#' @import methods
#' @export dpMean
#' @exportClass dpMean
#'
#' @include mechanism.R
#' @include mechanism-laplace.R
#' @include mechanism-bootstrap.R

dpMean <- setRefClass(
    Class = 'dpMean',
    contains = c('mechanismLaplace', 'mechanismBootstrap')
)

dpMean$methods(
    initialize = function(mechanism, varType, variable, n, rng=NULL, epsilon=NULL,
                          accuracy=NULL, imputeRng=NULL, alpha=0.05, nBoot=20, ...) {
        .self$name <- 'Differentially private mean'
        .self$mechanism <- mechanism
        .self$varType <- varType
        .self$variable <- variable
        .self$n <- checkN(n)
        .self$alpha <- alpha
        .self$rng <- checkRange(rng, varType)
        .self$sens <- diff(.self$rng) / n
        
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- laplaceGetEpsilon(.self$sens, .self$accuracy, alpha)
        } else {
            checkEpsilon(epsilon)
            .self$epsilon <- epsilon
            .self$accuracy <- laplaceGetAccuracy(.self$sens, .self$epsilon, alpha)
        }
        
        if (is.null(imputeRng)) {
            .self$imputeRng <- .self$rng
        } else {
            .self$imputeRng <- checkImputationRange(imputationRange=imputeRng, rng=.self$rng, varType=.self$varType)
        }
        
        .self$bootFun <- bootMean
        .self$nBoot <- nBoot
})


dpMean$methods(
    release = function(data, ...) {
        x <- data[, variable]
        .self$result <- export(mechanism)$evaluate(mean, x, .self$sens, .self$postProcess, ...)
})

dpMean$methods(
    postProcess = function(out) {
        out$variable <- variable
        if (mechanism == 'mechanismLaplace') {
            out$accuracy <- accuracy
            out$epsilon <- epsilon
            out$interval <- meanGetCI(out$release, epsilon, .self$sens, alpha)
        } 
        if (varType == 'logical') {
            if (mechanism == 'mechanismBootstrap') {
                baggedEstimate <- mean(out$release)
                out$stdDev <- meanPostStandardDeviation(baggedEstimate)
                out$median <- meanPostMedian(baggedEstimate)
                out$histogram <- meanPostHistgram(baggedEstimate)
            } else {
                out$stdDev <- meanPostStandardDeviation(out$release)
                out$median <- meanPostMedian(out$release)
                out$histogram <- meanPostHistogram(out$release, n)
            }
        }
        return(out)
})
