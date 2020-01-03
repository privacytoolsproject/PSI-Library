#' Gaussian mechanism
#'
#' @import methods
#' @export mechanismGaussian
#' @exportClass mechanismGaussian
#'
#' @include mechanism.R

mechanismGaussian <- setRefClass(
    Class = 'mechanismGaussian',
    contains = 'mechanism'
)

mechanismGaussian$methods(
    evaluate = function(fun, x, sens, ...) {
        x <- censorData(x, .self$varType, .self$rng, rngFormat=.self$rngFormat)
        x <- fillMissing(x, .self$varType, imputeRng=.self$rng, categories=.self$bins)
        funArgs <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
        inputVals = c(list(x=x), funArgs)
        trueVal <- do.call(fun, inputVals)
        scale <- sens * sqrt(2 * log(1.25 / .self$delta)) / .self$epsilon
        release <- trueVal + dpNoise(n=length(trueVal), scale=scale, dist='gaussian')
        out <- list('release' = release)
        return(out)
})
