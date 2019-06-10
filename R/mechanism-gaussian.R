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
    evaluate = function(fun, x, sens, postFun, stability, ...) {
        x <- censordata(x, .self$var.type, .self$rng)
        x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$bins)
        fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
        input.vals = c(list(x=x), fun.args)
        true.val <- do.call(fun, input.vals)
        scale <- sens * sqrt(2 * log(1.25 / .self$delta)) / .self$epsilon
        # release <- true.val + dpNoise(n=length(true.val), scale=scale, dist='gaussian')
        noiseVector <- dpNoise(n=length(true.val), scale=scale, dist='gaussian')
        if (stability) {
            # if this is for a stability mechanism, only add noise to bins
            # that are NOT equal to zero
            release <- ifelse(true.val == 0, true.val, true.val + noiseVector)
        } else {
            # if this is NOT stability mechanism, add noise to all bins
            release <- true.val + noiseVector
        }
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})
