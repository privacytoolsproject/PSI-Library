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
    evaluate = function(fun, x, sens, postFun, ...) {
        x <- censorData(x, .self$var.type, .self$rng)
        x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$bins)
        fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
        input.vals = c(list(x=x), fun.args)
        true.val <- do.call(fun, input.vals)
        scale <- sens * sqrt(2 * log(1.25 / .self$delta)) / .self$epsilon
        release <- true.val + dpNoise(n=length(true.val), scale=scale, dist='gaussian')
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})
