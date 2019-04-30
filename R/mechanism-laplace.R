#' Laplace mechanism
#'
#' @import methods
#' @export mechanismLaplace
#' @exportClass mechanismLaplace
#'
#' @include mechanism.R

mechanismLaplace <- setRefClass(
    Class = 'mechanismLaplace',
    contains = 'mechanism'
)

mechanismLaplace$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismLaplace$methods(
    evaluate = function(fun, x, sens, postFun, ...) {
        x <- censordata(x, .self$var.type, .self$rng, .self$bins)
        x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$bins)
        field.vals <- .self$getFunArgs(fun)
        ellipsis.vals <- getFuncArgs(list(...), fun)
        true.val <- do.call(fun, c(list(x=x), field.vals, ellipsis.vals))
        scale <- sens / .self$epsilon
        release <- true.val + dpNoise(n=length(true.val), scale=scale, dist='laplace')
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})
