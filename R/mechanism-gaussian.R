#' Gaussian mechanism
#'
#' @import methods
#' @export mechanismGaussian
#' @exportClass mechanismGaussian

mechanismGaussian <- setRefClass(
    Class = 'mechanismGaussian',
    contains = 'mechanism'
)

mechanismGaussian$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismGaussian$methods(
    evaluate = function(fun, x, sens, postFun) {
        x <- censordata(x, .self$var.type, .self$rng)
        if (.self$var.type %in% c('numeric', 'integer', 'logical')) {
            if (NCOL(x) > 1) {
                x <- fillMissing2d(x, .self$var.type, .self$impute.rng)
            } else {
                x <- fillMissing(x, .self$var.type, .self$impute.rng[1], .self$impute.rng[2])
            }
        } else {
            x <- fillMissing(x, .self$var.type, categories=.self$bins)
        }
        field.vals <- .self$getFunArgs(fun)
        true.val <- do.call(fun, c(list(x=x), field.vals))
        scale <- sens * sqrt(2 * log(1.25 / .self$delta)) / .self$epsilon
        release <- true.val + dpNoise(n=length(true.val), scale=scale, dist='gaussian')
        out <- list('release' = release)
        out <- postFun(out)
        return(out)
})
