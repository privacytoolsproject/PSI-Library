#' Exponential mechanism
#'
#' @import methods
#' @export mechanismExponential
#' @exportClass mechanismExponential
#'
#' @include mechanism.R

mechanismExponential <- setRefClass(
    Class = 'mechanismExponential',
    contains = 'mechanism'
)

mechanismExponential$methods(
    getFunArgs = function(fun) { 
        callSuper(fun)
})

mechanismExponential$methods(
    evaluate = function(fun, x, sens, postFun, ...) {
        x <- censordata(x, .self$var.type, rng=.self$rng, levels=.self$bins)
        x <- fillMissing(x, .self$var.type, rng=.self$rng, categories=.self$bins)
        field.vals <- .self$getFunArgs(fun)
        true.val <- do.call(fun, c(list(x=x), field.vals))
        quality <- true.val - max(true.val)
        probs <- ifelse(true.val == 0, 0, exp((.self$epsilon * quality) / (2 * sens)))
        gap <- as.numeric(true.val[.self$k] - true.val[.self$k + 1])
        if (gap < (-2 / epsilon * log(delta))) {
            out <- list('release' = NULL)
        } else {
            release <- sample(names(true.val), size=.self$k, prob=probs)
            out <- list('release' = release)
            out <- postFun(out, gap)
        }
        return(out)
})
