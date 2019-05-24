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
    evaluate = function(fun, x, sens, postFun, ...) {
        x <- censordata(x, .self$var.type, rng=.self$rng, levels=.self$bins)
        x <- fillMissing(x, .self$var.type, rng=.self$rng, categories=.self$bins)
        fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
        input.vals = c(list(x=x), fun.args)
        true.val <- do.call(fun, input.vals)  # Concern: are we confident that the environment this is happening in is getting erased.
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
