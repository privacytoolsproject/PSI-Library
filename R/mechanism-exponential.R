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
        x <- censordata(x, .self$var.type, levels=.self$bins)
        x <- fillMissing(x, .self$var.type, categories=.self$bins)   # Problem that this needs to work over other types than categorical
        field.vals <- .self$getFunArgs(fun)
        ellipsis.vals <- getFuncArgs(list(...), fun)
        true.val <- do.call(fun, c(list(x=x), field.vals, ellipsis.vals))
        quality <- true.val - max(true.val)
        likelihoods <- exp((.self$epsilon * quality) / (2 * sens))   # Problem that this zeroes out for quality differences > 800
        probs <- likelihoods/sum(likelihoods)
        release <- sample(names(true.val), size=.self$k, prob=probs) # Problem that this needs to use openssl randomness
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})
