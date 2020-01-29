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
    evaluate = function(fun, x, sens, ...) {
        x <- censorData(x, .self$varType, rng=.self$rng, levels=.self$bins)
        x <- fillMissing(x, .self$varType, rng=.self$rng, categories=.self$bins)
        fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
        inputVals = c(list(x=x), fun.args)
        trueVal <- do.call(fun, inputVals)
        quality <- trueVal - max(trueVal)
        probs <- ifelse(trueVal == 0, 0, exp((.self$epsilon * quality) / (2 * sens)))
        gap <- as.numeric(trueVal[.self$k] - trueVal[.self$k + 1])
        if (gap < (-2 / .self$epsilon * log(.self$delta))) { # THIS IS NOT DIFFERENTIALLY PRIVATE SINCE THE GAP IS NOT DP!
            out <- list('release' = NULL)
        } else {
            release <- sample(names(trueVal), size=.self$k, prob=probs)
            out <- list('release' = release)
            
        }
        return(out)
})
