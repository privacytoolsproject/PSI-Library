#' Generalized Laplace mechanism
#' 
#' Laplace mechanism as applied to a generic object which has a specified
#' add() function that describes how the noise should be added. 

mechanismGenLaplace <- setRefClass(
  Class = 'GenLaplace',
  contains = 'mechanism'
)

mechanismGenLaplace$methods(
  evaluate = function(fun, x, sens, postFun, ...) {
    x <- censordata(x, .self$var.type, .self$rng, .self$bins)
    x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$bins)
    fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
    input.vals <- c(list(x=x), fun.args)
    true.val <- do.call(fun, input.vals)
    scale <- sens / .self$epsilon
    release <- true.val$add(dpNoise(n=length(true.val), scale=scale, dist='laplace'))
    out <- list('release' = release)
    out <- postFun(out, ...)
    return(out)
  }
)