#' Laplace mechanism
#'
#' @import methods
#' @export mechanismLaplace
#'
#' @include mechanism.R

mechanismLaplace <- setRefClass(
    Class = 'mechanismLaplace',
    contains = 'mechanism'
)

mechanismLaplace$methods(
  #' Laplace Mechanism
  #' 
  #' Differentially private evaluation of input function "fun" with sensitivity "sens" on input data 
  #' "x" using the Laplace mechanism.
  #' 
  #' @name LaplaceMechanism
  #' @references C. Dwork, A. Roth The Algorithmic Foundations of Differential Privacy, Chapter 3.3 The Laplace Mechanism p.30-37. August 2014.
  #'
  #' @param fun function of input x to add Laplace noise to.
  #' @param x input that function fun will be evaluated on. 
  #' @param sens sensitivity of fun. Sensitivity is defined in above citation.
  #' @param postFun post-processing function. Takes differentially private release as input
  #'  and returns some form of output in principal based on the differentially private release. 
  #' @param ... any additional (optional) parameters
  #'
  #' @return result of post-processing on input function "fun" evaluated on database "x", assuming sensitivity of fun is "sens".
  #' @export
  #'
  # TODO: add examples 
  evaluate = function(fun, x, sens, postFun, ...) {
    x <- censordata(x, .self$var.type, .self$rng, .self$bins)
    x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$bins)
    fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
    input.vals = c(list(x=x), fun.args)
    true.val <- do.call(fun, input.vals)
    scale <- sens / .self$epsilon
    release <- true.val + dpNoise(n=length(true.val), scale=scale, dist='laplace')
    out <- list('release' = release)
    out <- postFun(out, ...)
    return(out)
  }

)