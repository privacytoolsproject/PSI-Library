#' Snapping mechanism
#'
#' @import methods
#' @export mechanismSnapping
#' @exportClass mechanismSnapping
#'
#' @include mechanism.R

mechanismSnapping <- setRefClass(
    Class = 'mechanismSnapping',
    contains = 'mechanism'
)

mechanismSnapping$methods(
  #' Snapping Mechanism
  #'
  #' Differentially private evaluation of input function "fun" with sensitivity "sens" on input data
  #' "x" using the Snapping mechanism.
  #'
  #' @name Snapping Mechanism
  #' @references Ilya Mironov,  On Significance of the Least Significant Bits for Differential Privacy. 2012
  #'
  #' @param fun function of input x for which we will add noise to the output.
  #' @param x input that function fun will be evaluated on.
  #' @param sens sensitivity of fun. Sensitivity is defined in above citation.
  #' @param B bound such that the user is confident that fun(x) is in the interval [-B, B].
  #' @param postFun post-processing function. Takes differentially private release as input
  #'  and returns some form of output in principal based on the differentially private release.
  #' @param ... any additional (optional) parameters
  #'
  #' @return result of post-processing on input function "fun" evaluated on database "x", assuming sensitivity of fun is "sens".
  #'
  evaluate = function(fun, x, sens, B, postFun, ...) {
    x <- censordata(x, .self$var.type, .self$rng, .self$bins)
    x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$impute.bins)
    fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
    input.vals = c(list(x=x), fun.args)
    true.val <- do.call(fun, input.vals)  # Concern: are we confident that the environment this is happening in is getting erased.
    scale <- sens / .self$epsilon
    release <- true.val + snappingNoise(true.val, n = length(true.val), sens, .self$epsilon, B)
    out <- list('release' = release)
    out <- postFun(out, ...)
    return(out)
  }

)
