#' Snapping mechanism
#'
#' @import methods
#' @export mechanismSnapping
#' @exportClass mechanismSnapping
#'
#' @include mechanism.R

mechanismSnapping <- setRefClass(
    Class = 'mechanismSnapping',
    contains = 'mechanism',
    fields = list(min_B = 'numeric')
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
  #' @param postFun post-processing function. Takes differentially private release as input
  #'  and returns some form of output in principal based on the differentially private release.
  #' @param ... any additional (optional) parameters
  #'
  #' @return result of post-processing on input function "fun" evaluated on database "x", assuming sensitivity of fun is "sens".
  #'
  evaluate = function(fun, x, sens, postFun, ...) {
    x <- censorData(x, .self$varType, .self$rng, .self$bins, .self$rngFormat)
    x <- fillMissing(x, .self$varType, imputeRng=.self$rng, categories=.self$imputeBins)
    fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
    inputVals = c(list(x=x), fun.args)
    trueVal <- do.call(fun, inputVals)  # Concern: are we confident that the environment this is happening in is getting erased.
    scale <- sens / .self$epsilon
    release <- trueVal + snappingNoise(trueVal, n=length(trueVal), sens, .self$epsilon, .self$min_B)
    out <- list('release' = release)
    out <- postFun(out, ...)
    return(out)
  }
)
