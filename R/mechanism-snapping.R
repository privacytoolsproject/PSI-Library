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
  #' @param postFun post-processing function. Takes differentially private release as input
  #'  and returns some form of output in principal based on the differentially private release.
  #' @param ... any additional (optional) parameters
  #'
  #' @return result of post-processing on input function "fun" evaluated on database "x", assuming sensitivity of fun is "sens".
  #'
  evaluate = function(fun, x, sens, postFun, ...) {
    x <- censordata(x, .self$var.type, .self$rng, .self$bins)
    x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$impute.bins)
    fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
    input.vals = c(list(x=x), fun.args)
    true.val <- do.call(fun, input.vals)  # Concern: are we confident that the environment this is happening in is getting erased.
    scale <- sens / .self$epsilon

    # TODO: think more about vectors of noise -- might need to loop over true.val and set
    #       B separately for each

    # TODO: Below is attempt to appropriately set B for each type of function for which the snapping mechanism might be used.
    #       Need to check further with someone (Ira maybe?) about how best to do this.
    if (identical(fun, mean)) {
        B <- max(abs(.self$rng[1]), abs(.self$rng[2]))
    } else if (identical(fun, var)) {
        B <- (.self$rng[2] - .self$rng[1])^2 / 4
    } else if (identical(fun, fun.covar)) {
        # TODO: figure this out
    }
    else if (identical(fun, fun.hist) || identical(fun, boot.hist)) {
        B <- length(x)
    } else {
        stop(paste0('The Snapping Mechanism does not support the ', fun, ' function.'))
    }
    n = length(true.val)
    release <- true.val + snappingNoise(true.val, n, sens, .self$epsilon, B)
    out <- list('release' = release)
    out <- postFun(out, ...)
    return(out)
  }

)
