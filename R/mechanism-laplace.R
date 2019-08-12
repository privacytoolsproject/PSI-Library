#' Laplace mechanism
#'
#' @import methods
#' @export mechanismLaplace
#' @exportClass mechanismLaplace
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
  #' @name Laplace Mechanism
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
  #'
  #' @examples
  #' # histogram example
  #' 
  #' # the function in `statistic-histogram.R` that creates a histogram from input data
  #' histogram_function <- fun.hist 
  #' # the data frame that holds the data we want to analyze, in this case the data is called "PUMS5extract10000"
  #' data <- data(PUMS5extract10000) 
  #' # the variable for which we want a histogram
  #' variable <- "age"
  #' # the sensitivity for the histogram, the default sensitivity for histograms is 2 
  #' sens <- 2 
  #' # the post-processing function to use to format the histogram release correctly
  #' post_processing_function <- dpHistogram$postProcess 
  #' 
  #' laplace_histogram <- mechanismLaplace$evaluate(histogram_function, data[, variable], sens, post_processing_function)
  #' 
  #' # mean example
  #' 
  #' mean_function <- mean
  #' # the sensitivity for a differntially private mean in calculated as the difference in the data range divided by the number of data points
  #' sens <- diff(rng) / n 
  #' # the post-processing function to use to format the mean release correctly
  #' post_processing_function <- dpMean$postProcess 
  #' # `data` and `variable` same as above
  #' 
  #' laplace_mean <- mechanismLaplace$evaluate(mean_function, data[, variable], sens, post_processing_function)
  #' 
  evaluate = function(fun, x, sens, postFun, ...) {
    x <- censordata(x, .self$var.type, .self$rng, .self$bins)
    x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$impute.bins)
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