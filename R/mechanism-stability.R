#' Stability mechanism
#'
#' @import methods
#' @export mechanismStability
#' @exportClass mechanismStability
#'
#' @include mechanism.R

mechanismStability <- setRefClass(
    Class = 'mechanismStability',
    contains = 'mechanism'
)

mechanismStability$methods(
    #' Stability Mechanism
    #' 
    #' Differentially private evaluation of input function "fun" with sensitivity "sens" on input data 
    #' "x" using the Laplace mechanism.
    #' 
    #' @name Stability Mechanism
    #' @references C. Dwork, A. Roth The Algorithmic Foundations of Differential Privacy, Chapter 7.3 Stabiilty and Privacy p.150-157. August 2014.
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
    # TODO: add examples 
    evaluate = function(fun, x, sens, postFun, ...) {
        # before calculating the histogram statistic, confirm that delta is less than (1/n)
        # if deta is greater than or equal to (1/n), return an error message to the user
        if (.self$delta >= (1 / .self$n)) stop("Delta must be less than 1/n")
        
        # if the variable is numeric or integer and the stability mechanism is being used,
        # then the stability mechanism needs to determine the bins to maintain privacy.
        # Get the range of the data, then get the number of bins from the input number of
        # bins or the input granularity
        dataRange <- NULL
        numHistogramBins <- NULL
        imputationRange <- NULL
        histogramBins <- NULL
        if (.self$var.type %in% c('numeric', 'integer')) {
            dataRange <- range(x)
            numHistogramBins <- ifelse(is.null(.self$n.bins), .self$n / .self$granularity, .self$n.bins)
            histogramBins <- seq(dataRange[1], dataRange[2], length.out=(numHistogramBins + 1))
            # set the imputation range to the detected data range to maintain privacy,
            # a user could have entered an imputation range without entering a range
            imputationRange <- dataRange
        }
        
        x <- censordata(x, .self$var.type, dataRange, histogramBins)
        x <- fillMissing(x, .self$var.type, impute.rng=imputationRange, categories=levels(x)) # levels(x) will be NULL for numeric variables, a vector of bins for character variables
        fun.args <- getFuncArgs(fun, inputList=list(bins=histogramBins), inputObject=.self)
        input.vals <- c(list(x=x), fun.args)
        true.val <- do.call(fun, input.vals)  # Concern: are we confident that the environment this is happening in is getting erased.
        
        # remove empty bins before noise is added (per definition of stability mechanism)
        true.val <- true.val[true.val > 0]
        
        scale <- sens / .self$epsilon
        release <- true.val + dpNoise(n=length(true.val), scale=scale, dist='laplace')
        
        # calculate the accuracy threshold, below which histogram buckets should be removed
        accuracyThreshold <- 1+2*log(2/delta)/epsilon
        # remove buckets below the threshold
        release <- release[release > accuracyThreshold]
        
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
    }
    
)