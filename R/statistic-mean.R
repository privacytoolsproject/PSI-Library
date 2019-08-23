#' Differentially private mean
#'
#' @param mechanism Character, the privacy mechanism. For \code{dpMean}, one
#'   of \code{c('mechanismLaplace', 'mechanismBootstrap')}.
#' @param varType Character, the R variable type. One of \code{c('numeric',
#'   'integer', 'logical')}.
#' @param Variable Character, variable name.
#' @param n Integer, number of observations
#' @param rng Numeric, a priori estimate of the range
#' @param epsilon Numeric, privacy cost parameter
#' @param accuracy Numeric, accuracy guarantee given \code{epsilon}
#' @param imputeRng Numeric, range within which missing values are imputed. If \code{NULL},
#'   the range provided in \code{rng} is used.
#' @param alpha Numeric, the level of statistical significance. Default 0.05.
#' @param nBoot Integer, the number of bootstrap replications if using the bootstrap
#'   bootstrap mechanism, ignored otherwise. Default 20.
#'
#' @import methods
#' @export dpMean
#' @exportClass dpMean
#'
#' @include mechanism.R
#' @include mechanism-laplace.R
#' @include mechanism-bootstrap.R

dpMean <- setRefClass(
    Class = 'dpMean',
    contains = c('mechanismLaplace', 'mechanismBootstrap')
)

dpMean$methods(
    initialize = function(mechanism, varType, variable, n, rng=NULL, epsilon=NULL,
                          accuracy=NULL, imputeRng=NULL, alpha=0.05, nBoot=20, ...) {
        .self$name <- 'Differentially private mean'
        .self$mechanism <- mechanism
        .self$varType <- varType
        .self$variable <- variable
        .self$n <- checkNValidity(n)
        .self$alpha <- alpha
        .self$rng <- checkRange(rng, varType)
        .self$sens <- meanSensitivity(.self$rng, .self$n)
        
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- laplaceGetEpsilon(.self$sens, .self$accuracy, alpha)
        } else {
            checkEpsilon(epsilon)
            .self$epsilon <- epsilon
            .self$accuracy <- laplaceGetAccuracy(.self$sens, .self$epsilon, alpha)
        }
        
        if (is.null(imputeRng)) {
            .self$imputeRng <- .self$rng
        } else {
            .self$imputeRng <- checkImputationRange(imputationRange=imputeRng, rng=.self$rng, varType=.self$varType)
        }
        
        .self$bootFun <- bootMean
        .self$nBoot <- nBoot
})


dpMean$methods(
    #' Differentially private mean release
    #' 
    #' @name dpMeanRelease
    #' @param data Dataframe with a column named .self$variable, where
    #'  that column has data of type .self$varType and which is bounded by 
    #'  .self$rng.
    #'
    #' Assigns to .self$result a dataframe that describes the differentially private 
    #' mean, calculated by some mechanism as defined in .self$mechanism, of that
    #' column of the dataframe and any post-processing on that output. This 
    #' postprocessing is done in @seealso{dpMean$postProcess}
    #' 
    #' Note that the actual differentially private release is calculated in a call to the
    #' differentially private mechanism .self$mechanism's \code{evaluate} function within 
    #' the \code{dpMean$release} function. 
    release = function(data, ...) {
        x <- data[, variable]
        .self$result <- export(mechanism)$evaluate(mean, x, .self$sens, .self$postProcess, ...)
})

dpMean$methods(
    #' Post-processing on differentially private mean, called within the \code{.self$mechanism$evaluate}
    #' function, which in turn is called within \code{dpMean$release}. 
    #' 
    #' @name dpMeanPostProcess
    #' @param out Input dataframe that describes the differentially private release that was created in
    #' \code{.self$mechanism$evaluate}. This dataframe will have at least one pre-existing attribute, 
    #' out$release, which is a numeric value of length one that is the differentially private mean.
    #' 
    #' This function is able to calculate a confidence interval on the output. If the variable was logical,
    #' standard deviation, median, and a histogram may be additionally computed as post-processing steps.
    #' 
    #' Additionally, known portions of the input such as the variable name 
    #' and the epsilon value may be appended with no extra privacy loss.
    #'
    #' @return Dataframe \code{out}, updated to include post-processed values.
    postProcess = function(out) {
        out$variable <- variable
        if (mechanism == 'mechanismLaplace') {
            out$accuracy <- accuracy
            out$epsilon <- epsilon
            out$interval <- meanGetCI(out$release, epsilon, .self$sens, alpha)
        } 
        if (varType == 'logical') {
            if (mechanism == 'mechanismBootstrap') {
                baggedEstimate <- mean(out$release)
                out$stdDev <- meanPostStandardDeviation(baggedEstimate)
                out$median <- meanPostMedian(baggedEstimate)
                out$histogram <- meanPostHistgram(baggedEstimate)
            } else {
                out$stdDev <- meanPostStandardDeviation(out$release)
                out$median <- meanPostMedian(out$release)
                out$histogram <- meanPostHistogram(out$release, n)
            }
        }
        return(out)
})
