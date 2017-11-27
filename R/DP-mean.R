#' Postprocessed standard deviation for logical variables 
#' 
#' Calculates the standard deviation of the differentially private mean from a 
#' logical variable.
#'
#' @param release Differentially private release of a mean for a logical 
#'    variable.
#'    
#' @return Standard deviation of the logical variable
#' @rdname mean.postStandardDeviation
mean.postStandardDeviation <- function(release) {
    sd <- sqrt(release * (1 - release))
    return(sd)
}


#' Postprocessed median for logical variables
#'
#' Calculates the median of the differentially private mean from a 
#' logical variable.
#'
#' @param release Differentially private release of a mean for a logical 
#'    variable.
#'
#' @return Median of the logical variable
#' @rdname mean.postMedian
mean.postMedian <- function(release) {
    m <- ifelse(release < 0.5, 0, 1)
    return(m)
}


#' Postprocessed histogram for logical variables
#'
#' Generate counts for levels of a logical variable based on the release
#'
#' @param release Numeric private mean
#' @param n Integer indicating number of observations
#'
#' @return Data frame, histogram of the logical variable
#' @rdname mean.postHistogram
mean.postHistogram <- function(release, n) {
    ones <- round(release * n)
    histogram <- data.frame(matrix(c(n - ones, ones), ncol=2))
    names(histogram) <- c(0, 1)
    return(histogram)
}


#' Mean accuracy
#' 
#' Function to find the accuracy guarantee of a mean release at a given epsilon 
#' value. 
#'    
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the vector calculating the mean for.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @export mean.getAccuracy
#' @return Accuracy guarantee for mean release given epsilon.
#' @rdname mean.getAccuracy
mean.getAccuracy <- function(epsilon, n, alpha=0.05, rng) {
	range <- max(rng)-min(rng)
    accuracy <- log(1 / alpha)*range/ (n * epsilon)
    return(accuracy)
}


#' Mean epsilon
#' 
#' Function to find the epsilon value necessary to meet a desired level of 
#' accuracy for a mean release.
#'
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param n A numeric vector of length one specifying the number of
#'    observations in the vector calculating the mean for.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @export mean.getParameters
#' @return The scalar epsilon necessary to guarantee the needed accuracy.
#' @rdname mean.getParameters
mean.getParameters <- function(accuracy, n, alpha=0.05, rng) {
	range <- max(rng)-min(rng)
    epsilon <- log(1 / alpha)*range / (n * accuracy)
    return(epsilon)
}


#' Mean confidence interval
#' 
#' Return the confidence interval for the differentially private mean release given the
#' accuracy.
#'
#' @param release Differentially private release of a mean.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param sensitivity The difference of the range of the data divided 
#'    by \code{n}.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @return Confidence bounds for differentially private mean release.
mean.getCI <- function(release, epsilon, sensitivity, alpha=0.05) {
    z <- qlap((1 - (alpha / 2)), b=(sensitivity / epsilon))
    interval <- c(release - z, release + z)
    return(interval)
}


#' JSON doc for mean
#' 
#' Produce a JSON doc for differentially private means.
#' 
#' @param output.json Should the output be converted to JSON format. Default
#' to \code{TRUE}.
#'
#' @return JSON for mean function
#' @rdname mean.getJSON
mean.getJSON <- function(output.json=TRUE) {
    out <- list()
    out$statistic <- 'Mean'
    out$description <- 'Differentially Private Mean'
    out$mechanisms <- c('Laplace')
    out$variableTypes <- list('numeric' = list(), 'categorical' = list())
    out$variableTypes$numeric$rTypes <- c('numeric', 'integer')
    out$variableTypes$numeric$fields <- list(
        'n' = 'Number of observations',
        'range' = 'Ordered pair indicating effective lower and upper bounds'
    )
    out$variableTypes$categorical$rTypes <- c('logical')
    out$variableTypes$categorical$fields <- list(
        'n' = 'Number of observations',
        'range' = 'Should be (0, 1)'
    )
    if (output.json) {
        out <- jsonlite::toJSON(out, pretty=TRUE)
    }
    return(out)
}



#' Mean function
#'
#' This function is used when the bootstrap mechanism is used
#'
#' @param xi Vector of values
#' @param n Number of observations
#' @return Mean

boot.mean <- function(xi, n) {
    return(sum(xi) / n)
}


#' Differentially private mean
#'
#' @import methods
#' @export dpMean
#' @exportClass dpMean
#'
#' @include mechanisms.R

dpMean <- setRefClass(
    Class = 'dpMean',
    contains = c('mechanismLaplace', 'mechanismBootstrap')
)

dpMean$methods(
    initialize = function(mechanism, var.type, n, rng, epsilon=NULL, accuracy=NULL, 
                          impute.rng=NULL, alpha=0.05, boot.fun=boot.mean, ...) {
        .self$name <- 'Differentially private mean'
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$n <- n
        .self$alpha <- alpha
        .self$rng <- rng
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- mean.getParameters(accuracy, n, alpha, rng)
        } else {
            .self$epsilon <- epsilon
            .self$accuracy <- mean.getAccuracy(epsilon, n, alpha, rng)
        }
        if (is.null(impute.rng)) {
            .self$impute.rng <- rng
        } else {
            .self$impute.rng <- impute.rng
        }
        .self$boot.fun <- boot.fun
})


dpMean$methods(
    release = function(x, ...) {
        sens <- diff(rng) / n
        .self$result <- export(mechanism)$evaluate(mean, x, sens, .self$postProcess, ...)
})

dpMean$methods(
    postProcess = function(out) {
        if (mechanism == 'mechanismLaplace') {
            out$accuracy <- accuracy
            out$epsilon <- epsilon
            out$interval <- mean.getCI(out$release, epsilon, (diff(rng) / n), alpha)
        } 
        if (var.type == 'logical') {
            if (mechanism == 'mechanismBootstrap') {
                bagged.estimate <- mean(out$release)
                out$std.dev <- mean.postStandardDeviation(bagged.estimate)
                out$median <- mean.postMedian(bagged.estimate)
                out$histogram <- mean.postHistgram(bagged.estimate)
            } else {
                out$std.dev <- mean.postStandardDeviation(out$release)
                out$median <- mean.postMedian(out$release)
                out$histogram <- mean.postHistogram(out$release, n)
            }
        }
        return(out)
})
