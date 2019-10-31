#' Sensitivity of mean
#'
#' For a detailed derivation of the sensitivity, see /extra_docs/sensitivities/mean_sensitivity.pdf
#' 
#' @param rng Numeric vector of length two; first entry is minimal bound on the database entries, second is maximal bound on the database entries.
#' @param n Numeric vector of length one; the number of datapoints in the database.
#'
#' @return Numeric vector of length one; a maximal bound on the sensitivity of the population variance.
#'
#' @examples
#' meanSensitivity(c(0,10),5) #should return 2
meanSensitivity <- function(rng, n){
  return(diff(rng)/n)
}

#' Postprocessed standard deviation for logical variables 
#' 
#' Calculates the standard deviation of the differentially private mean from a 
#' logical variable.
#'
#' @param release Differentially private release of a mean for a logical 
#'    variable.
#'    
#' @return Standard deviation of the logical variable
#' @rdname meanPostStandardDeviation

meanPostStandardDeviation <- function(release) {
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
#' @rdname meanPostMedian

meanPostMedian <- function(release) {
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
#' @rdname meanPostHistogram

meanPostHistogram <- function(release, n) {
    ones <- round(release * n)
    histogram <- data.frame(matrix(c(n - ones, ones), ncol=2))
    names(histogram) <- c(0, 1)
    return(histogram)
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
#' @return Confidence bounds for differentially private mean release.
#'
#' @export meanGetCI
#' @rdname meanGetCI

meanGetCI <- function(release, epsilon, sensitivity, alpha=0.05) {
    z <- qLap((1 - (alpha / 2)), b=(sensitivity / epsilon))
    interval <- c(release - z, release + z)
    return(interval)
}


#' JSON doc for mean
#' 
#' Produce a JSON doc for differentially private means.
#' 
#' @param outputJSON Should the output be converted to JSON format. Default
#' to \code{TRUE}.
#'
#' @return JSON for mean function
#' @rdname meanGetJSON

meanGetJSON <- function(outputJSON=TRUE) {
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
    if (outputJSON) {
        out <- jsonlite::toJSON(out, pretty=TRUE)
    }
    return(out)
}



#' Bootstrap mean function
#'
#' This function is used when the bootstrap mechanism is used
#'
#' @param xi Vector of values
#' @param n Number of observations
#' @return Mean

bootMean <- function(xi, n) {
    return(sum(xi) / n)
}


