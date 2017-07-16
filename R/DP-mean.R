#' Function to evaluate the mean and specify parameters for mean functions
#'
#' @param x A numeric, logical, or integer vector.
#' @param var.type A character vector specifying variable type of \code{x}. 
#'    Should be of length one and should contain either 'numeric', 
#'    'logical', or 'integer'.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param sensitivity A numeric vector representing the probability of an arbitrary
#'    leakage of information from \code{x}. Should be of length one 
#'    and should be a very small value, such as 2^-30.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @return A list with fields `name` specifying the statistic and `stat` with 
#'    the value of the statistic.
dp.mean <- function(x, var.type, n, sensitivity, epsilon) {
    out <- list('name' = 'mean',
                'stat' = mean(x),
                'var.type' = var.type,
                'n' = n,
                'sensitivity' = sensitivity,
                'epsilon' = epsilon)
    return(out)
}

#' Function for differentially private release of mean
#'
#' @param x A numeric, logical, or integer vector.
#' @param var.type A character vector specifying variable type of \code{x}. 
#'    Should be of length one and should contain either 'numeric', 
#'    'logical', or 'integer'.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param rng A numeric vector specifying an a priori estimate of the range
#'    of \code{x}. Should be of length two. 
#' @return Differentially private mean of vector \code{x}.
#' @examples
#' n <- 1000
#' x_num <- runif(n)
#' x_int <- as.integer(round(x_num * 100))
#' x_bool <- x_num >= 0.5
#' x_dich <- ifelse(x_bool, 3.483, -9.657)
#' r_num <- mean.release(x=x_num, var.type='numeric', epsilon=0.5, n=n, range=c(0, 1))
#' r_int <- mean.release(x=x_int, var.type='integer', epsilon=0.5, n=n, range=c(5, 95))
#' r_bool <- mean.release(x=x_bool, var.type='logical', epsilon=0.5, n=n, range=c(0, 1))
#' r_dich <- mean.release(x=x_dich, var.type='logical', epsilon=0.5, n=n, range=c(-9.657, 3.483))
#' @export
mean.release <- function(x, var.type, n, epsilon, rng) {
    var.type <- check_variable_type(var.type, in_types=c('numeric', 'integer', 'logical'))
    postlist <- list('accuracy' = 'getAccuracy',
                     'epsilon' = 'getParameters',
                     'interval' = 'getCI')
    if (var.type == 'logical') {
        rng <- c(0, 1)
        postlist <- c(postlist, list('std' = 'postStandardDeviation',
                                     'median' = 'postMedian'))
    }
    rng <- checkrange(rng)
    sensitivity <- diff(rng) / n
    release <- mechanism.laplace(fun=dp.mean, x=x, var.type=var.type, rng=rng,
                                 sensitivity=sensitivity, epsilon=epsilon, n=n,
                                 postlist=postlist)
    return(release)
}

#' Postprocessed Standard Deviation for Logical Variables
#'
#' @param release Differentially private release of a mean for a logical 
#'    variable.
#' @return Standard deviation of \code{release}.
mean.postStandardDeviation <- function(release) {
    sd <- sqrt(release * (1 - release))
    return(sd)
}


#' Postprocessed Median for Logical Variables
#'
#' @param release Differentially private release of a mean for a logical 
#'    variable.
#' @return Median of \code{release}.
mean.postMedian <- function(release) {
    m <- ifelse(release < 0.5, 0, 1)
    return(m)
}


#' Function to find the accuracy guarantee of a mean release at a given epsilon 
#'    value. 
#'    
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the vector calculating the mean for.
#' @param alpha A numeric vector specifying the statistical significance level.
#' @return Accuracy guarantee for mean release given epsilon.
mean.getAccuracy <- function(epsilon, n, alpha=0.05) {
    accuracy <- log(1 / alpha) / (n * epsilon)
    return(accuracy)
}


#' Function to find the epsilon value necessary to meet a desired level of 
#'    accuracy.
#'
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param n A numeric vector of length one specifying the number of
#'    observations in the vector calculating the mean for.
#' @param alpha A numeric vector specifying the statistical significance level.
#' @return The scalar epsilon necessary to guarantee the needed accuracy. 
mean.getParameters <- function(accuracy, n, alpha=0.05) {
    epsilon <- log(1 / alpha) / (n * accuracy)
    return(epsilon)
}

#' Describe Here
#'
#' @param release something here
#' @param epsilon something here
#' @param sensitivity something here
#' @param n something here
#' @param range something here
#' @param alpha something here
#' @return Confidence bounds for differentially private release
mean.getCI <- function(release, epsilon, sensitivity, alpha=0.05) {
    z <- qlap((1 - (alpha / 2)), b=(sensitivity / epsilon))
    interval <- c(release - z, release + z)
    return(interval)
}


#' JSON doc for histogram
#'
#' @return JSON for histogram function
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


# --------------------------------------------------------- #
# --------------------------------------------------------- #
# Reference class implementation of mean

dpMean <- setRefClass(
    Class = 'dpMean',
    contains = c('mechanismLaplace')
)

dpMean$methods(
    initialize = function(mechanism, var.type, n, epsilon, rng, alpha=0.05) {
        .self$name <- 'Differentially private mean'
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$n <- n
        .self$epsilon <- epsilon
        .self$rng <- rng
        .self$alpha <- alpha
})


dpMean$methods(
    release = function(x) {
        sens <- diff(rng) / n
        .self$result <- export(mechanism)$evaluate(mean, x, sens, .self$postProcess)
})

dpMean$methods(
    postProcess = function(out) {
        out$accuracy <- mean.getAccuracy(epsilon, n, alpha)
        out$epsilon <- mean.getParameters(out$accuracy, n, alpha)
        out$interval <- mean.getCI(out$release, epsilon, (diff(rng) / n), alpha)
        if (var.type == 'logical') {
            out$std.dev <- mean.postStandardDeviation(out$release)
            out$median <- mean.postMedian(out$release)
        }
        return(out)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
