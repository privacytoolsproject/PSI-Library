#' Function to evaluate the mean and specify parameters for mean functions
#'
#' @param x Numeric vector
#' @return List with fields `name` specifying the statistic and `stat` with the value of the statistic

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
#' @param x Numeric vector 
#' @param epsilon Numeric
#' @param n The number of samples
#' @param range An a priori estimate of the range
#' @return Differentially private release of mean of vector x
#'
#' @examples
#'
#' n <- 1000
#' x_num <- runif(n)
#' x_int <- as.integer(round(x_num * 100))
#' x_bool <- x_num >= 0.5
#' x_dich <- ifelse(x_bool, 3.483, -9.657)
#' 
#' r_num <- mean.release(x=x_num, var.type='numeric', epsilon=0.5, n=n, range=c(0, 1))
#' r_int <- mean.release(x=x_int, var.type='integer', epsilon=0.5, n=n, range=c(5, 95))
#' r_bool <- mean.release(x=x_bool, var.type='logical', epsilon=0.5, n=n, range=c(0, 1))
#' r_dich <- mean.release(x=x_dich, var.type='logical', epsilon=0.5, n=n, range=c(-9.657, 3.483))

mean.release = function(x, var.type, n, epsilon, rng) {
    var.type <- check_variable_type(var.type, in_types=c('numeric', 'integer', 'logical'))
    if (var.type == 'logical') { rng = c(0, 1) }
    rng <- checkrange(rng)
    sensitivity <- diff(rng) / n
    release <- mechanism.laplace(
        fun=dp.mean,
        x=x,
        var.type=var.type,
        rng=rng,
        sensitivity=sensitivity,
        epsilon=epsilon,
        n=n)
    return(release)
}


#' @param epsilon Privacy parameter epsilon
#' @param n Number of observations
#' @param alpha The statistical significance level
#' @return Accuracy guarantee for mean release given epsilon

mean.getAccuracy = function(epsilon, n, alpha=0.05) {
    accuracy <- log(1 / alpha) / (n * epsilon)
    return(accuracy)
}


#' @param accuracy The accuracy we need to guarantee (percent)
#' @param n The number of samples
#' @param alpha The statistical signifcance level
#' @return The scalar epsilon necessary to guarantee the accuracy needed

mean.getParameters = function(accuracy, n, alpha=0.05) {
    epsilon <- log(1 / alpha) / (n * accuracy)
    return(epsilon)
}


#' @param release
#' @param epsilon
#' @param sensitivity
#' @param n
#' @param range
#' @param alpha
#' @return Confidence bounds for differentially private release

mean.getCI = function(release, epsilon, sensitivity, n, rng, alpha=0.05) {
    z <- qlap((1 - alpha), b=(sensitivity / epsilon))
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
