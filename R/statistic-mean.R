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
#' @export mean.getCI
#' @rdname mean.getCI

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
#' @param mechanism Character, the privacy mechanism. For \code{dpMean}, one
#'   of \code{c('mechanismLaplace', 'mechanismBootstrap', 'mechanismSnapping')}.
#' @param var.type Character, the R variable type. One of \code{c('numeric',
#'   'integer', 'logical')}.
#' @param Variable Character, variable name.
#' @param n Integer, number of observations
#' @param rng Numeric, a priori estimate of the range
#' @param epsilon Numeric, privacy cost parameter
#' @param accuracy Numeric, accuracy guarantee given \code{epsilon}
#' @param impute.rng Numeric, range within which missing values are imputed. If \code{NULL},
#'   the range provided in \code{rng} is used.
#' @param alpha Numeric, the level of statistical significance. Default 0.05.
#' @param n.boot Integer, the number of bootstrap replications if using the bootstrap
#'   bootstrap mechanism, ignored otherwise. Default 20.
#' @param gamma Numeric, used to provide high probability (1-gamma) guarantee that clamping bound does not bind. Default 0.05. Used only for \code{mechanismSnapping}.
#'
#' @import methods
#' @export dpMean
#' @exportClass dpMean
#'
#' @include mechanism-snapping.R
#' @include mechanism-laplace.R
#' @include mechanism-bootstrap.R

dpMean <- setRefClass(
    Class = 'dpMean',
    contains = c('mechanismLaplace', 'mechanismBootstrap', 'mechanismSnapping')
)

dpMean$methods(
    initialize = function(mechanism, var.type, variable, n, rng=NULL, epsilon=NULL,
                          accuracy=NULL, impute.rng=NULL, alpha=0.05, n.boot=20, gamma = 0.05, ...) {
        .self$name <- 'Differentially private mean'
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$variable <- variable
        .self$n <- check_n_validity(n)
        .self$alpha <- alpha
        .self$rng <- checkrange(rng, var.type)
        .self$sens <- diff(.self$rng) / n
        .self$gamma <- gamma

        min_B <- max(abs(.self$rng[1]), abs(.self$rng[2]))

        if (mechanism == 'mechanismSnapping') {
            reticulate::source_python(system.file('python', 'cc_snap.py', package = 'PSIlence'))
            # create dummy version of Snapping Mechanism object in order to get epsilon/accuracy guarantees
            snap_obj <- Snapping_Mechanism(mechanism_input = 0,
                                            sensitivity = .self$sens,
                                            epsilon = epsilon,
                                            accuracy = accuracy,
                                            alpha = alpha,
                                            min_B = min_B,
                                            gamma = gamma)
        }

        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            if (mechanism == 'mechanismLaplace'){
                .self$epsilon <- laplace.getEpsilon(.self$sens, .self$accuracy, alpha)
            } else if (mechanism == 'mechanismSnapping'){
                .self$epsilon <- snap_obj$epsilon
            }
        } else {
            checkepsilon(epsilon)
            .self$epsilon <- epsilon
            if (mechanism == 'mechanismLaplace'){
                .self$accuracy <- laplace.getAccuracy(.self$sens, .self$epsilon, alpha)
            } else if (mechanism == 'mechanismSnapping'){
                .self$accuracy <- snap_obj$accuracy
            }
        }
        if (is.null(impute.rng)) {
            .self$impute.rng <- .self$rng
        } else {
            .self$impute.rng <- checkImputationRange(imputationRange=impute.rng, rng=.self$rng, var.type=.self$var.type)
        }
        .self$boot.fun <- boot.mean
        .self$n.boot <- n.boot
})


dpMean$methods(
    release = function(data, ...) {
        x <- data[, variable]
        .self$result <- export(mechanism)$evaluate(mean, x, .self$sens, .self$postProcess, ...)
})

dpMean$methods(
    postProcess = function(out) {
        out$variable <- variable
        if (mechanism == 'mechanismLaplace') {
            out$accuracy <- accuracy
            out$epsilon <- epsilon
            out$interval <- mean.getCI(out$release, epsilon, .self$sens, alpha)
        } else if (mechanism == 'mechanismSnapping') {
            out$accuracy <- accuracy
            out$epsilon <- epsilon
        }
        if (var.type == 'logical') {
            if (mechanism == 'mechanismBootstrap') {
                bagged.estimate <- mean(out$release)
                out$std.dev <- mean.postStandardDeviation(bagged.estimate)
                out$median <- mean.postMedian(bagged.estimate)
                out$histogram <- mean.postHistogram(bagged.estimate)
            } else {
                out$std.dev <- mean.postStandardDeviation(out$release)
                out$median <- mean.postMedian(out$release)
                out$histogram <- mean.postHistogram(out$release, n)
            }
        }
        return(out)
})
