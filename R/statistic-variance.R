#' Postprocessed variance standard deviation
#'
#' Function to extract standard deviation from noisy estimate of variance.
#'
#' @param release Differentially private release of variance.
#'
#' @return Noisy estimate of the standard deviation of \code{release}.
#' @rdname variance.postStandardDeviation

variance.postStandardDeviation <- function(release) {
    std <- sqrt(release)
    return(std)
}


#' Differentially private variance
#'
#' @param mechanism Character, the privacy mechanism. For \code{dpVariance}, one
#'                  of \code{c('mechanismLaplace', 'mechanismSnapping')}.
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
#' @param gamma Numeric, used to provide high probability (1-gamma) guarantee that clamping bound does not bind. Default 0.05. Used only for \code{mechanismSnapping}.
#' @import methods
#' @export dpVariance
#' @exportClass dpVariance
#'
#' @include mechanism.R
#' @include mechanism-laplace.R
#' @include mechanism-snapping.R

dpVariance <- setRefClass(
    Class = 'dpVariance',
    contains = c('mechanismLaplace', 'mechanismSnapping'),
    fields = list(gamma = 'numeric')
)

dpVariance$methods(
    initialize = function(mechanism, var.type, variable, n, rng=NULL, epsilon, accuracy = NULL,
                          impute.rng=NULL, alpha=0.05, gamma=0.05, ...) {
        .self$name <- 'Differentially private variance'
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$variable <- variable
        .self$n <- check_n_validity(n)
        .self$rng <- checkrange(rng, var.type)
        .self$sens <- (n - 1) / n^2 * diff(rng)^2
        .self$alpha <- alpha
        .self$gamma <- gamma

        .self$min_B <- (abs(.self$rng[1] - .self$rng[2])/2)^2

        if (mechanism == 'mechanismSnapping') {
            reticulate::source_python(system.file('python', 'cc_snap.py', package = 'PSIlence'))
            # create dummy version of Snapping Mechanism object in order to get epsilon/accuracy guarantees
            snap_obj <- Snapping_Mechanism(mechanism_input = 0,
                                            sensitivity = .self$sens,
                                            epsilon = epsilon,
                                            accuracy = accuracy,
                                            alpha = .self$alpha,
                                            min_B = .self$min_B,
                                            gamma = .self$gamma)
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
            }}

        if (is.null(impute.rng)) {
            .self$impute.rng <- rng
        } else {
            .self$impute.rng <- checkImputationRange(impute.rng, .self$rng, .self$var.type)
        }
})

dpVariance$methods(
    release = function(data, ...) {
        x <- data[, variable]
        .self$result <- export(mechanism)$evaluate(var, x, .self$sens, .self$postProcess, ...)
})

dpVariance$methods(
    postProcess = function(out) {
        out$variable <- variable
        out$std <- sqrt(out$release)
        return(out)
})
