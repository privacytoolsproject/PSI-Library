#' Release additional model coefficients from DP covariance matrix
#' 
#' Function to extract regression coefficients using the differentially private covariance matrix 
#'    via the sweep operator. This is the function to obtain coefficients for additional models
#'    after the \code{covariance.release()} function has already been called. 
#' 
#' @param formula Formula, regression formula used on data
#' @param release Numeric, private release of covariance matrix
#' @param n Integer, indicating number of observations
#' @export
coefficient.release <- function(formula, release, n) {
  intercept <- ifelse('intercept' %in% names(release), TRUE, FALSE)
  coefficients <- linear.reg(formula, release, n, intercept)
  release <- list(name='Linear regression', 
                  n=n, 
                  formula=formula, 
                  coefficients=coefficients)
  return(release)
}


#' Function to get the sensitivity of the covariance matrix
#'
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data frame.
#' @param rng A numeric matrix of 2-tuples with the lower and upper bounds for
#'    each of P variables in the data frame, dimensions Px2.
#' @param intercept A logical vector of length one indicating whether an 
#'    intercept should be added prior to evaluating the inner product x'x.
#' @return The sensitivity of the data frame for which the covariance matrix
#'   is being calculated.
covariance.sensitivity <- function(n, rng, intercept) {
    diffs <- apply(rng, 1, diff)
    if (intercept) { diffs <- c(1, diffs) }
    sensitivity <- c()
    for (i in 1:length(diffs)) {
        for (j in i:length(diffs)) {
            s <- ((n - 1) / n) * diffs[i] * diffs[j]
            sensitivity <- c(sensitivity, s)
        }
    }
    return(sensitivity)
}


#' Function to convert unique private covariances into symmetric matrix
#'
#' @param release Differentially private release of elements in lower triangle 
#'    of covariance matrix.
#' @param columns A character vector indicating columns in the private 
#'    covariance to be included in the output. Length should be equal to the 
#'    number of columns the user wants to include.
#' @return A symmetric differentially private covariance matrix.
#'
#' This function is used by the mechanism in post-processing and not intended for interactive post-processing
covariance.formatRelease <- function(release, columns) {
    out.dim <- length(columns)
    out.matrix <- matrix(0, nrow=out.dim, ncol=out.dim)
    out.matrix[lower.tri(out.matrix, diag=TRUE)] <- release
    out.matrix[upper.tri(out.matrix, diag=FALSE)] <- t(out.matrix)[upper.tri(out.matrix, diag=FALSE)]
    release <- data.frame(out.matrix)
    rownames(release) <- names(release) <- columns
    return(release)
}


#' Function to perform linear regression using private release of covariance matrix
#'
#' @param release Differentially private release of elements in lower triangle 
#'    of covariance matrix.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data frame.
#' @param intercept Logical indicating whether the intercept is included in 
#'    \code{release}.
#' @param formula A list of the regression equations to be performed on the 
#'    covariance matrix.
#' @return Linear regression coefficients and standard errors for all specified
#'    \code{formula}.
covariance.postLinearRegression <- function(release, n, intercept, formula) {
    out.summaries <- vector('list', length(formula))
    for (f in 1:length(formula)) {
        out.summaries[[f]] <- linear.reg(formula[[f]], release, n, intercept)
    }
    return(out.summaries)
}


#' Lower triangle of ovariance matrix
#'
#' This function is called by \code{dpCovariance} and
#' produces the true value to be perturbed.
#'
#' @param x Data frame
#' @param columns Columns to include in output
#' @param intercept Logical, should the intercept be included?

fun.covar <- function(x, columns, intercept) {
    data <- x[, columns]
    if (intercept) { data <- cbind(1, data) }
    covariance <- t(as.matrix(data)) %*% as.matrix(data)
    lower <- lower.tri(covariance, diag=TRUE)
    covariance <- covariance[lower]
    return(covariance)
}


#' Differentially private covariance matrix
#'
#' @import methods
#' @export dpCovariance
#' @exportClass dpCovariance
#'
#' @include mechanisms.R

dpCovariance <- setRefClass(
    Class = 'dpCovariance',
    contains = 'mechanismLaplace'
)

dpCovariance$methods(
    initialize = function(mechanism, var.type, n, epsilon, columns, rng, impute.rng=NULL, 
                          intercept=FALSE, formula=NULL, delta=1e-5) {
        .self$name <- 'Differentially private covariance matrix'
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$n <- n
        .self$epsilon <- epsilon
        .self$delta <- delta
        .self$rng <- rng
        if (is.null(impute.rng)) {
            .self$impute.rng <- rng
        } else {
            .self$impute.rng <- impute.rng
        }
        .self$formula <- formula
        .self$intercept <- intercept
        if (.self$intercept) { 
            .self$columns <- c('intercept', columns)
        } else {
            .self$columns <- columns
        }
})

dpCovariance$methods(
    release = function(x) {
        sens <- covariance.sensitivity(n, rng, intercept)
        .self$result <- export(mechanism)$evaluate(fun=fun.covar, x=x, sens=sens, postFun=.self$postProcess, 
                                                   columns=columns, formula=formula, intercept=intercept)
})

dpCovariance$methods(
    postProcess = function(out, columns, formula, intercept) {
        out$release <- covariance.formatRelease(out$release, columns)
        out$variable <- columns
        if (!is.null(formula)) {
            out$linear.regression <- covariance.postLinearRegression(out$release, n, intercept, formula)
        }
        return(out)
})
