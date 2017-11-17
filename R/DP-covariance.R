#' Differentially Private Covariance Matrix
#' 
#' Function to evaluate the covariance matrix and specify parameters for covariance functions.
#' 
#' @param x A numeric data frame with at least two columns.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param rng A numeric matrix of 2-tuples with the lower and upper bounds for
#'    each of P variables in \code{x}, dimensions Px2. 
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param columns A character vector indicating columns in \code{x} to be 
#'    included in the covariance matrix. Length should be equal to the number 
#'    of columns the user wants to include.
#' @param intercept A logical vector of length one indicating whether an 
#'    intercept should be added prior to evaluating the inner product x'x.
#' @param formulae The regression equations the user would like to perform on
#'    the covariance matrix. The equations should be of class 'formula'. The 
#'    user may specify as many equations as desired.
#' @return A list with fields `name` specifying the statistic and `stat` with 
#'    the lower triangle of the covariance matrix.
#' @export
dp.covariance <- function(x, n, rng, epsilon, columns, intercept, formulae) {

    # subset and optionally append an intercept
    data <- x[, columns]
    if (intercept) {
        data <- cbind(1, data)
        columns <- c('intercept', columns)
    }

    # the true statistic, return only the lower triangle
    covariance <- t(as.matrix(data)) %*% as.matrix(data)
    lower <- lower.tri(covariance, diag=TRUE)
    covariance <- covariance[lower]

    # specif
    output <- list('name' = 'covariance',
                   'stat' = covariance,
                   'n' = n,
                   'rng' = rng,
                   'epsilon' = epsilon,
                   'columns' = columns,
                   'intercept' = intercept)
    if (!is.null(formulae)) {
        output[['formulae']] <- formulae
    }
    return(output)
}


#' Function to release a differentially private inner product of input matrix.
#'
#' @param x A numeric data frame with at least two columns.
#' @param var.type A character vector specifying variable type of \code{x}. 
#'    Should be of length one and should contain either 'numeric', 
#'    'logical', or 'integer'.
#' @param n A numeric vector of length one specifying the number of
#'    observations in \code{x}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param rng A numeric matrix of 2-tuples with the lower and upper bounds for
#'    each of P variables in \code{x}, dimensions Px2. 
#' @param columns A character vector indicating columns in \code{x} to be 
#'    included in the covariance matrix. Length should be equal to the number 
#'    of columns the user wants to include.
#' @param impute.rng Numeric matrix of ranges for each column within which 
#'    missing values are imputed. If \code{NULL}, defaults to \code{rng}.
#' @param delta A numeric vector representing the probability of an arbitrary
#'    leakage of information from \code{x}. Should be of length one 
#'    and should be a very small value. Default to 10^-6.
#' @param intercept A logical vector of length one indicating whether an 
#'    intercept should be added prior to evaluating the inner product x'x.
#'    Default to FALSE.
#' @param formulae The regression equations the user would like to perform on
#'    the covariance matrix. The equations should be of class 'formula'. The 
#'    user may specify as many equations as desired. Default to FALSE.
#' @param mechanism A character vector specifying the mechanism used to apply 
#'    noise to the covariance matrix. Should be of length one and contain 
#'    either 'laplace', 'gaussian', or 'wishart'. Default to 'laplace'.
#' @return Differentially private covariance matrix of \code{x}.
#' @examples
#' 
#' data(PUMS5extract10000, package = "PSIlence")
#' data <- data.frame(income = PUMS5extract10000$income, education = PUMS5extract10000$educ)
#' range.income <- c(-10000, 713000)
#' range.education <- c(1, 16)
#' range <- rbind(range.income, range.education)
#' covariance.release(x = data, var.type = 'numeric', n = 10000, epsilon = 0.2, rng = range, 
#'                    columns = c("income", "education"), formulae = income ~ education)
#' @export
covariance.release <- function(x, var.type, n, epsilon, rng, columns, impute.rng=NULL, 
                               delta=0.000001, intercept=FALSE, formulae=NULL, mechanism='laplace') {

    sensitivity <- covariance.sensitivity(n, rng, intercept)
    if (is.null(formulae)) {
        formulae <- as.list(formulae)
    } else {
        if (inherits(formulae, 'formula')) { formulae <- list(formulae) }
    }
    if (is.null(impute.rng)) { impute.rng <- rng }

    # pass to mechanism
    postlist <- list('release' = 'formatRelease')
    if (length(formulae) != 0) {
        postlist <- c(postlist, list('linear.regression' = 'postLinearRegression'))
    }
    # Note: shouldn't pass function by scopeing
    if (mechanism == 'laplace') {
        release <- mechanism.laplace(fun=dp.covariance, x=x, var.type='numeric', n=n,
                                     rng=rng, epsilon=(epsilon / length(sensitivity)), columns=columns,
                                     sensitivity=sensitivity, intercept=intercept, formulae=formulae,
                                     impute.rng=impute.rng, postlist=postlist)
    } else if (mechanism == 'gaussian') {
        release <- mechanism.gaussian(fun=dp.covariance, x=x, var.type='numeric', n=n, rng=rng, 
                                      sensitivity=sensitivity, epsilon=(epsilon / length(sensitivity)), 
                                      columns=columns, delta=delta, intercept=intercept, formulae=formulae, 
                                      impute.rng=impute.rng, postlist=postlist)
    } else {
        stop('no noise mechanism defined')
    }
    return(release)
}

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
  intercept <- ifelse('intercept' %in% names(release), T, F)
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
       
# --------------------------------------------------------- #
# --------------------------------------------------------- #
# Reference class for covariance with Laplace noise

fun.covar <- function(x, columns, intercept) {
    data <- x[, columns]
    if (intercept) { data <- cbind(1, data) }
    covariance <- t(as.matrix(data)) %*% as.matrix(data)
    lower <- lower.tri(covariance, diag=TRUE)
    covariance <- covariance[lower]
    return(covariance)
}

#'
#' @include mechanisms.R

dpCovariance <- setRefClass(
    Class = 'dpCovariance',
    contains = c('mechanismLaplace')
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
        if (!is.null(formula)) {
            out$linear.regression <- covariance.postLinearRegression(out$release, n, intercept, formula)
        }
        return(out)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
