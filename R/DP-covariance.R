#' Function to evaluate the covariance matrix from input matrix and specify parameters for post-processing
#'
#' @param x Numeric data frame
#' @param n Integer indicating the number of observations
#' @param rng Numeric matrix of 2-tuples with the lower and upper bounds for each of P variables, dimensions Px2
#' @param epsilon Numeric differential privacy parameter epsilon
#' @param columns Character vector indicating columns in \code{x}
#' @param intercept Logical indicating whether an intercept should be added prior to evaluating the inner product x'x

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


#' Function to release a differentially private inner product of input matrix
#'
#' @param x Numeric data frame, dimensions NxP
#' @param var.type Numeric variable types
#' @param n Integer indicating the number of observations
#' @param epsilon Numeric differential privacy parameter
#' @param rng Numeric matrix of 2-tuples with the lower and upper bounds for each of P variables, dimensions Px2
#' @param columns Character vector indicating columns in \code{x}
#' @param delta something here
#' @param intercept Logical indicating whether an intercept should be added prior to evaluating the inner product x'x, default to FALSE
#' @export

covariance.release <- function(x, var.type, n, epsilon, rng, columns, delta=0.000001, intercept=FALSE, formulae=NULL, mechanism='laplace') {

    sensitivity <- covariance.sensitivity(n, rng, intercept)
    if (is.null(formulae)) {
        formulae <- as.list(formulae)
    } else {
        if (inherits(formulae, 'formula')) { formulae <- list(formulae) }
    }

    # pass to mechanism
    postlist <- list('release' = 'formatRelease')
    if (length(formulae)!=0) {
        postlist <- c(postlist, list('linear.regression' = 'postLinearRegression'))
    }
    if (mechanism=='laplace') {
      release <- mechanism.laplace(fun=dp.covariance, x=x, var.type='numeric', n=n,
                                   rng=rng, epsilon=(epsilon / length(sensitivity)), columns=columns,
                                   sensitivity=sensitivity, intercept=intercept, formulae=formulae,
                                   postlist=postlist)
    } else if (mechanism=='gaussian') {
      release <- mechanism.gaussian(fun=dp.covariance, x=x, var.type='numeric', n=n, rng=rng, sensitivity=sensitivity, 
                                    epsilon=(epsilon / length(sensitivity)), columns=columns,
                                    delta=delta, intercept=intercept, formulae=formulae, postlist=postlist)
    } else {stop('no noise mechanism defined')}
    return(release)
}


#' Function to get the sensitivity of the covariance matrix
#'
#' @param n Integer indicating the number of observations
#' @param rng Numeric matrix of 2-tuples with the lower and upper bounds for each of P variables, dimensions Px2
#' @param intercept Logical indicating whether an intercept should be included

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
#' @param release Numeric private release of elements in lower triangle of covariance matrix
#' @param columns Character vector with column names
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
#' @param release Numeric private release of covariance matrix
#' @param n Integer indicating number of observations
#' @param intercept Logical indicating whether the intercept is included
#' @param formulae List of R formulae

covariance.postLinearRegression <- function(release, n, intercept, formulae) {
    out.summaries <- vector('list', length(formulae))
    for (f in 1:length(formulae)) {
        out.summaries[[f]] <- linear.reg(formulae[[f]], release, n, intercept)
    }
    return(out.summaries)
}
       
# --------------------------------------------------------- #
# --------------------------------------------------------- #
# Reference class for covariance with Laplace noise

fun.covar <- function(x, columns) {
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
    initialize = function(mechanism, var.type, n, epsilon, rng) {
        .self$name <- 'Differentially private covariance matrix'
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$n <- n
        .self$epsilon <- epsilon
        .self$rng
})

dpCovariance$methods(
    release = function(x, columns, formulae=NULL, intercept=FALSE) {
        if (intercept) { columns <- c('intercept', columns) }
        sens <- covariance.sensitivity(rng, n, intercept)
        .self$result <- export(mechanism)$evaluate(fun=fun.covar, x=x, sensitivity=sens, postFun=.self$postProcess,
                                                   columns=columns, formulae=formulae, intercept=intercept)
})

dpCovariance$methods(
    postProcess = function(out, columns, formulae, intercept) {
        out$release <- covariance.formatRelease(out$release, columns)
        out$linear.regression <- covariance.postLinearRegression(out$release, n, intercept, formulae)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
