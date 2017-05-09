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
#' @param intercept Logical indicating whether an intercept should be added prior to evaluating the inner product x'x, default to FALSE

covariance.release <- function(x, var.type, n, epsilon, rng, columns, intercept=FALSE, formulae=NULL) {

    # get the vector of sensitivities from the ranges
    diffs <- apply(rng, 1, diff)
    if (intercept) { diffs <- c(-1, diffs) }
    sensitivity <- c()
    for (i in 1:length(diffs)) {
        for (j in i:length(diffs)) {
            s <- ((n - 1) / n) * diffs[i] * diffs[j]
            sensitivity <- c(sensitivity, s)
        }
    }

    # pass to mechanism
    postlist <- list('release' = 'formatRelease')
    if (!is.null(formulae)) {
        postlist <- c(postlist, list('linear.regression' = 'postLinearRegression'))
    }
    release <- mechanism.laplace(fun=dp.covariance, x=x, var.type='numeric', n=n,
                                 rng=rng, epsilon=(epsilon / length(sensitivity)), columns=columns,
                                 sensitivity=sensitivity, intercept=intercept, formulae=as.list(formulae),
                                 postlist=postlist)
    return(release)
}


#' Function to convert unique private covariances into symmetric matrix
#'
#' @param release Numeric private release of elements in lower triangle of covariance matrix
#' @param intercept Logical indicating whether the intercept is included
#' @param columns Character vector with column names
#'
#' This function is used by the mechanism in post-processing and not intended for interactive post-processing

covariance.formatRelease <- function(release, intercept, columns) {
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
        out.summaries[[f]] <- regress(formulae[[f]], release, n, intercept)
    }
    return(out.summaries)
}


#' Function to perform regression using the covariance matrix via the sweep operator
#'
#' @param formula Formula
#' @param release Numeric private release of covariance matrix
#' @param n Integer indicating number of observations
#' @param intercept Logical indicating whether the intercept is included

regress <- function(formula, release, n, intercept) {
    xy.locs <- extract.indices(formula, release, intercept)
    x.loc <- xy.locs$x.loc
    y.loc <- xy.locs$y.loc
    loc.vec <- rep(TRUE, (length(x.loc) + 1))
    loc.vec[y.loc] <- FALSE
    sweep <- amsweep((as.matrix(release) / n), loc.vec)
    coefs <- sweep[y.loc, x.loc]
    se <- sqrt(sweep[y.loc, y.loc] * diag(solve(release[x.loc, x.loc])))
    coefs <- data.frame(cbind(coefs, se))
    coefs <- format(round(coefs, 5), nsmall=5)
    rownames(coefs) <- xy.locs$x.label
    names(coefs) <- c('Estimate', 'Std. Error')
    return(coefs)
}


#' Function to obtain indices in data frame for dependent & independent variables from a formula
#'
#' @param formula Formula
#' @param data Data frame, in this case being the data frame of a private covariance matrix
#' @param intercept Logical indicating whether the intercept is included

extract.indices <- function(formula, data, intercept) {
    t <- terms(formula, data=data)
    y.loc <- attr(t, 'response')
    x.loc <- which(names(data) %in% attr(t, 'term.labels'))
    x.label <- names(data)[x.loc]
    if (intercept) {
        intercept.loc <- which(names(data) == 'intercept')
        x.loc <- c(intercept.loc, x.loc)
        x.label <- append(x.label, 'Intercept', after=(intercept.loc - 1))
        if (intercept.loc <= y.loc) { y.loc <- y.loc + 1 }
    }
    return(list('y.loc' = y.loc,
                'x.loc' = x.loc,
                'x.label' = x.label))
}


if (interactive()) {

    set.seed(681521)
    n <- 10000
    x1 <- rnorm(n, mean=73, sd=17)
    x2 <- rpois(n, lambda=4)
    x3 <- as.integer(rnorm(n, mean=25, sd=4))
    y <- 1.276 + 0.145 * x1 - 0.565 * x2 + 0.013 * x3 + rnorm(n, sd=1.5)
    df <- as.data.frame(cbind(y, x1, x2, x3))
    cols <- names(df)

    rng <- rbind(range(y), range(x1), range(x2), range(x3))
    dtypes <- 'numeric'
    eps <- 1.4

    # covariance example
    release <- covariance.release(df, dtypes, n, eps, rng, cols)
    observe <- t(as.matrix(df)) %*% as.matrix(df)

    # example with regression
    form1 <- as.formula('y ~ x1 + x2 + x3')
    form2 <- as.formula('y ~ x1 + x2')
    release2 <- covariance.release(df, dtypes, n, eps, rng, cols, intercept=TRUE, formulae=c(form1, form2))
    df2 <- cbind(1, df)
    names(df2) <- c('intercept', names(df))
    observe2 <- t(as.matrix(df2)) %*% as.matrix(df2)
}
