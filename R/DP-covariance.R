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

    if (is.null(formulae)) {
        formulae <- as.list(formulae)
    } else {
        if (inherits(formulae, 'formula')) { formulae <- list(formulae) }
    }

    # pass to mechanism
    postlist <- list('release' = 'formatRelease')
    if (!is.null(formulae)) {
        postlist <- c(postlist, list('linear.regression' = 'postLinearRegression'))
    }
    release <- mechanism.laplace(fun=dp.covariance, x=x, var.type='numeric', n=n,
                                 rng=rng, epsilon=(epsilon / length(sensitivity)), columns=columns,
                                 sensitivity=sensitivity, intercept=intercept, formulae=formulae,
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

    source('mechanisms.R')
    source('utilities.R')

    set.seed(681521)
    eps.list <- seq(0.1, 1.0, length.out=10)
    n.list <- c(500, 1000, 2000, 5000, 10000)

    out.priv.b0 <- matrix(NA, nrow=length(eps.list), ncol=length(n.list))
    out.true.b0 <- matrix(NA, nrow=length(eps.list), ncol=length(n.list))

    out.priv.b1 <- matrix(NA, nrow=length(eps.list), ncol=length(n.list))
    out.true.b1 <- matrix(NA, nrow=length(eps.list), ncol=length(n.list))
    
    out.priv.se <- matrix(NA, nrow=length(eps.list), ncol=length(n.list))

    for (i in 1:length(eps.list)) {
        eps <- eps.list[i]
        for (j in 1:length(n.list)) { 
            n <- n.list[j]
            x <- rnorm(n, sd=2)
            y <- 0.76 + 0.145 * x + rnorm(n)
            df <- as.data.frame(cbind(y, x))
            cols <- names(df)
            rng <- rbind(range(y), range(x))
            dtypes <- 'numeric'
            form1 <- as.formula('y ~ x')

            est.true <- lm(form1, df)$coefficients
            out.true.b0[i, j] <- est.true[1]
            out.true.b1[i, j] <- est.true[2]

            n.sims <- 1000
            est0.priv.out <- vector(mode='numeric', length=n.sims)
            est1.priv.out <- vector(mode='numeric', length=n.sims)
            se.priv.out <- vector(mode='numeric', length=n.sims)
            for (k in 1:n.sims) {
                release <- covariance.release(df, dtypes, n, eps, rng, cols, intercept=TRUE, formulae=form1)
                est.priv <- release$linear.regression[[1]]
                est.priv.b <- as.numeric(est.priv[, 1])
                est.priv.se <- as.numeric(est.priv[, 2])
                est0.priv.out[k] <- est.priv.b[1]
                est1.priv.out[k] <- est.priv.b[2]
                se.priv.out[k] <- !is.na(est.priv.se)[2]
            }
            out.priv.b0[i, j] <- mean(est0.priv.out)
            out.priv.b1[i, j] <- mean(est1.priv.out)
            out.priv.se[i, j] <- sum(se.priv.out) / n.sims
        }
    }

    abs.error.b0 <- abs(out.true.b0 - out.priv.b0)
    abs.error.b1 <- abs(out.true.b1 - out.priv.b1)
    abs.error.b0 <- round(data.frame(abs.error.b0), 4)
    abs.error.b1 <- round(data.frame(abs.error.b1), 4)
    names(abs.error.b0) <- as.character(n.list)
    rownames(abs.error.b0) <- as.character(eps.list)
    names(abs.error.b1) <- as.character(n.list)
    rownames(abs.error.b1) <- as.character(eps.list)
    out.priv.se <- data.frame(out.priv.se)
    names(out.priv.se) <- as.character(n.list)
    rownames(out.priv.se) <- as.character(eps.list)

}
