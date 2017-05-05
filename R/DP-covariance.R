#' Function to evaluate the covariance matrix from input matrix and specify parameters for post-processing
#'
#' @param x Numeric data frame
#' @param n Integer indicating the number of observations
#' @param rng Numeric matrix of 2-tuples with the lower and upper bounds for each of P variables, dimensions Px2
#' @param epsilon Numeric differential privacy parameter epsilon
#' @param columns Character vector indicating columns in \code{x}
#' @param intercept Logical indicating whether an intercept should be added prior to evaluating the inner product x'x

dp.covariance <- function(x, n, rng, epsilon, columns, intercept) {

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

    # specify the output
    return(list('stat' = covariance,
                'n' = n,
                'rng' = rng,
                'epsilon' = epsilon,
                'columns' = columns,
                'intercept' = intercept))
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

covariance.release <- function(x, var.type, n, epsilon, rng, columns, intercept=FALSE) {

    # get the vector of sensitivities from the ranges
    diffs <- apply(rng, 1, diff)
    if (intercept) { diffs <- c(1, diffs) }
    sensitivity <- c()
    for (i in 1:length(diffs)) {
        for (j in i:length(diffs)) {
            sensitivity <- ((n - 1) / n) * diffs[i] * diffs[j]
        }
    }

    # pass to mechanism
    postlist <- NULL
    release <- mechanism.laplace(fun=dp.covariance, x=x, var.type='numeric', n=n,
                                 rng=rng, epsilon=(epsilon / length(sensitivity)), columns=columns,
                                 sensitivity=sensitivity, intercept=intercept, postlist=postlist)

    # convert to symmetric matrix
    out.dim <- ifelse(intercept, length(columns) + 1, length(columns))
    out.matrix <- matrix(0, nrow=out.dim, ncol=out.dim)
    out.matrix[lower.tri(out.matrix, diag=TRUE)] <- release$release
    out.matrix[upper.tri(out.matrix, diag=FALSE)] <- t(out.matrix)[upper.tri(out.matrix, diag=FALSE)]
    release$release <- data.frame(out.matrix)
    rownames(release$release) <- names(release$release) <- release$columns
    return(release)
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
    eps <- 0.5

    #release <- covariance.release(df, dtypes, n, eps, rng, cols)
    #observe <- t(as.matrix(df)) %*% as.matrix(df)

    release2 <- covariance.release(df, dtypes, n, eps, rng, cols, intercept=TRUE)
    df2 <- cbind(1, df)
    names(df2) <- c('intercept', names(df))
    observe2 <- t(as.matrix(df2)) %*% as.matrix(df2)

    source('mechanisms.R')
    source('utilities.R')


    extract.locs <- function(formula, data, intercept) {
        t <- terms(formula, data=data)
        x.locs <- which(names(data) %in% attr(t, 'term.labels'))
        x.labels <- names(data)[x.locs]
        y.loc <- attr(t, 'response')
        if (intercept) {
            intercept.loc <- which(names(data) == 'intercept')
            x.locs <- c(intercept.loc, x.locs)
            x.labels <- append(x.labels, 'intercept', after=(intercept.loc - 1))
            if (intercept.loc <= y.loc) { y.loc <- y.loc + 1 }
        }
        return(list('y.loc' = y.loc,
                    'x.locs' = x.locs,
                    'x.labels' = x.labels))
    }

    regress <- function(form, release, n, intercept) {
        xy.locs <- extract.locs(form, release, intercept)
        x.locs <- xy.locs$x.locs
        y.loc <- xy.locs$y.loc
        loc.vec <- rep(TRUE, (length(x.locs) + 1))
        loc.vec[y.loc] <- FALSE
        coefs <- amsweep((as.matrix(release) / n), loc.vec)[y.loc, x.locs]
        coefs <- data.frame(coefs)
        rownames(coefs) <- xy.locs$x.labels
        names(coefs) <- 'Estimate'
        return(coefs)
    }

    form <- as.formula('y ~ x1 + x2 + x3')
    estimates <- regress(form, release2$release, release2$n, release2$intercept)

    x <- c(1, 3, 4, 5)
    y <- 2
    loc <- rep(TRUE, 5)
    loc[y] <- FALSE
    swept <- amsweep((as.matrix(release2$release) / n), loc)[y, x]

    linreg <- lm(form, data=df)
}
