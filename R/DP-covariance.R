#' Function to evaluate the covariance matrix from input matrix and specify parameters for post-processing
#'
#' @param x Numeric data frame
#' @param n Integer indicating the number of observations
#' @param rng Numeric 2-tuple, lower and upper bounds for standardized data
#' @param epsilon Numeric differential privacy parameter epsilon
#' @param columns Character vector indicating columns in \code{x}, if NULL then use all columns
#' @param trim.thresh Numeric, default 0.95

dp.covariance <- function(x, n, rng, epsilon, columns, intercept, trim, trim.thresh) {

    # subset and optionally append an intercept
    data <- x[, columns]
    if (intercept) {
        data <- cbind(1, data)
    }

    # trim the extreme values
    if (trim) {
        dist <- order(apply(data, 1, function(x) { max(abs(x)) } ))
        n.trim <- as.integer(trim.thresh * n)
        data <- data[1:n.trim, ]
    }

    # the true statistic, return only the lower triangle
    covariance <- cov(data)
    lower <- lower.tri(covariance, diag=TRUE)
    covariance <- covariance[lower]

    # specify the output
    return(list('stat' = covariance,
                'n' = n,
                'rng' = rng,
                'epsilon' = epsilon,
                'trim.thresh' = trim.thresh))
}



covariance.release <- function(x, var.type, n, epsilon, rng, columns, delta=2e-16, intercept=FALSE, trim=TRUE, trim.thresh=0.95) {

    if (intercept) {
        columns <- c('intercept', columns)
        rng <- rbind(c(1, 1), rng)
    }
    sensitivity <- 2 * length(columns)
    postlist <- NULL

    # pass to mechanism
    release <- mechanism.laplace(fun=dp.covariance, x=x, var.type='numeric', n=n,
                                 rng=rng, epsilon=(n * epsilon), columns=columns,
                                 sensitivity=sensitivity, intercept=intercept, trim=trim,
                                 trim.thresh=trim.thresh, postlist=postlist)

    # convert to symmetric matrix
    out.dim <- ifelse(intercept, length(columns) + 1, length(columns))
    out.matrix <- matrix(0, nrow=out.dim, ncol=out.dim)
    out.matrix[lower.tri(out.matrix, diag=TRUE)] <- release$release
    out.matrix[upper.tri(out.matrix, diag=FALSE)] <- t(out.matrix)[upper.tri(out.matrix, diag=FALSE)]
    release$release <- data.frame(out.matrix)
    rownames(release$release) <- names(release$release) <- columns
    return(release)
}


if (interactive()) {

    set.seed(2)
    n <- 10000
    x1 <- rnorm(n, mean=73, sd=17)
    x2 <- rpois(n, lambda=4)
    x3 <- as.integer(rnorm(n, mean=25, sd=4))
    y <- 0.145 * x1 - 0.565 * x2 + 0.013 * x3 + rnorm(n, sd=1.5)
    df <- as.data.frame(cbind(y, x1, x2, x3))
    cols <- names(df)

    rng <- rbind(range(y), range(x1), range(x2), range(x3))
    dtypes <- 'numeric'
    eps <- 0.5
    release <- covariance.release(df, dtypes, n, eps, rng, cols)
}
