#' Function to evaluate the covariance matrix from input matrix and specify parameters for post-processing
#'
#' @param x Numeric data frame
#' @param n Integer indicating the number of observations
#' @param rng Numeric 2-tuple, lower and upper bounds for standardized data
#' @param epsilon Numeric differential privacy parameter epsilon
#' @param columns Character vector indicating columns in \code{x}, if NULL then use all columns
#' @param trim.thresh Numeric, default 0.95

dp.covariance <- function(x, n, rng, epsilon, columns, intercept, trim, trim.thresh) {

    # subset, standardize, and censor
    data <- x[, columns]
    data <- scale(data)
    data <- censor(data, rng)

    # optionally append an intercept
    if (intercept) {
        data <- cbind(1, data)
    }

    # trim the extreme values
    if (trim) {
        dist <- order(apply(data, 1, function(x) { max(abs(x)) } ))
        n.trim <- as.integer(trim.thresh * n)
        data <- data[1:n.trim, ]
    }

    # means and st devs to unscale, use this instead of info in scale() to account for trimming??
    means <- apply(data, 2, mean)
    std.devs <- apply(data, 2, sd)

    # the true statistic, return only the lower triangle
    covariance <- cov(data)
    lower <- lower.tri(covariance, diag=TRUE)
    covariance <- covariance[lower]

    # specify the output
    return(list('stat' = covariance,
                'n' = n,
                'rng' = rng,
                'epsilon' = epsilon,
                'trim.thresh' = trim.thresh,
                'means' = means,
                'std.devs' = std.devs))
}


# need to allow for an intercept

covariance.release <- function(x, var.type, n, epsilon, rng, columns, delta=2e-16, intercept=FALSE, trim=TRUE, trim.thresh=0.95) {

    if (intercept) { columns <- c('intercept', columns) }

    rng <- checkrange(rng)
    sensitivity <- sqrt((rng[2] - rng[1])^2)  # L2 bound
    postlist <- NULL

    release <- mechanism.gaussian(fun=dp.covariance, x=x, var.type='numeric', n=n,
                                  rng=rng, epsilon=epsilon, delta=delta, columns=columns,
                                  sensitivity=sensitivity, intercept=intercept, trim=trim,
                                  trim.thresh=trim.thresh, postlist=postlist)

    release$corr <- release$release

    # unscale release and convert to symmetric matrix
    out.dim <- ifelse(intercept, length(columns) + 1, length(columns))
    out.matrix <- matrix(0, nrow=out.dim, ncol=out.dim)
    out.matrix[lower.tri(out.matrix, diag=TRUE)] <- release$release
    for (row in 1:nrow(out.matrix)) {
        for (col in row:ncol(out.matrix)) {
            out.matrix[row, col] <- out.matrix[row, col] * release$std.devs[row] * release$std.devs[col]

        }
    }
    out.matrix[upper.tri(out.matrix, diag=FALSE)] <- t(out.matrix)[upper.tri(out.matrix, diag=FALSE)]
    release$release <- data.frame(out.matrix)
    rownames(release$release) <- names(release$release) <- columns

    return(release)
}


#' Gaussian mechanism

mechanism.gaussian <- function(fun, x, var.type, rng, sensitivity, epsilon, delta, postlist=NULL, ...) {

    epsilon <- checkepsilon(epsilon)

    mechanism.args <- c(as.list(environment()), list(...))
    out <- do.call(fun, getFuncArgs(mechanism.args, fun))
    noise <- rnorm(length(out$stat), mean=0, sd=(sensitivity * sqrt(2 * log(1.25 / delta) * epsilon)))
    print(noise)
    out$release <- out$stat + noise
    out <- out[names(out) != 'stat']

    # post-processing
    if (!is.null(postlist)) {
        out <- postprocess(out, postlist, ...)
    }
    return(out)
}


# if rng is expressed in standard deviations, then trim the orginal data frame by column?
# does censordata fn work here?
censor <- function(z, rng) {
    rng <- sort(rng)
    z[z < rng[1]] <- rng[1]
    z[z > rng[2]] <- rng[2]
    return(z)
}


if (interactive()) {

    set.seed(2)

    source('utilities.R')

    n <- 10000
    x1 <- rnorm(n, mean=73, sd=17)
    x2 <- rpois(n, lambda=4)
    x3 <- as.integer(rnorm(n, mean=25, sd=4))
    y <- 34.5 + 0.145 * x1 - 0.565 * x2 + 0.013 * x3 + rnorm(n, sd=1.5)
    df <- as.data.frame(cbind(y, x1, x2, x3))
    cols <- names(df)

    rng <- c(-3.0, 3.0)
    dtypes <- 'numeric'
    eps <- 0.5

    check <- covariance.release(df, dtypes, n, eps, rng, cols)

    #######################

    stop()

    #######################

    library(mvtnorm)

    # fake data
    N <- 5000
    d <- 10
    vcov <- matrix(0, nrow=d, ncol=d)
    diag(vcov) <- 1
    vcov[1:3, 1:3] <- c(1.0, 0.3, 0.15, 0.3, 1.0, 0.1, 0.15, 0.1, 1.0)

    # draw some data & bound
    df <- as.data.frame(rmvnorm(n=N, mean=rep(0, d), sigma=vcov))
    intercept <- rep(1, N)
    lower.bound <- -4.0
    upper.bound <- 4.0
    df[df < lower.bound] <- lower.bound
    df[df > upper.bound] <- upper.bound
    df <- cbind(1, df)
    names(df) <- c('intercept', 'y', paste0('x', 1:(d-1)))
    features <- names(df)

    # trim details
    trim.thresh <- 0.95
    n.trim <- N * trim.thresh

    # covariance routine
    df$dist <- apply(df, 1, function(x) { max(abs(x)) } )
    df <- df[order(df$dist, decreasing=TRUE), features]
    df.trim <- as.matrix(df[1:n.trim, features])

    # inner product
    inner <- function(x) {
        return(t(x) %*% x)
    }

    # covariance
    vc <- inner(as.matrix(df))
    vc.trim <- inner(df.trim)

    # largest possible sensitivity
    l2.bound <- sqrt((upper.bound - lower.bound)^2)
    delta <- 2e-16
    epsilon <- 0.1

    # symmetric noise matrix
    the.matrix <- matrix(NA, nrow=nrow(vc), ncol=ncol(vc))
    n.draws <- (nrow(the.matrix) + 1) * nrow(the.matrix) / 2
    draws <- rnorm(n.draws, mean=0, sd=(l2.bound^2 * sqrt(2 * log(1.25 / delta) / epsilon)))
    the.matrix[upper.tri(the.matrix, diag=TRUE)] <- draws
    the.matrix[lower.tri(the.matrix, diag=FALSE)] <- t(the.matrix)[lower.tri(the.matrix, diag=FALSE)]

    # jackknife the vc matrix for local sensitivity
    vc.lower <- vc.upper <- inner(df.trim[1, , drop=FALSE])
    for (row in 2:nrow(df.trim)) {
        row.inner <- inner(df.trim[row, , drop=FALSE])
        below.flag <- row.inner < vc.lower
        above.flag <- row.inner > vc.upper
        vc.lower[below.flag] <- row.inner[below.flag]
        vc.upper[above.flag] <- row.inner[above.flag]
    }
    jack.vc.diff <- vc.upper - vc.lower

    # symmetric noise
    trim.matrix <- matrix(NA, nrow=nrow(vc.trim), ncol=ncol(vc.trim))
    for(row in 1:nrow(vc.trim)) {
        for (col in row:ncol(vc.trim)) {
            noise <- rnorm(1, mean=0, sd=jack.vc.diff[row, col] * sqrt(2 * log(1.25 / delta) / epsilon))
            trim.matrix[row, col] <- trim.matrix[col, row] <- noise
        }
    }

    priv <- vc + the.matrix
    priv.trim <- vc.trim + trim.matrix
}
