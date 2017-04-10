#' Function to evaluate the covariance matrix from a design matrix
#'
#' @param z Numeric design matrix
#' @param n Integer indicating the number of observations
#' @param epsilon Numeric differential privacy parameter epsilon
#' @param trim.thresh Numeric, default 0.95

dp.covariance <- function(z, n, epsilon, trim.thresh=0.95) {
    return(list('z' = z,
                'n' = n,
                'epsilon' = epsilon,
                'trim.thresh' = trim.thresh))
}


if (interactive()) {

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

    the.matrix <- matrix(NA, nrow=nrow(vc), ncol=ncol(vc))
    for (row in 1:nrow(vc)) {
        n.draws <- nrow(vc) - row + 1
        draws <- rnorm(n.draws, mean=0, sd=(l2.bound^2 * sqrt(2 * log(1.25 / delta) / epsilon)))
        the.matrix[row, row:nrow(vc)] <- the.matrix[row:nrow(vc), row] <- draws
    }

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
