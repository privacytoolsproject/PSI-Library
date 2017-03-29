#' Function to evaluate the covariance matrix from a design matrix
#'
#' @param z Numeric design matrix
#' @param n Integer indicating the number of observations
#' @param epsilon Numeric differential privacy parameter epsilon
#' @param trim.thresh Numeric, default 0.95

dp.covariance <- function(z, n, epsilon, trim.thresh=0.95) {

    # features
    features <- names(z)

    # euclidean norm
    z$dist <- apply(z, 1, function(x) { sqrt(sum(x^2)) } )

    # sort increasing by distance
    z <- z[order(z$dist), ]

    # evaluate sensitivity
    threshold <- n * trim.thresh
    sensitivity <- z$dist[threshold + 1] - z$dist[threshold - 1]

    # add Laplace noise to radius
    threshold <- threshold + rlap(b=(sensitivity / epsilon))

    # trim z
    z <- z[1:threshold, features]

    # noisy z'z
    cov.mat <- cov(z)
    noise <- rlap(b=(sensitivity / epsilon), size=(nrow(cov.mat) * ncol(cov.mat)))
    dim(noise) <- dim(cov.mat)
    release <- cov.mat + noise

    out <- list(
        'threshold' = threshold,
        'sensitivity' = sensitivity,
        'noise' = noise,
        'release' = release)
    return(out)
}
