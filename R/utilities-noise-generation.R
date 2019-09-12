#' Get accuracy for Laplace statistics
#'
#' Function to find the accuracy guarantee of a statistic release at a given epsilon
#' value.
#'
#' @param sensitivity the sensitivity of the statistic
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param alpha A numeric vector specifying the statistical significance level.
#'
#' @return Accuracy guarantee for statistic release given epsilon.

laplace.getAccuracy <- function(sensitivity, epsilon, alpha=0.05) {
    accuracy <- log(1 / alpha) * (sensitivity / epsilon)
    return(accuracy)
}

#' Get accuracy for Snapping statistics
#'
#' Function to find the accuracy guarantee of a statistic release at a given epsilon
#' value.
#'
#' @param sensitivity the sensitivity of the statistic
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param alpha A numeric vector specifying the statistical significance level.
#'
#' @return Accuracy guarantee for statistic release given epsilon.

snapping.getAccuracy <- function(sensitivity, epsilon, alpha=0.05) {
    accuracy <- (1 + log(1 / alpha)) * (sensitivity / epsilon)
    return(accuracy)
}

#' Get epsilon for Laplace statistics
#'
#' Function to find the epsilon value necessary to meet a desired level of
#' accuracy for a statistic release.
#'
#' @param sensitivity the sensitivity of the statistic
#' @param accuracy A numeric vector representing the accuracy needed to
#'    guarantee (percent).
#' @param alpha A numeric vector specifying the statistical significance level.
#'
#' @return The scalar epsilon necessary to guarantee the needed accuracy.

laplace.getEpsilon <- function(sensitivity, accuracy, alpha=0.05) {
    epsilon <- log(1 / alpha) * (sensitivity / accuracy)
    return(epsilon)
}

#' Get epsilon for Snapping statistics
#'
#' Function to find the epsilon value necessary to meet a desired level of
#' accuracy for a statistic release.
#'
#' @param sensitivity the sensitivity of the statistic
#' @param accuracy A numeric vector representing the accuracy needed to
#'    guarantee (percent).
#' @param alpha A numeric vector specifying the statistical significance level.
#'
#' @return The scalar epsilon necessary to guarantee the needed accuracy.

snapping.getEpsilon <- function(sensitivity, accuracy, alpha=0.05) {
    epsilon <- (1 + log(1 / alpha)) * (sensitivity / accuracy)
    return(epsilon)
}

#' Differentially Private Uniform Draw
#'
#' Draw cryptographically secure random variates from a uniform distribution.
#'
#' @param n An integer giving number of variates needed.
#' @param seed An integer indicating a seed for R's PNRG, defaults to \code{NULL}.
#'
#' @return Random numeric vector of length \code{n} containing values between
#'    zero and one.
#'
#' Draws secure random variates from the uniform distribution through \code{openssl}.
#' If a seed is provided, the \code{runif} function is used to draw the random variates.
#' @examples
#'
#' uniform_secure <- dpUnif(n=1000)
#' uniform_repeatable <- dpUnif(n=1, seed=75436)
#' @seealso \code{\link{dpNoise}}
#' @rdname dpUnif
#' @export
dpUnif <- function(n, seed=NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
        return(runif(n))
    }
    return(openssl::rand_num(n))
}


#' Differentially Private Noise Generator
#'
#' Compile noise from a cryptographically secure random variates to achieve
#'    differentially private statistics.
#'
#' @param n An integer giving number of variates needed.
#' @param scale Numeric, the scale for the distribution.
#' @param dist A character specifying the distribution from which to draw the
#'    noise.
#' @param shape An integer giving the shape parameter for the gamma
#'    distribution. Default to \code{NULL}.
#' @param seed An integer indicating a seed for R's PNRG, defaults
#'    to \code{NULL}.
#'
#' @return Cryptographically secure noise vector or matrix.
#' @examples
#'
#' laplace_noise <- dpNoise(n=1000, scale=1, dist='laplace')
#' gaussian_noise <- dpNoise(n=1000, scale=1, dist='gaussian')
#' laplace_noise_repeatable <- dpNoise(n=1, scale=1, dist='laplace', seed=96845)
#' @seealso \code{\link{dpUnif}}
#' @rdname dpNoise
#' @export
dpNoise <- function(n, scale, dist, shape=NULL, seed=NULL) {
    u <- dpUnif(n, seed)
    if (dist == 'laplace') {
        return(qlap(u, b=scale))
    } else if (dist == 'gaussian') {
        return(qnorm(u, sd=scale))
    } else if (dist == 'gamma') {
        return(qgamma(u, scale=scale, shape=shape))
    } else {
        stop(sprintf('Distribution "%s" not understood', dist))
    }
}

#' Noise Generator for Snapping Mechanism
#'
#' Generate noise for use in Snapping Mechanism.
#'
#' @param true_val True value of the parameter.
#' @param n An integer giving number of variates needed.
#' @param sens Sensitivity of the function for which we are releasing a private estimate.
#' @param epsilon Desired level of differential privacy
#' @param B Bound such that the true value should lie within [-B, B]
#'
#' noise <- snappingNoise(true_val = 50, sens = 5e-4, epsilon = 1e-3, B = 200)
#' @rdname snappingNoise
#' @export
snappingNoise <- function(true_val, n, sens, epsilon, B) {
    # source snapping mechanism file
    reticulate::source_python(system.file('python', 'cc_snap.py', package = 'PSIlence'))

    # initialize noise vector
    noise <- c()

    for (i in 1:n) {
        # intialize snapping mechanism object
        # TODO: should sens, epsilon, B be able to vary by value of true_val?
        snapping_mech <- Snapping_Mechanism(mechanism_input = true_val[[i]],
                                            sensitivity = sens,
                                            epsilon = epsilon,
                                            B = B)

        # calculate and return noise (as well as epsilon and accuracy)
        noise[[i]] <- snapping_mech$get_snapped_noise()
    }
    return(noise)
}


#' Random draw from Laplace distribution
#'
#' @param mu numeric, center of the distribution
#' @param b numeric, spread
#' @param size integer, number of draws
#'
#' @return Random draws from Laplace distribution
#' @examples
#'
#' rlap(size=1000)
#' @export
rlap = function(mu=0, b=1, size=1) {
    p <- runif(size) - 0.5
    draws <- mu - b * sgn(p) * log(1 - 2 * abs(p))
    return(draws)
}


#' Probability density for Laplace distribution
#'
#' @param x numeric, value
#' @param mu numeric, center of the distribution
#' @param b numeric, spread
#'
#' @return Density for elements of x
#' @examples
#'
#' x <- seq(-3, 3, length.out=61)
#' dlap(x)
#' @export
dlap <- function(x, mu=0, b=1) {
    dens <- 0.5 * b * exp(-1 * abs(x - mu) / b)
    return(dens)
}


#' LaPlace Cumulative Distribution Function
#'
#' Determines the probability a draw from a LaPlace distribution is less than
#'    or equal to the specified value.
#'
#' @param x Numeric, the value(s) at which the user wants to know the CDF height.
#' @param mu Numeric, the center of the LaPlace distribution, defaults to 0.
#' @param b Numeric, the spread of the LaPlace distribution, defaults to 1.
#'
#' @return Probability the LaPlace draw is less than or equal to \code{x}.
#' @examples
#'
#' x <- 0
#' plap(x)
#' @rdname plap
#' @export
plap <- function(x, mu=0, b=1) {
    cdf <- 0.5 + 0.5 * sgn(x - mu) * (1 - exp(-1 * (abs(x - mu) / b)))
    return(cdf)
}


#' Quantile function for Laplace distribution
#'
#' @param p Numeric, vector of probabilities
#' @param mu numeric, center of the distribution
#' @param b numeric, spread
#' @return Quantile function
#' @examples
#' probs <- c(0.05, 0.50, 0.95)
#' qlap(probs)
#' @export
qlap <- function(p, mu=0, b=1) {
    q <- ifelse(p < 0.5, mu + b * log(2 * p), mu - b * log(2 - 2 * p))
    return(q)
}


#' Sign function
#'
#' Function to determine what the sign of the passed values should be.
#'
#' @param x numeric, value or vector or values
#' @return The sign of passed values
#' @examples
#' sgn(rnorm(10))
#' @export
sgn <- function(x) {
    return(ifelse(x < 0, -1, 1))
}
