#' LaPlace Noise Mechanism for Differential Privacy
#' 
#' Mechanism to add noise from a LaPlace distribution for releasing 
#'    differentially private queries.
#'
#' @param fun A user supplied function or string vector of length one naming a 
#'    function.
#' @param x A vector of data to run \code{fun} on.
#' @param var.type Data type of vector \code{x}.
#' @param rng An a priori estimate of the range of \code{x}.
#' @param sensitivity The difference of \code{rng} divided by \code{n}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param postlist A list with names, function pairs for post-processing 
#'    statistics.
#' @param ... Other arguments passed to \code{fun}.
#' 
#' @return Differentially private release of function \code{fun} on 
#'    data \code{x}.
#' @examples
#' 
#' n <- 1000
#' rng <- c(0,1)
#' x <- runif(n, min=min(rng), max=max(rng))
#' sensitivity <- diff(rng) / n
#' mechanism.laplace(dp.mean, x, 'numeric', rng, sensitivity, 0.5, n=n)
#' @seealso \code{\link{mechanism.gaussian}}
#' @rdname mechanism.laplace
#' @export
mechanism.laplace <- function(fun, x, var.type, rng, sensitivity, epsilon, postlist=NULL, ...) {

    # checks & transformations
    epsilon <- checkepsilon(epsilon)
    if (var.type == 'logical') { x <- make_logical(x) }
    if (var.type %in% c('numeric', 'integer', 'logical')) {
        rng <- checkrange(rng)
        x <- censordata(x, var.type, rng=rng)
    } else {
        x <- censordata(x, var.type, levels=list(...)$bins)
    }

    # evaluate the noisy statistic
    mechanism.args <- c(as.list(environment()), list(...))
    out <- do.call(fun, getFuncArgs(mechanism.args, fun))
    out$release <- out$stat + dpNoise(n=length(out$stat), scale=(sensitivity / epsilon), dist='laplace')
    out <- out[names(out) != 'stat']

    # post-processing
    if (!is.null(postlist)) {
        out <- postprocess(out, postlist, ...)
    }
    return(out)
}

#' Exponential Noise Mechanism for Differential Privacy
#' 
#' Mechanism to add exponential noise for releasing differentially private 
#'    queries.
#' 
#' @param fun A user supplied function or string vector of length one naming a 
#'    function.
#' @param x A vector of data to run \code{fun} on.
#' @param var.type Data type of vector \code{x}.
#' @param sensitivity The sensitivity of \code{x}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param k The number of desired releases.
#' @param postlist A list with names, function pairs for post-processing 
#'    statistics.
#' @param ... Other arguments passed to \code{fun}.
#' 
#' @return Differentially private release of function \code{fun} on 
#'    data \code{x}.
#' @examples
#' 
#' n <- 1000
#' epsilon <- 0.5
#' delta <- 1e-7
#' observed.levels <- bins <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
#' probs <- c(0.40, 0.25, 0.15, 0.10, 0.04, 0.03, 0.02, 0.01)
#' x <- sample(observed.levels, size=n, prob=probs, replace=TRUE)
#' mechanism.exponential(fun=dp.heavyhitters, x=x, var.type = 'character', sensitivity = 2,
#'                      epsilon = epsilon, k=3, bins=bins, n=n, delta=delta)
#' @rdname mechanism.exponential
#' @export
mechanism.exponential <- function(fun, x, var.type, sensitivity, epsilon, k, postlist=NULL, ...) {

    epsilon <- checkepsilon(epsilon)
    x <- censordata(x, var.type, rng, levels=list(...)$bins)

    mechanism.args <- c(as.list(environment()), list(...))
    out <- do.call(fun, getFuncArgs(mechanism.args, fun))
    quality <- out$stat - max(out$stat)
    probs <- ifelse(out$stat == 0, 0, exp((epsilon * quality) / (2 * sensitivity)))
    out$release <- sample(names(out$stat), k, prob=probs)
    out <- out[names(out) != 'stat']

    if (!is.null(postlist)) {
        out <- postprocess(out, postlist, ...)
    }
    return(out)
}

#' Gaussian Noise Mechanism for Differential Privacy
#' 
#' Mechanism to add noise from a Gaussian distribution for releasing 
#'    differentially private queries.
#'
#' @param fun A user supplied function or string vector of length one naming a 
#'    function.
#' @param x A vector of data to run \code{fun} on.
#' @param var.type Data type of vector \code{x}.
#' @param rng An a priori estimate of the range of \code{x}.
#' @param sensitivity The difference of \code{rng} divided by \code{n}.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param delta The probability of an arbitrary leakage of information from 
#'    \code{x}. Should be of length one and should be a very small value. 
#'    Default to 10^-6.
#' @param postlist A list with names, function pairs for post-processing 
#'    statistics.
#' @param ... Other arguments passed to \code{fun}.
#' 
#' @return Differentially private release of function \code{fun} on 
#'    data \code{x}.
#' @examples
#' 
#' n <- 1000
#' rng <- c(0, 1)
#' x <- runif(n, min=min(rng), max=max(rng))
#' sensitivity <- diff(rng) / n
#' mechanism.gaussian(fun=dp.mean, x=x, var.type='numeric', rng=rng, sensitivity=sensitivity, 
#'    epsilon=0.5, delta=0.000001, n=n)
#' @seealso \code{\link{mechanism.laplace}}
#' @rdname mechanism.gaussian
#' @export
mechanism.gaussian <- function(fun, x, var.type, rng, sensitivity, epsilon, delta, postlist=NULL, ...) {

    # checks
    epsilon <- checkepsilon(epsilon)
    if (var.type %in% c('numeric', 'integer', 'logical')) {
        rng <- checkrange(rng)
        x <- censordata(x, var.type, rng=rng)
    } else {
        x <- censordata(x, var.type, levels=list(...)$bins)
    }

    # evaluate the noisy statistic
    mechanism.args <- c(as.list(environment()), list(...))
    out <- do.call(fun, getFuncArgs(mechanism.args, fun))
    scale <- sensitivity * sqrt(2 * log(1.25 / delta)) / epsilon
    out$release <- out$stat + dpNoise(n=length(out$stat), scale=scale, dist='gaussian')
    out <- out[names(out) != 'stat']

    # post-processing
    if (!is.null(postlist)) {
        out <- postprocess(out, postlist, ...)
    }
    return(out)
}

#' Cycle through available postprocessing functions for a released statistic
#'
#' @param out A list containing differentially private released statistic, and 
#'    mechanism and statistic names.
#' @param postlist A list with names, function pairs for post-processing statistics.
#' @param ... Other arguments.
#' 
#' @return Postprocessed statistics.
#' @rdname postprocess
postprocess <- function(out, postlist, ...) {
    available.attrs <- c(out, list(...))
    for (process in names(postlist)) {
        get.name <- paste0(out$name, ".", postlist[[process]])
        if (exists(get.name, mode='function')) {
            available.attrs[[process]] <- out[[process]] <- do.call(get.name, getFuncArgs(available.attrs, get.name))
        } else {
            out[[process]] <- 'Function not provided'
        }
    }
    return(out)
}

# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# base mechanism class

mechanism <- setRefClass(
    Class = 'mechanism',
    fields = list(
        mechanism = 'character',
        name = 'character',
        var.type = 'character',
        n = 'numeric',
        epsilon = 'numeric',
        delta = 'numeric',
        rng = 'ANY',
        result = 'ANY',
        alpha = 'numeric',
        accuracy = 'numeric',
        bins = 'ANY',
        n.bins = 'ANY',
        error = 'numeric'
    )
)

mechanism$methods(
    getFields = function() {
        f <- names(getRefClass()$fields())
        out <- setNames(vector('list', length(f)), f)
        for (fd in f) {
            out[[fd]] <- .self[[fd]]
        }
        return(out)
})

mechanism$methods(
    getFunArgs = function(fun) {
        f <- .self$getFields()
        spec <- list()
        for (arg in names(f)) {
            if (arg %in% names(formals(fun))) {
                spec[[arg]] <- f[[arg]]
            }
        }
        return(spec)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
#' Laplace mechanism
#' 
#' @export mechanismLaplace
#' @exportClass mechanismLaplace

mechanismLaplace <- setRefClass(
    Class = 'mechanismLaplace',
    contains = 'mechanism'
)

mechanismLaplace$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismLaplace$methods(
    evaluate = function(fun, x, sens, postFun, ...) {
        xc <- censordata(x, .self$var.type, .self$rng)
        field.vals <- .self$getFunArgs(fun)
        true.val <- do.call(fun, c(list(x=x), field.vals))
        scale <- sens / .self$epsilon
        release <- true.val + dpNoise(n=length(out$stat), scale=scale, dist='laplace')
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
# Exponential mechanism

mechanismExponential <- setRefClass(
    Class = 'mechanismExponential',
    contains = 'mechanism'
)

mechanismExponential$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismExponential$methods(
    evaluate = function(fun, x, sens, postFun, ...) {
        xc <- censordata(x, .self$var.type, levels=.self$bins)
        field.vals <- .self$getFunArgs(fun)
        true.val <- do.call(fun, c(list(x=x), field.vals))
        quality <- true.val - max(true.val)
        probs <- ifelse(true.value == 0, 0, exp((.self$epsilon * quality) / (2 * sens)))
        release <- sample(names(true.value), .self$k, prob=probs)
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
#' Gaussian mechanism

mechanismGaussian <- setRefClass(
    Class = 'mechanismGaussian',
    contains = 'mechanism'
)

mechanismGaussian$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismGaussian$methods(
    evaluate = function(fun, x, sens, postFun) {
        xc <- censordata(x, .self$var.type, .self$rng)
        field.vals <- .self$getFunArgs(fun)
        true.val <- do.call(fun, c(list(x=x), field.vals))
        scale <- sens * sqrt(2 * log(1.25 / .self$delta)) / .self$epsilon
        release <- true.val + dpNoise(n=length(out$stat), scale=scale, dist='gaussian')
        out <- list('release' = release)
        out <- postFun(out)
        return(out)
})

# --------------------------------------------------------- #
# --------------------------------------------------------- #
