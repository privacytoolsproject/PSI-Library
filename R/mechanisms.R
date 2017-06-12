#' LaPlace mechanism for releasing differentially private queries
#'
#' @param fun A user supplied function, or string naming a function.  Must accept only one argument, named \code{x}.
#' @param x A vector of data to run the function on.
#' @param var.type Data type of vector x
#' @param n The number of samples
#' @param rng An a priori estimate of the range
#' @param sensitivity numeric
#' @param epsilon numeric
#' @param ... Other arguments passed to \code{fun}
#' @return Differentially private release of function \code{fun} on data \code{x}
#' @examples
#' n <- 1000
#' rng <- c(0,1)
#' x <- runif(n, min=min(rng), max=max(rng))
#' sensitivity <- diff(rng) / n
#' mechanism.laplace(dp.mean, x, 'numeric', rng, sensitivity, 0.5, n)

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
    out$release <- out$stat + rlap(mu=0, b=(sensitivity / epsilon), size=length(out$stat))
    out <- out[names(out) != 'stat']

    # post-processing
    if (!is.null(postlist)) {
        out <- postprocess(out, postlist, ...)
    }
    return(out)
}


#' Exponential mechanism

mechanism.exponential <- function(fun, x, var.type, sensitivity, epsilon, k, postlist=NULL, ...) {

    epsilon <- checkepsilon(epsilon)
    x <- censordata(x, var.type, levels=list(...)$bins)

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


#' Gaussian mechanism

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
    noise <- rnorm(length(out$stat), mean=0, sd=(sensitivity * sqrt(2 * log(1.25 / delta)) / epsilon))
    out$release <- out$stat + noise
    out <- out[names(out) != 'stat']

    # post-processing
    if (!is.null(postlist)) {
        out <- postprocess(out, postlist, ...)
    }
    return(out)
}

#' Cycle through available postprocessing functions for a released statistic
#'
#' @param out list containing differentially private released statistic, and mechanism and statistic names
#' @param var.type Data type of vector x
#' @param rng An a priori estimate of the range
#' @param sensitivity numeric
#' @param epsilon numeric
#' @param postlist List with name, function pairs for post-processing statistics
#' @param ... Other arguments passed to \code{fun}
#' @param Original list with released statistic, appended with available postprocessed releases

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
        rng = 'ANY',
        result = 'ANY',
        alpha = 'numeric',
        accuracy = 'numeric',
        bins = 'ANY',
        n.bins = 'ANY'
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
# Laplace mechanism

mechanismLaplace <- setRefClass(
    Class = 'mechanismLaplace',
    contains = 'mechanism'
)

mechanismLaplace$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismLaplace$methods(
    evaluate = function(fun, x, sens, postFun) {
        xc <- censordata(x, .self$var.type, .self$rng)
        field.vals <- .self$getFunArgs(fun)
        true.val <- do.call(fun, c(list(x=x), field.vals))
        scale <- sens / .self$epsilon
        release <- true.val + rlap(b=scale, size=length(true.val))
        z <- qlap((1 - (.self$alpha / 2)), b=scale)
        interval <- c(release - z, release + z)
        out <- list('release' = release, 'interval' = interval)
        out <- postFun(out)
        return(out)
})
