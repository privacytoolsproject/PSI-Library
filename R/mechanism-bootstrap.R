#' Bootstrap replication for a function
#'
#' @param x Vector
#' @param n Number of observations
#' @param sensitivity Sensitivity of the function
#' @param epsilon Numeric differential privacy parameter
#' @param fun Function to evaluate
#' @return Value of the function applied to one bootstrap sample
#' @import stats
#' @export

bootstrap.replication <- function(x, n, sensitivity, epsilon, fun) {
    partition <- rmultinom(n=1, size=n, prob=rep(1 / n, n))
    max.appearances <- max(partition)
    probs <- sapply(1:max.appearances, dbinom, size=n, prob=(1 / n))
    stat.partitions <- vector('list', max.appearances)
    for (i in 1:max.appearances) {
        variance.i <- (i * probs[i] * (sensitivity^2)) / (2 * epsilon)
        stat.i <- fun(x[partition == i])
        noise.i <- dpNoise(n=length(stat.i), scale=sqrt(variance.i), dist='gaussian')
        stat.partitions[[i]] <- i * stat.i + noise.i
    }
    stat.out <- do.call(rbind, stat.partitions)
    return(apply(stat.out, 2, sum))
}


#' Bootstrap mechanism
#'
#' @import methods
#' @export mechanismBootstrap
#' @exportClass mechanismBootstrap
#'
#' @include mechanism.R

mechanismBootstrap <- setRefClass(
    Class = 'mechanismBootstrap',
    contains = 'mechanism'
)

mechanismBootstrap$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismBootstrap$methods(
    bootStatEval = function(xi) {
        field.vals <- .self$getFunArgs(boot.fun)
        stat <- do.call(boot.fun, c(list(xi=xi), field.vals))
        return(stat)
})

mechanismBootstrap$methods(
    bootSE = function(release, n.boot, sens) {
        se <- sd(release)
        c.alpha <- qchisq(0.01, df=(n.boot - 1))
        conservative <- sqrt(max(c(se^2 - (c.alpha * sens^2 * n.boot) / (2 * epsilon * (n.boot - 1)), 0)))
        naive <- sqrt(max(c(se^2 - (sens^2 * n.boot) / (2 * epsilon), 0)))
        return(list('sd' = se,
                    'conservative' = conservative,
                    'naive' = naive))
})

mechanismBootstrap$methods(
    evaluate = function(fun, x, sens, postFun) {
        x <- censordata(x, .self$var.type, .self$rng)
        x <- fillMissing(x, .self$var.type, .self$impute.rng[0], .self$impute.rng[1])
        epsilon.part <- epsilon / .self$n.boot
        release <- replicate(.self$n.boot, bootstrap.replication(x, n, sens, epsilon.part, fun=.self$bootStatEval))
        std.error <- .self$bootSE(release, .self$n.boot, sens)
        out <- list('release' = release, 'std.error' = std.error)
        out <- postFun(out)
        return(out)
})
