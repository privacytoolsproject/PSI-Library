#' Bootstrap replication for a function
#'
#' @param x Vector
#' @param n Number of observations
#' @param sensitivity Sensitivity of the function
#' @param epsilon Numeric differential privacy parameter
#' @param fun Function to evaluate
#' @param inputObject the Bootstrap mechanism object on which the input function will be evaluated
#' @return Value of the function applied to one bootstrap sample
#' @import stats
#' @export

# There are 2 options for handling empty partitions:

# 1: skip it entirely, and say the total number of partitions is just the number of partitions that are not empty

bootstrap.replication <- function(x, n, sensitivity, epsilon, fun, inputObject, ...) {
    partition <- rmultinom(n=1, size=n, prob=rep(1 / n, n))
    # make a sorted vector of the partitions of the data
    # because it is not guaranteed that every partition from 1:max.appearances will have a value in it
    # so we need to loop through only the partitions that have data
    validPartitions <- sort(unique(partition[,1]))
    # we do not want the 0 partition, so we remove it from the list
    validPartitions <- validPartitions[2:length(validPartitions)]
    # print the unique values of the partition, to track which entries may result in NaN
    print(validPartitions)
    probs <- sapply(1:length(validPartitions), dbinom, size=n, prob=(1 / n))
    stat.partitions <- vector('list', length(validPartitions))
    for (i in 1:length(validPartitions)) {
        currentPartition <- validPartitions[i]
        variance.currentPartition <- (currentPartition * probs[i] * (sensitivity^2)) / (2 * epsilon)
        stat.currentPartition <- inputObject$bootStatEval(x[partition == currentPartition], fun, ...)
        noise.currentPartition <- dpNoise(n=length(stat.currentPartition), scale=sqrt(variance.currentPartition), dist='gaussian')
        stat.partitions[[i]] <- currentPartition * stat.currentPartition + noise.currentPartition
    }
    stat.out <- do.call(rbind, stat.partitions)
    # return(apply(stat.out, 2, sum))
    returnedBootstrappedResult <- apply(stat.out, 2, sum)
    return(returnedBootstrappedResult)
}

# 2: treat it as a partition with a mean of 0 and keep it in the calculation, adding noise and adding it to the final calculation

# bootstrap.replication <- function(x, n, sensitivity, epsilon, fun, inputObject, ...) {
#     partition <- rmultinom(n=1, size=n, prob=rep(1 / n, n))
#     # make a sorted vector of the partitions of the data
#     # because it is not guaranteed that every partition from 1:max.appearances will have a value in it
#     validPartitions <- validPartitions <- sort(unique(partition[,1]))
#     # print the unique values of the partition, to track which entries may result in NaN
#     print(validPartitions)
#     max.appearances <- max(partition)
#     probs <- sapply(1:max.appearances, dbinom, size=n, prob=(1 / n))
#     stat.partitions <- vector('list', max.appearances)
#     for (i in 1:max.appearances) {
#         variance.i <- (i * probs[i] * (sensitivity^2)) / (2 * epsilon)
#         if (i %in% validPartitions) {
#             stat.i <- inputObject$bootStatEval(x[partition == currentPartition], fun, ...)
#             noise.i <- dpNoise(n=length(stat.i), scale=sqrt(variance.i), dist='gaussian')
#             stat.partitions[[i]] <- i * stat.i + noise.i
#         } else {
#             stat.i <- 0
#             noise.i <- dpNoise(n=length(stat.i), scale=sqrt(variance.i), dist='gaussian')
#             stat.partitions[[i]] <- i * stat.i + noise.i
#         }
#     }
#     stat.out <- do.call(rbind, stat.partitions)
#     return(apply(stat.out, 2, mean))
# }


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
    bootStatEval = function(xi, fun, ...) {
        fun.args <- getFuncArgs(boot.fun, inputList=list(...), inputObject=.self)
        input.vals = c(list(x=xi), fun.args)
        stat <- do.call(boot.fun, input.vals)
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
    evaluate = function(fun, x, sens, postFun, ...) {
        x <- censordata(x, .self$var.type, .self$rng)
        x <- fillMissing(x, .self$var.type, .self$impute.rng[0], .self$impute.rng[1])
        epsilon.part <- epsilon / .self$n.boot
        release <- replicate(.self$n.boot, bootstrap.replication(x, n, sens, epsilon.part, fun=fun, inputObject = .self, ...))
        std.error <- .self$bootSE(release, .self$n.boot, sens)
        out <- list('release' = release, 'std.error' = std.error)
        out <- postFun(out)
        return(out)
})
