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

bootstrapReplication <- function(x, n, sensitivity, epsilon, fun, inputObject, ...) {
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

# 2: treat it as a partition with a statistic of value 0 and keep it in the calculation, adding noise and adding it to the final calculation

# bootstrapReplication <- function(x, n, sensitivity, epsilon, fun, inputObject, ...) {
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
    bootStatEval = function(xi, fun,...) {
        funArgs <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
        inputVals = c(list(x=xi), funArgs)
        stat <- do.call(bootFun, inputVals)
        return(stat)
})

mechanismBootstrap$methods(
    bootSE = function(release, nBoot, sens) {
        se <- sd(release)
        cAlpha <- qchisq(0.01, df=(nBoot - 1))
        conservative <- sqrt(max(c(se^2 - (cAlpha * sens^2 * nBoot) / (2 * epsilon * (nBoot - 1)), 0)))
        naive <- sqrt(max(c(se^2 - (sens^2 * nBoot) / (2 * epsilon), 0)))
        return(list('sd' = se,
                    'conservative' = conservative,
                    'naive' = naive))
})

mechanismBootstrap$methods(
    evaluate = function(fun, x, sens, postFun, ...) {
        x <- censorData(x, .self$varType, .self$rng, rngFormat=.self$rngFormat)
        x <- fillMissing(x, .self$varType, .self$imputeRng[0], .self$imputeRng[1])
        epsilonPart <- epsilon / .self$nBoot
        release <- replicate(.self$nBoot, bootstrapReplication(x, .self$n, sens, epsilonPart, fun, .self))
        stdError <- .self$bootSE(release, .self$nBoot, sens)
        out <- list('release' = release, 'stdError' = stdError)
        out <- postFun(out)
        return(out)
})
