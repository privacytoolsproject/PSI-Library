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

bootstrapReplication <- function(x, n, sensitivity, epsilon, fun) {
    partition <- rmultinom(n=1, size=n, prob=rep(1 / n, n))
    maxAppearances <- max(partition)
    probs <- sapply(1:maxAppearances, dbinom, size=n, prob=(1 / n))
    statPartitions <- vector('list', maxAppearances)
    for (i in 1:maxAppearances) {
        iVariance <- (i * probs[i] * (sensitivity^2)) / (2 * epsilon)
        iStat <- fun(x[partition == i])
        iNoise <- dpNoise(n=length(iStat), scale=sqrt(iVariance), dist='gaussian')
        statPartitions[[i]] <- i * iStat + iNoise
    }
    statOut <- do.call(rbind, statPartitions)
    return(apply(statOut, 2, sum))
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
    bootStatEval = function(xi,...) {
        funArgs <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
        inputVals = c(list(x=x), funArgs)
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
    evaluate = function(fun, x, sens) {
        x <- censorData(x, .self$varType, .self$rng)
        x <- fillMissing(x, .self$varType, .self$imputeRng[0], .self$imputeRng[1])
        epsilonPart <- epsilon / .self$nBoot
        release <- replicate(.self$nBoot, bootstrapReplication(x, n, sens, epsilonPart, fun=.self$bootStatEval))
        stdError <- .self$bootSE(release, .self$nBoot, sens)
        out <- list('release' = release, 'stdError' = stdError)
        return(out)
})
