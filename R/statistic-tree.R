tree <- function(x, depth){
  i <- 0
  while(i<=depth){
    print(funHist(x, varType='numeric', bins=2^depth))
    i <- i+1
  }
} 

#' Accuracy for a differentially private binary tree
#'
#' @param epsilon Numeric differential privacy parameter
#' @param rng Numeric a priori estimate of the variable range
#' @param gran Numeric granularity
#' @param alpha Numeric level of statistical significance, default 0.05
#' @return Accuracy guarantee for the tree given epsilon
#' @export treeGetAccuracy
#' @rdname treeGetAccuracy
 
# treeGetAccuracy <- function(epsilon, rng, gran, alpha=0.05) {
#     universeSize <- diff(rng) / gran + 1
#     accuracy <- (2 * sqrt(2) / epsilon) * sqrt(log(2 / alpha)) * log2(universeSize)^(1.5)
#     return(accuracy)
# }


#' Epsilon for a differentially private binary tree
#'
#' @param accuracy Numeric accuracy needed
#' @param rng Numeric a priori estimate of the variable range
#' @param gran Numeric granularity
#' @param alpha Numeric level of statistical significance, default 0.05
#' @return Epsilon necessary to guarantee the given accuracy
#' @export treeGetParameters
#' @rdname treeGetParameters

# treeGetParameters <- function(accuracy, rng, gran, alpha=0.05) {
#     universeSize <- diff(rng) / gran + 1
#     epsilon <- (2 * sqrt(2) / accuracy) * sqrt(log(2 / alpha)) * log2(universeSize)^(1.5)
#     return(epsilon)
# }


#' Function to truncate negative noisy node counts at zero
#'
#' @param release The differentially private noisy binary tree
#' @return Noisy binary tree truncated at zero

# treePostFormatRelease <- function(release) {
#     release <- round(release)
#     release[release < 0] <- 0
#     return(release)
# }


#' Function to derive CDF from efficient terminal node counts
#'
#' @param release Efficient differentially private binary tree
#' @param rng An a priori estimate of the range of the vector
#'      being represented as a binary tree
#' @param terminalIndex Vector of indices corresponding to the terminal
#'      leaf nodes of the binary tree
#' @return Differentially private estimate of the empirical cumulative
#'      distribution function

treePostCDF <- function(release, rng, terminalIndex) {
    terminal <- release[terminalIndex]
    stepSize <- diff(rng) / length(terminal)
    cdfSteps <- seq(rng[1], rng[2], stepSize)
    cdf <- c(0, cumsum(terminal) / sum(terminal))
    cdf <- data.frame(list('val' = cdfSteps, 'cdf' = cdf))
    return(cdf)
}


#' Function to evaluate the mean using the DP CDF
#'
#' @param cdf Differentially private estimate of the empirical cumulative
#'      distribution function
#' @param rng Numeric a priori estimate of the range
#' @param gran Granularity
#' @return Differentially private estimate of the mean

treePostMean <- function(cdf, rng) {
    ecdf <- cdf$cdf
    pdf <- sapply(2:length(ecdf), function(i) ecdf[i] - ecdf[i - 1])
    p <- c(ecdf[1], pdf) * cdf$val
    return(sum(p))
}


#' Function to evaluate the median using the DP CDF
#'
#' @param cdf Differentially private estimate of the empirical cumulative
#'      distribution function
#' @return Differentially private estimate of the median

treePostMedian <- function(cdf) {
    outMedian <- treePostPercentiles(cdf, 0.5)$value
    return(outMedian)
}


#' Quantile function using the DP CDF
#'
#' @param cdf Differentially private estimate of the empirical cumulative
#'      distribution function
#' @param percentiles Vector of probabilities given to the quantile function
#' @return Differnetially private estimate of the values corresponding to
#'      the provided probabilities

treePostPercentiles <- function(cdf, percentiles) {
    absArgMin <- function(q, cdf) {
        target <- abs(q - cdf$cdf)
        out <- cdf$val[which(target == min(target))]
        return(c(q, mean(out)))
    }
    outValues <- lapply(percentiles, absArgMin, cdf)
    outValues <- data.frame(do.call(rbind, outValues))
    names(outValues) <- c('percentile', 'value')
    return(outValues)
}


#' Function to efficiently estimate noisy node counts
#'
#' @param release The truncated differentially private noisy binary tree
#'      in vector form
#' @param treeData Data frame with binary tree attributes, including depth
#'      and indicators of parent and adjacent nodes. Note that
#'      \code{nrow(treeData) == length(release)}
#' @param n Number of observations
#' @param nNodes Number of nodes in the binary tree, also \code{length(release)}
#' @param variance The variance of the noise used to perturb tree nodes
#' @param terminalIndex Vector of indices corresponding to the terminal
#'      leaf nodes of the binary tree
#' @return Efficient differentially private binary tree

treePostEfficient <- function(release, treeData, n, variance, terminalIndex) {
    nNodes <- length(release)
    sigma <- sqrt(variance)
    invSigmaSq <- 1 / variance
    tree <- cbind(treeData, release)
    names(tree)[ncol(tree)] <- 'noisy'
    tree <- estBottomUp(tree, min(terminalIndex), nNodes, sigma, invSigmaSq)
    tree <- estTopDown(tree, n, nNodes, sigma, invSigmaSq)
    tree <- estEfficiently(tree, n, nNodes, sigma, invSigmaSq)
    return(round(tree$est.efficient))
}


#' Differentially private binary tree
#'
#' @param mechanism Character, the privacy mechanism.
#' @param varType Character, the R variable type. One of \code{'numeric'} or
#'   \code{'integer'}.
#' @param Variable Character, variable name.
#' @param n Integer, number of observations.
#' @param rng Numeric, a priori estimate of the range.
#' @param gran Numeric, the granularity of the variable.
#' @param epsilon Numeric, privacy cost parameter.
#' @param accuracy Numeric, accuracy guarantee given \code{epsilon}.
#' @param imputeRng Numeric, range within which missing values are imputed. If \code{NULL},
#'   the range provided in \code{rng} is used.
#' @param percentiles Numeric, the percentiles to evaluate in post-processing. Optional, 
#'    default \code{NULL}.
#' @param alpha Numeric, the level of statistical significance. Default 0.05.
#'
#' @import methods
#' @export dpTree
#' @exportClass dpTree
#'
#' @include mechanism.R
#' @include mechanism-laplace.R

dpTree <- setRefClass(
    Class = 'dpTree',
    contains = 'mechanismLaplace',
    fields = list(
      globalEps = 'numeric'
    )
    
)

dpTree$methods(
    initialize = function(mechanism, varType, variable, n, rng=NULL, gran, epsilon=NULL,
                          accuracy=NULL, imputeRng=NULL, percentiles=NULL, alpha=0.05, ...) {
        .self$name <- 'Differentially private binary tree'
        .self$mechanism <- checkMechanism(mechanism, "mechanismLaplace")
        .self$varType <- checkVariableType(varType, c('numeric', 'integer', 'logical', 'character'))
        .self$variable <- variable
        .self$n <- checkN(n)
        .self$rng <- checkRange(rng) # CHANGE
        .self$gran <- checkN(gran, emptyOkay=TRUE) #should be positive whole number
        .self$alpha <- checkNumeric(alpha)
        #.self$sens <- 2 * log2(diff(rng) / gran + 1)
        
        checkVariableType(variable, "character")
        
        if (is.null(epsilon)) {
            .self$accuracy <- checkAccuracy(accuracy)
            .self$epsilon <- treeGetParameters(accuracy, rng, gran, alpha)
        } else {
            .self$epsilon <- checkEpsilon(epsilon)
            .self$accuracy <- treeGetAccuracy(epsilon, rng, gran, alpha)
        }
        if (is.null(imputeRng)) {
            .self$imputeRng <- rng
        } else {
            .self$imputeRng <- imputeRng
        }
        .self$percentiles <- percentiles
})

dpTree$methods(
    release = function(data) {
        x <- data[, variable]
        variance <- 2 * sens / epsilon
        universeSize <- floor(diff(rng) / gran + 1)
        depth <- ceiling(log2(universeSize))
        # for (i in 1:depth){
        #   export(mechanism)$evaluate(funHist, x, sens, (function(out) return(out))) #don't bother with post-processing for now
        # }
        #terminalIndex <- seq(2^(depth - 1), 2^depth - 1)
        #.self$result <- export(mechanism)$evaluate(.self$treeFun, x, sens, .self$postProcess, 
        #                                           variance=variance, universeSize=universeSize, 
        #                                           depth=depth, terminalIndex=terminalIndex, self=.self)
})

#dpTree$methods(
    #treeFun = function(x, universeSize, depth) {
        #tree <- binaryTree(x, n, rng, gran, universeSize, depth)
        #.self$treeData <- tree[, which(names(tree) != 'count')]
        #return(tree$count)
        
#})

dpTree$methods(
    postProcess = function(out, ...) {
        # out$variable <- variable
        # out$release <- treePostFormatRelease(out$release)
        # ellipsisVals <- getFuncArgs(list(...), treePostEfficient)
        # out$release <- do.call(treePostEfficient, c(list(release=out$release, treeData=treeData, n=n), ellipsisVals))
        # ellipsisVals <- getFuncArgs(list(...), treePostCDF)
        # out$cdf <- do.call(treePostCDF, c(list(release=out$release, rng=rng), ellipsisVals))
        # out$mean <- treePostMean(out$cdf, rng)
        # out$median <- treePostMedian(out$cdf)
        # out$accuracy <- .self$accuracy
        # out$epsilon <- .self$epsilon
        # if (!is.null(percentiles)) {
        #     out$percentiles <- treePostPercentiles(out$cdf, percentiles)
        # }
        # return(out)
})
