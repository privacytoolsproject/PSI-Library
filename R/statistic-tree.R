#' Calculate bins of each leaf in the tree by row
#'
#' @param rng Total range of data to be binned into tree in form (min, max) A tuple of numerics of length 2.
#' @param depth Depth of the tree, considering the root to have depth 0.
#' @param n Number of data points to be binned.
#'
#' Note: This can be calculated entirely publically. The fact that 'n' is input is unnecessary,
#' but is a residual affect of the fact that the helper function called can also be used if the user
#' specifies their desired granularity of the bins instead of supplying a range, which may be something 
#' we want to eventually implement here.
#' 
#' @return A list of the bins used in the tree, from the first level of the tree to the leaves.
treeBins <- function(rng, depth, n){
  binsByLevel <- list()
  i <- 1
  while(i <= depth){
    nBins <- 2^i
    bins <- determineNumericIntegerBins(rng, n, nBins, NULL) #NULL is passed as granularity since that is unnecessary here
    binsByLevel <- append(binsByLevel, list(bins))
    i <- i+1
  }
  return(binsByLevel)
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
      globalEps = 'numeric',
      depth = 'numeric',
      binsByLevel = 'list',
      bins = 'numeric'
    )
    
)

dpTree$methods(
    initialize = function(varType, variable, n, depth, rng=NULL, globalEps=NULL,
                          accuracy=NULL, imputeRng=NULL, alpha=0.05, ...) {
        .self$name <- 'Differentially private binary tree'
        .self$mechanism <- "mechanismLaplace"
        .self$varType <- checkVariableType(varType, c('numeric', 'integer'))
        .self$variable <- variable
        .self$n <- checkN(n)
        .self$depth <- checkN(depth)
        .self$rng <- checkRange(rng, .self$varType, 'vector')
        .self$rngFormat <- 'vector'
        .self$alpha <- checkNumeric(alpha)
        
        checkVariableType(typeof(variable), "character")
        .self$sens <- 2
        
        # Option 1: Specify global epsilon value
        if (!is.null(globalEps)){
          .self$globalEps <- checkEpsilon(globalEps)
          .self$epsilon <- .self$globalEps/depth
          .self$accuracy <- laplaceGetAccuracy(.self$sens, .self$epsilon)
        }
        # Option 2: Specify epsilon value for each row of tree
        else if (!is.null(epsilon)){
          .self$epsilon <- checkEpsilon(epsilon)
          .self$globalEps <- checkEpsilon(epsilon*.self$depth)
          .self$accuracy <- laplaceGetAccuracy(.self$sens, .self$epsilon)
        }
        # Option 3: Specify an accuracy value
        else if (!is.null(accuracy)){
          .self$accuracy <- checkAccuracy(accuracy)
          .self$epsilon <- laplaceGetEpsilon(.self$sens, .self$accuracy, .self$alpha)
          .self$globalEps <- .self$epsilon * .self$depth
        }
        
        .self$binsByLevel <- treeBins(rng, depth, n)
        
        if (is.null(imputeRng)) {
            .self$imputeRng <- rng
        } else {
            .self$imputeRng <- imputeRng
        }
})

dpTree$methods(
    release = function(data) {
        x <- data[, variable]
        counts <- list(n) #n is public so the root of tree need not be noisy. Double nested just to match fact that later elements are also lists.
        names(counts[[1]]) <- paste("[",toString(.self$rng),"]") # adding bin range to the root node
        #counts[.self$rng]
        i <- 1
        while(i <= .self$depth){
          .self$bins <- .self$binsByLevel[[i]]   #Bins of ith row (note this is publically computable)
          noisyCount <- export(mechanism)$evaluate(funHist, x, .self$sens, (function(out) return(out)))
          # In evaluate, identity function is passed as postProcess function since we want to postprocess
          # on all of the noisy counts together.
          counts <- append(counts, list(noisyCount$release))
          i <- i+1
        }
        #Note: postprocessing is called here instead of in the evaluate function
        out <- list('release' = counts)
        .self$result <- .self$postProcess(out)
})

dpTree$methods(
    postProcess = function(out) {
      
        out$epsilon <- .self$epsilon # epsilon used for each of the node calculations
        out$globalEps <- .self$globalEps # global epsilon used in total
        out$accuracy <- .self$accuracy
        out$variable <- variable
        out$bins <- .self$binsByLevel
        
        # release an optimal version of the tree that leverages the multiple levels of information to produce higher accuracy counts
        out$optimalPostProcess <- optimalPostProcess(out$release, .self$epsilon)
        # release a cdf of the tree with highest granularity possible (equivalent to granularity of the bins in the leaf level of the tree)
        out$postCDF <- treePostCDF(out$optimalPostProcess$optimalTree, out$bins)
        # release the median from the cdf, or a best estimate if the 50th percentile is not an option from the cdf calculation 
        out$postMedian <- cdfMedian(out$postCDF)
        # release an estimate of the mean from the tree, using the leaf bins and counts. 
        out$postMean <- treeMean(out$optimalPostProcess$optimalTree, out$bins)
        
        return(out)
})
