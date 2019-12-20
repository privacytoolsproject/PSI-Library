#' Sample-and-Aggregate mechanism
#'
#' @import methods
#' @export mechanismSampleAndAggregate
#' @exportClass mechanismSampleAndAggregate
#'
#' @include mechanism.R

mechanismSampleAndAggregate <- setRefClass(
    Class = 'mechanismSampleAndAggregate',
    contains = 'mechanism'
)

mechanismSampleAndAggregate$methods(
  evaluate = function(fun, x, sens, postFun, ...) {
    # CC: keeping this in for now because it is in all the other mechanisms, think more about it later
    x <- censorData(x, .self$varType, .self$rng, .self$bins, .self$rngFormat)
    x <- fillMissing(x, .self$varType, imputeRng=.self$rng, categories=.self$imputeBins)

    ### CC: attempt to add convert x into y, a subsampled/aggregated version ###
    # create partition of indices
    index_partition <- create_partition_of_indices(n_obs = length(x), n_subsets = n_subsets)

    # apply function over each subset of the partition
    fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
    trueVals <- list('vector', length = n_subsets)
    for (i in c(1:n_subsets) {
        # create inputVals on subset
        inputVals <- c(list(x = x[index_partition[[i]]]), fun.args)
        # release non-private statistic calculated on subset
        trueVals[[i]] <- do.call(fun, inputVals)
    }

    # trueVal <- do.call(aggregationMechanism, inputVals)
    # scale <- sens / .self$epsilon
    # release <- trueVal + dpNoise(n=length(trueVal), scale=scale, dist='laplace')
    # out <- list('release' = release)
    # out <- postFun(out, ...)
    # return(out)
  }

)