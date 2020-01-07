#' Sample and Aggregate "Statistic"
#'
#' @param innerFun Character, the function calculated within each subset of the partition.
#' @param aggregationFun Character, the function calculated over the set of results from
#'   each subset of the partition.
#' @param mechanism Character, the privacy mechanism. For \code{dpSampleAndAggregate}, one
#'   of \code{c('mechanismLaplace')}.
#' @param numSubsets, Numeric, the number of subsets desired in the partition
#' @param varType Character, the R variable type. One of \code{c('numeric',
#'   'integer', 'logical')}.
#' @param variable Character, variable name.
#' @param n Integer, number of observations
#' @param innerFunSens Numeric, sensitivity of inner function.
#' @param aggregationFunSens Numeric, sensitivity of aggregation function.
#' @param subsetSize Numeric, minimum size of a subset of the partition.
#' @param rng Numeric, a priori estimate of the range
#' @param epsilon Numeric, privacy cost parameter
#' @param accuracy Numeric, accuracy guarantee given \code{epsilon}
#' @param imputeRng Numeric, range within which missing values are imputed. If \code{NULL},
#'   the range provided in \code{rng} is used.
#' @param alpha Numeric, the level of statistical significance. Default 0.05.
#'
#' @import methods
#' @export dpSampleAndAggregate
#' @exportClass dpSampleAndAggregate
#'
#'
#' NOTE: will need to include all mechanisms and statistics used in the "mechanism" and "aggregationFun" fields, respectively
#' NOTE: when a new DP statistic is added as an aggregation function, it will need to have an option to
#'       accept sens as an argument, rather than setting it as a function of rng and n. Also, every DP statistic should follow
#'       the naming convention dp<x>. The argument value for aggregationFun should correspond to the <x> -- so for the dpMean
#'       statistic, the argument value is 'Mean'.
#' @include mechanism.R
#' @include mechanism-laplace.R
#' @include statistic-mean.R

dpSampleAndAggregate <- setRefClass(
    Class = 'dpSampleAndAggregate',
    contains = c('mechanismLaplace',
                 'dpMean'), # NOTE: "contains" vector will need to contain every mechanism and
                           # aggregationFun used. The aggregationFuns all correspond to DP statistics
    fields = list(innerFun = 'character', aggregationFun = 'character', numSubsets = 'numeric',
                  innerFunSens = 'numeric', aggregationFunSens = 'numeric', subsetSize = 'numeric')
)

dpSampleAndAggregate$methods(
    initialize = function(innerFun, aggregationFun, mechanism, numSubsets,
                          varType, variable, n, innerFunSens=NULL,
                          aggregationFunSens=NULL, subsetSize=NULL, rng=NULL, epsilon=NULL,
                          accuracy=NULL, imputeRng=NULL, alpha=0.05, ...) {
        ### establish acceptable inner functions, aggregation functions, and privacy mechanisms ###
        acceptableInnerFuns <- c('mean', 'median', 'var')
        acceptableAggregationFuns <- c('Mean')
        acceptableMechanisms <- c('mechanismLaplace')

        ### perform parameter setup/checking ###
        .self$name <- 'Differentially private sample and aggregate'
        .self$innerFun <- checkFunction(innerFun, acceptableInnerFuns, type = 'Inner')
        .self$aggregationFun <- checkFunction(aggregationFun, acceptableAggregationFuns, type = 'Aggregation')
        .self$mechanism <- checkMechanism(mechanism, acceptableMechanisms)
        .self$numSubsets <- checkNumSubsets(numSubsets)
        checkVariableType(typeof(variable), c('character'))
        .self$varType <- checkVariableType(varType, c('numeric', 'integer', 'logical'))
        .self$variable <- variable
        .self$n <- checkN(n)
        .self$subsetSize <- floor(.self$n / .self$numSubsets) #NOTE: this assumes that the data partitioning
                                                              #      is done such that no two subsets vary
                                                              #      in size by more than one element
        if (.self$subsetSize < 1) {
            stop("numSubsets cannot be larger than n")
        }
        .self$alpha <- checkNumeric(alpha)
        .self$rngFormat <- 'vector'
        .self$rng <- checkRange(rng, .self$varType, .self$rngFormat, expectedLength=1)

        ### set/get sensitivities for inner and aggregation functions ###
        # set inner function sensitivity
        .self$innerFunSens <- get_inner_fun_sens(.self, innerFun)

        # set aggregation function sensitivity
        .self$aggregationFunSens <- get_aggregation_fun_sens(.self, aggregationFun)

        ### finish parameter setup ###
        if (is.null(epsilon)) {
            .self$accuracy <- checkAccuracy(accuracy, expectedLength=1)
            # NOTE: need to update this for mechanisms other than Laplace if/when they are added
            .self$epsilon <- laplaceGetEpsilon(.self$aggregationFunSens, .self$accuracy, .self$alpha)
        } else {
            checkEpsilon(epsilon)
            .self$epsilon <- checkEpsilon(epsilon, expectedLength=1)
            # NOTE: need to update this for mechanisms other than Laplace if/when they are added
            .self$accuracy <- laplaceGetAccuracy(.self$aggregationFunSens, .self$epsilon, .self$alpha)
        }

        if (is.null(imputeRng)) {
            .self$imputeRng <- .self$rng
        } else {
            .self$imputeRng <- checkImputationRange(imputationRange=imputeRng, rng=.self$rng, varType=.self$varType)
        }
})


dpSampleAndAggregate$methods(
    #' Differentially private sample and aggregate release
    #'
    #' @name dpSampleAndAggregateRelease
    #' @param data Dataframe with a column named .self$variable, where
    #'  that column has data of type .self$varType and which is bounded by
    #'  .self$rng.
    #'
    #' Assigns to .self$result a dataframe that describes the differentially private
    #' statistic, calculated by some mechanism as defined in .self$mechanism, of that
    #' column of the dataframe.
    #'
    #' Note that the actual differentially private release is calculated in a call to the
    #' differentially private mechanism .self$mechanism's \code{evaluate} function within
    #' the \code{dpSampleAndAggregate$release} function.
    release = function(data, ...) {
        # establish and preprocess data
        x <- data[, variable]
        x <- censorData(x, .self$varType, .self$rng, .self$bins, .self$rngFormat)
        x <- fillMissing(x, .self$varType, imputeRng=.self$rng, categories=.self$imputeBins)

        # create partition of indices
        index_partition <- create_partition_of_indices(n_obs = length(x), n_subsets = .self$numSubsets)

        # apply function over each subset of the partition
        innerFun.args <- getFuncArgs(innerFun, inputList=list(...), inputObject=.self)
        subset_vals <- list('vector', length = numSubsets)
        for (i in c(1:numSubsets)) {
            # create inputVals on subset
            input_vals <- c(list(x = x[index_partition[[i]]]), innerFun.args)
            # release non-private statistic calculated on subset
            subset_vals[[1]][i] <- do.call(innerFun, input_vals)
        }
        # create data frame of function values over subsets
        subset_vals_df <- data.frame(subset_vals[[1]])
        colnames(subset_vals_df) <- .self$variable

        # ensure that data are correct data type
        for (col in colnames(subset_vals_df)) {
            subset_vals_df[, col] <- eval(parse(text = paste0('as.', .self$varType)))(as.character(subset_vals_df[, col]))
        }

        # apply aggregation statistic and privacy mechanism
        dp_stat <- eval(parse(text = paste0('dp', .self$aggregationFun)))$new(mechanism=.self$mechanism, varType=.self$varType, variable = variable,
                                                                            epsilon = .self$epsilon, n = .self$numSubsets, rng = .self$rng, sens = .self$aggregationFunSens,
                                                                            ...)
        # release aggregation statistic and save result to dpSampleAndAggregate object
        dp_stat$release(subset_vals_df)
        .self$result <- dp_stat$result
    }
)