#' UnbiasedPrivacy "Statistic"
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
#'BEGIN UNBIASED PRIVACY
#' @param statistic Function that calculates quantity of interest (should change this to a character name)
#' @param B Integer, Number of bootstraps to run via BLB algorithm
#' @param n Integer, Split size
#' @param P Integer of partitions 
#' @param lambda Numeric, Bounding parameter for the QOI
#' @param lambda_var Numeric, Bounding parameter for the variance
#' @param delta Numeric, Privacy parameter
#' @param epsilon Numeric, Privacy budget for the QOI
#' @param epsilon_alpha Numeric, Privacy budget for estimating alpha^{dp}
#' @param censoring_cutoff Numeric, Maximum amount of censoring to allow
#' @param bias_cutoff Maximum amount of censoring to allow withou doing bias correction
#' @param parallelize Whether to parallelize the BLB calculations (EK -  removed)
#' @param ... Parameters necessary for \code{statistic}
#' 
#'
#'
#' @import methods
#' @export dpUnbiasedPrivacy
#' @exportClass dpUnbiasedPrivacy
#'
#'
#'
#algorithmUDP <- function(data, statistic, B, n, P, lambda, lambda_var, delta, epsilon = 0.1, epsilon_alpha = 0.1, 
#                         censoring_cutoff = 0.6, bias_cutoff = 0.1, parallelize = F, ...) {
        
dpUnbiasedPrivacy <- setRefClass(
    Class = 'dpUnbiasedPrivacy',
   # contains = c('mechanismLaplace',
   #              'dpMean'),
   # NOTE: "contains" vector will need to contain every mechanism and
                           # aggregationFun used. The aggregationFuns all correspond to DP statistics
    fields = list(statistic, B, n, P, lambda, lambda_var, delta, epsilon = 0.1, epsilon_alpha = 0.1, 
                                           censoring_cutoff = 0.6, bias_cutoff = 0.1)
)

dpSampleAndAggregate$methods(
    initialize = function(statistic, B, n, P, lambda, lambda_var, delta, epsilon = 0.1, epsilon_alpha = 0.1, 
                          censoring_cutoff = 0.6, bias_cutoff = 0.1, ...) {
        ### establish acceptable inner functions, aggregation functions, and privacy mechanisms ###
        acceptableQOIs <- c('mean')  #TODO - validate this when parameter is a character
       
        ### perform parameter setup/checking ###
        .self$name <- 'Unbiased Privacy Algorithm'
        .self$statistic <- statistic
        .self$B <- B
        .self$n <- n
        .self$P <- P
        .self$lambda <- lambda
        .self$lambda_var <- lambda_var
        .self$delta <- delta
        .self$epsilon <- epsilon
        .self$epsilon_alpha <- epsilon_alpha
        .self$censoring_cutoff <- censoring_cutoff
        .self$bias_cutoff <- bias_cutoff
        
        #TODO - add validation functions for all these parameters
       
})

#algorithmUDP <- function(data, statistic, B, n, P, lambda, lambda_var, delta, epsilon = 0.1, epsilon_alpha = 0.1, 
#                         censoring_cutoff = 0.6, bias_cutoff = 0.1, parallelize = F, ...) {

dpUnbiasedPrivacy$methods(
    #' Differentially private sample and aggregate release
    #'
    #' @name dpUnbiasedPrivacyRelease
    #' @param data Dataframe with a column named .self$variable, where
    #'  that column has data of type .self$varType and which is bounded by
    #'  .self$rng.
    #'
    #' #TODO update description of release method
    #' Assigns to .self$result a dataframe that describes the differentially private
    #' statistic, calculated by some mechanism as defined in .self$mechanism, of that
    #' column of the dataframe.
    #'
    #' Note that the actual differentially private release is calculated in a call to the
    #' differentially private mechanism .self$mechanism's \code{evaluate} function within
    #' the \code{dpSampleAndAggregate$release} function.
    release = function(data, ...) {
       algorithmUDP(data, self$statistic, self$B, self$n, self$P, self$lambda, self$lambda_var,
                    self$delta, self$epsilon, self$epsilon_alpha, self$censoring_cutoff, self$bias_cutoff)
    }
)