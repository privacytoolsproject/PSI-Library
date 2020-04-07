#' UnbiasedPrivacy "Statistic"
#'
#'
#'
#' @param statistic Function that calculates quantity of interest (should change this to a character name?)
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
    fields = list(statistic = 'ANY',
                  name = 'character',
                  B = 'numeric',
                  n = 'numeric', 
                  P = 'numeric', 
                  lambda = 'numeric', 
                  lambda_var = 'numeric', 
                  delta ='numeric', 
                  epsilon = 'numeric', 
                  epsilon_alpha = 'numeric', 
                  censoring_cutoff = 'numeric', 
                  bias_cutoff = 'numeric')
)

dpUnbiasedPrivacy$methods(
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
    #' TODO: make the passing of form & coef (for statistic coefFn) more dynamic
    #' TODO: remove parallelize
    #' TODO: remove hardcoded ref to coefFn
    release = function(data, form, coefVal...) {
        # call from simulation code:
        #sim <- algorithmUDP(data = dat, statistic = coefFn, B = pr$R, n = pr$b, P = pr$P, lambda = l, lambda_var = 0.025, 
        #           delta = 0.01, epsilon = pr$e, epsilon_alpha = pr$e_alpha, 
        #            parallelize = F, censoring_cutoff = 0.9,
       #            bias_cutoff = 0.1, form = form, coef = coef)
      # browser()
       algorithmUDP(data = data, statistic = .self$statistic, B = .self$B, n= .self$n, P = .self$P, lambda =  .self$lambda, lambda_var = .self$lambda_var,
                    delta = .self$delta, epsilon =  .self$epsilon, epsilon_alpha =  .self$epsilon_alpha, 
                    parallelize = F, censoring_cutoff =  .self$censoring_cutoff, 
                    bias_cutoff =  .self$bias_cutoff, form = form, coefVal = coefVal)
    }
)