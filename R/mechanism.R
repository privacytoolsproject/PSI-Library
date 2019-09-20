#' Base mechanism class
#'
#' @import methods
#' @export mechanism
#' @exportClass mechanism
#'
#' @field mechanism Name of the mechanism
#' @field name Name of the statistic
#' @field variable Name of the variable
#' @field varType Variable type
#' @field varTypeOrig Variable type at instantiation
#' @field n Number of observations
#' @field epsilon Differential privacy parameter
#' @field delta Differential privacy parameter
#' @field rng A priori estimate of the variable range
#' @field result List with statistical output
#' @field alpha Level of statistical signficance
#' @field accuracy Accuracy guarantee of the estimate
#' @field bins Bins
#' @field nBins Number of bins
#' @field k Number of bins desired for the release
#' @field error Error
#' @field nBoot Number of bootstrap replications
#' @field bootFun Function passed to the bootstrap mechanism
#' @field imputeRng The range from which to impute missing values
#' @field impute Logical, impute categorical types?
#' @field formula R formula for regression models
#' @field columns Vector of column names
#' @field intercept Logical, is the intercept included?
#' @field stability Logical, use stability histogram
#' @field objective Objective function for regression models
#' @field granularity Granularity
#' @field percentiles Percentiles evaluated by binary tree
#' @field treeData Binary tree attributes needed for efficient estimation

mechanism <- setRefClass(
    Class = 'mechanism',
    fields = list(
        mechanism = 'character',
        name = 'character',
        variable = 'character',
        varType = 'character',
        n = 'numeric',
        epsilon = 'numeric',
        delta = 'ANY', # not 'numeric' because it can be NULL if the Stability mechanism is not being used
        rng = 'ANY', # not 'numeric' because it can be NULL (if the Stability mechanism is being used, or the variable is logical)
        result = 'ANY',
        alpha = 'numeric',
        accuracy = 'ANY', # not 'numeric' because can be NULL if epsilon is given
        k = 'numeric',
        error = 'numeric',
        nBoot = 'ANY', # not 'numeric' becuase can be NULL
        bootFun = 'function',
        imputeRng = 'ANY',
        impute = 'logical',
        formula = 'ANY',
        intercept = 'logical', 
        stability = 'logical',
        objective = 'function',
        sens = 'numeric'
))

mechanism$methods(
    getFields = function() {
        f <- names(getRefClass()$fields())
        out <- setNames(vector('list', length(f)), f)
        for (fd in f) {
            out[[fd]] <- .self[[fd]]
        }
        return(out)
})
