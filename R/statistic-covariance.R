#' Function to get the sensitivity of the lower triangle of the covariance matrix
#'
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data frame.
#' @param rng A numeric list of 2-tuples of the lower and upper bounds for
#'    each of the variables.
#' @param intercept A logical vector of length one indicating whether an
#'    intercept should be added prior to evaluating the inner product x'x.
#' @return The sensitivity of the data frame for which the covariance matrix
#'   is being calculated.
#'
#  A complete derivation of the sensitivity of the covariance may be found in
#' /extra_docs/sensitivities/covariance_sensitivity.pdf
#' @example Should output equal to rep(2/10000, 3)
#' range.sex <- range(PUMS5extract10000['sex'])
#' range.married <- range(PUMS5extract10000['married'])
#' range <- list(range.sex, range.married)
#' covarianceSensitivity(10000, range, FALSE)
covarianceSensitivity <- function(n, rng, intercept) {
    diffs <- sapply(rng, diff)
    if (intercept) { diffs <- c(0, diffs) }
    sensitivity <- c()
    const <- 2/n
    for (i in 1:length(diffs)) {
        for (j in i:length(diffs)) {
            s <- const * diffs[i] * diffs[j]
            sensitivity <- append(sensitivity, s)
        }
    }
    return(sensitivity)
}

#' Lower triangle of covariance matrix
#'
#' This function is called by \code{mechanismLaplace$evaluate}, which is called within the
#' differentially private covariance release function \code{dpCovariance$release}. It
#' produces the true covariance matrix of input \code{x} that is then perturbed within
#' \code{mechanismLaplace$evaluate}.
#'
#' Since the Laplace mechanism as instantiated within \code{mechanismLaplace$evaluate} expects
#' an output of the true function which can have a 1-dimensional vector of noise added to it,
#' the covariance function here outputs a flattened version of the lower triangle of the
#' covariance matrix, rather than a matrix. The lower triangle is sufficient since covariance
#' matrices are symmetric.
#'
#' The flattened output corresponds to the covariance matrix by proceeding
#' proceeding top-to-bottom down each column of the lower triangular matrix (including the diagonal).
#'
#' The traditional covariance matrix is then reconstructed from the noisy version of this output
#' as a post-processing step in \code{covarianceFormatRelease}.
#'
#' @param x Input data frame that covariance matrix will be calculated with.
#' @param intercept Logical, indicates if intercept column should be appended to x.
#' @return The lower triangle of the covariance matrix of x, flattened to a 1-dimensional array.
#'
#' @seealso dpCovariance$release
#' @seealso mechanismLaplace$evaluate

covar <- function(x, intercept) {
    if (intercept) { x <- cbind(1, x) }
    covariance <- cov(x)
    lower <- lower.tri(covariance, diag=TRUE)
    covariance <- covariance[lower]
    return(covariance)
}


#' Differentially private covariance matrix
#'
#' @import methods
#' @export dpCovariance
#' @exportClass dpCovariance
#'
#' @include mechanism.R
#' @include mechanism-laplace.R
#'
#' @field epsilon       Vector of epsilon values where the ith epsilon will be used for the ith covariance calculation in flattened
#'  lower triangle of covariance matrix. (The flattened output corresponds to the covariance matrix by proceeding
#'  proceeding top-to-bottom down each column of the lower triangular matrix including the diagonal.)
#' @field accuracy      Single accuracy value that will be used to compute each epsilon for each individual covariance calculation
#'   in the covariance matrix.
#' @field globalEps     Global epsilon to be split between all of the covariance calculations
#' @field epsilonDist   Vector of percentages (valued 0 to 1) that describes how global epsilon \code{globalEps} should be
#'   split for each covariance calculation.
#' @field  accuracyVals     Vector of accuracy values where the ith accuracy will be used for the ith covariance calculation in
#'  flattened lower triangle of covariance matrix. (The flattened output corresponds to the covariance matrix by proceeding
#'  proceeding top-to-bottom down each column of the lower triangular matrix including the diagonal.)
#' @field formula       R formula for regression models

dpCovariance <- setRefClass(
    Class = 'dpCovariance',
    contains = 'mechanismLaplace',
    fields = list(
        epsilonDist = 'numeric',
        globalEps = 'numeric',
        accuracyVals = 'numeric',
        formula = 'ANY'
    )
)

dpCovariance$methods(
  initialize = function(mechanism, varType, n, columns, rng, epsilon=NULL, globalEps=NULL, epsilonDist= NULL,
                        accuracy=NULL, accuracyVals=NULL, imputeRng=NULL, intercept=FALSE, formula=NULL,
                        alpha=0.05) {

    .self$name <- 'Differentially private covariance matrix'
    .self$mechanism <- checkMechanism(mechanism, "mechanismLaplace")
    .self$varType <- checkVariableType(varType, c('integer', 'double', 'numeric', 'logical'))  #NEED TO CHANGE TO ALLOW MULTIPLE VARTYPES


    .self$formula <- formula
    .self$intercept <- intercept
    .self$alpha <- alpha


    .self$n <- checkN(n)
    .self$rngFormat <- 'list'
    .self$rng <- checkRange(rng, .self$varType, .self$rngFormat, expectedLength=length(columns))
    .self$sens <- covarianceSensitivity(n, rng, intercept)


    checkVariableType(class(formula), c("formula", "character"), emptyOkay=TRUE)
    checkVariableType(typeof(intercept), "logical")
    checkVariableType(typeof(columns), "character")


    if (is.null(imputeRng)) { # NEED TO ALLOW FOR IMPUTING ON ONLY SOME OF THE RANGES
      .self$imputeRng <- rng
    } else {
      .self$imputeRng <- imputeRng
    }

    if (.self$intercept) {
      .self$columns <- c('intercept', columns)
    } else {
      .self$columns <- columns
    }

    # Distribute epsilon across all covariances that will be calculated
    outputLength <- lowerTriangleSize(.self$columns)
    # Option 1: Enter vector of epsilon values to be used for each covariance calculation in matrix.
    if (!is.null(epsilon)){
        .self$epsilon <- checkEpsilon(epsilon, expectedLength=outputLength)
        .self$globalEps <- sum(.self$epsilon)
    }
    # Option 2: Enter global epsilon value and vector of percentages specifying how to split global
    # epsilon between covariance calculations.
    else if (!is.null(globalEps) && !is.null(epsilonDist)){
        .self$globalEps <- checkEpsilon(globalEps)
        .self$epsilonDist <- checkEpsilonDist(epsilonDist, outputLength)
        .self$epsilon <- distributeEpsilon(.self$globalEps, epsilonDist=.self$epsilonDist)
        .self$accuracyVals <- laplaceGetAccuracy(.self$sens, .self$epsilon, .self$alpha)
    }
    # Option 3: Only enter global epsilon, and have it be split evenly between covariance calculations.
    else if (!is.null(globalEps)){
        .self$globalEps <- checkEpsilon(globalEps)
        .self$epsilon <- distributeEpsilon(.self$globalEps, nCalcs=outputLength)
        .self$accuracyVals <- laplaceGetAccuracy(.self$sens, .self$epsilon, .self$alpha)
    }
    # Option 4: Enter an accuracy value instead of an epsilon, and calculate individual epsilons with this accuracy.
    else if (!is.null(accuracy)){
        .self$accuracy = checkAccuracy(accuracy)
        .self$epsilon = laplaceGetEpsilon(.self$sens, .self$accuracy, .self$alpha)
        .self$globalEps = sum(.self$epsilon)
    }
    # Option 5: Enter vector of accuracy values, and calculate ith epsilon value from ith accuracy value
    else if (!is.null(accuracyVals)){
        .self$accuracyVals = checkAccuracyVals(accuracyVals, outputLength)
        .self$epsilon = laplaceGetEpsilon(.self$sens, .self$accuracyVals, .self$alpha)
        .self$globalEps = sum(.self$epsilon)
    }
  })

dpCovariance$methods(
  release = function(data) {
    x <- data[columns];
    .self$result <- export(mechanism)$evaluate(fun=covar, x=x, sens=sens, postFun=.self$postProcess,
                                               formula=formula, columns=columns, intercept=intercept)
  })

dpCovariance$methods(
  postProcess = function(out, columns, formula, intercept) {

    out$release <- covarianceFormatRelease(out$release, columns)
    out$variable <- columns
    if (!is.null(formula)) {
      out$linearRegression <- covariancePostLinearRegression(out$release, n, intercept, formula)
    }
    return(out)
  })
