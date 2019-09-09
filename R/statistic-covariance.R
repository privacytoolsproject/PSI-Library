#' Release additional model coefficients from DP covariance matrix
#' 
#' Function to extract regression coefficients using the differentially private covariance matrix 
#'    via the sweep operator. This is the function to obtain coefficients for additional models
#'    after the \code{covariance.release()} function has already been called. 
#' 
#' @param formula Formula, regression formula used on data
#' @param release Numeric, private release of covariance matrix
#' @param n Integer, indicating number of observations
#' @export

coefficientRelease <- function(formula, release, n) {
  intercept <- ifelse('intercept' %in% names(release), TRUE, FALSE)
  coefficients <- linearReg(formula, release, n, intercept)
  release <- list(name='Linear regression', 
                  n=n, 
                  formula=formula, 
                  coefficients=coefficients)
  return(release)
}


#' Function to get the sensitivity of the covariance matrix
#'
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data frame.
#' @param rng A numeric matrix of 2-tuples with the lower and upper bounds for
#'    each of P variables in the data frame, dimensions Px2.
#' @param intercept A logical vector of length one indicating whether an 
#'    intercept should be added prior to evaluating the inner product x'x.
#' @return The sensitivity of the data frame for which the covariance matrix
#'   is being calculated.
covarianceSensitivity <- function(n, rng, intercept) { # THIS IS INCORRECT
  diffs <- apply(rng, 1, diff)
  if (intercept) { diffs <- c(1, diffs) }
  sensitivity <- c()
  for (i in 1:length(diffs)) {
    for (j in i:length(diffs)) {
      s <- ((n - 1) / n) * diffs[i] * diffs[j]
      sensitivity <- c(sensitivity, s)
    }
  }
  return(sensitivity)
}


#' Function to convert unique private covariances into symmetric matrix
#'
#' @param release Differentially private release of elements in lower triangle 
#'    of covariance matrix.
#' @param columns A character vector indicating columns in the private 
#'    covariance to be included in the output. Length should be equal to the 
#'    number of columns the user wants to include.
#' @return A symmetric differentially private covariance matrix.
#'
#' This function is used by the mechanism in post-processing and not intended for interactive post-processing
covarianceFormatRelease <- function(release, columns) {
  outDim <- length(columns)
  outMatrix <- matrix(0, nrow=outDim, ncol=outDim)
  outMatrix[lower.tri(outMatrix, diag=TRUE)] <- release
  outMatrix[upper.tri(outMatrix, diag=FALSE)] <- t(outMatrix)[upper.tri(outMatrix, diag=FALSE)]
  release <- data.frame(outMatrix)
  rownames(release) <- names(release) <- columns
  return(release)
}


#' Function to perform linear regression using private release of covariance matrix
#'
#' @param release Differentially private release of elements in lower triangle 
#'    of covariance matrix.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data frame.
#' @param intercept Logical indicating whether the intercept is included in 
#'    \code{release}.
#' @param formula A list of the regression equations to be performed on the 
#'    covariance matrix.
#' @return Linear regression coefficients and standard errors for all specified
#'    \code{formula}.
covariancePostLinearRegression <- function(release, n, intercept, formula) {
  outSummaries <- vector('list', length(formula))
  for (f in 1:length(formula)) {
    outSummaries[[f]] <- linearReg(formula[[f]], release, n, intercept)
  }
  return(outSummaries)
}


#' Lower triangle of covariance matrix
#'
#' This function is called by \code{dpCovariance} and
#' produces the true value to be perturbed.
#'
#' @param x Data frame
#' @param columns Columns to include in output
#' @param intercept Logical, should the intercept be included?

funCovar <- function(x, intercept) {
  if (intercept) { x <- cbind(1, x) }
  covariance <- t(as.matrix(x)) %*% as.matrix(x)
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

dpCovariance <- setRefClass(
  Class = 'dpCovariance',
  contains = 'mechanismLaplace'
)

dpCovariance$methods(
  initialize = function(mechanism, varType, n, epsilon, columns, rng=NULL, imputeRng=NULL, 
                        intercept=FALSE, formula=NULL) {
    # NOTE THIS IS CURRENTLY NOTE SECURE DUE TO SENS AND EPSILON NOT DISTRIBUTED PROPERLY
    .self$name <- 'Differentially private covariance matrix'
    .self$mechanism <- checkMechanism(mechanism, "mechanismLaplace")
    .self$varType <- checkVariableType(varType, c('integer', 'double', 'numeric', 'logical'))  #NEED TO CHANGE TO ALLOW MULTIPLE VARTYPES  
    .self$n <- checkN(n)  #NEED TO ADD THING SAYING HOW THIS IS N PER COLUMN
    .self$epsilon <- checkEpsilon(epsilon) #NEED TO CHANGE
    .self$rng <- checkRange(rng, expectedLength=length(columns)) 
    .self$sens <- covarianceSensitivity(n, rng, intercept)
    
    .self$formula <- formula
    .self$intercept <- intercept
    
    checkVariableType(class(formula), "formula")
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
  })

dpCovariance$methods(
  release = function(data) {
    x <- data[columns];
    .self$result <- export(mechanism)$evaluate(fun=funCovar, x=x, sens=sens, postFun=.self$postProcess,
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
