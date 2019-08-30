#' Check if input is numeric type
#' 
#' Helper function that generates error message if non-numeric typed value \code{n}
#' is passed as input, and otherwise returns \code{n}.
#'
#' @param n Input value of arbitrary type.
#'
#' @return n or errors.
checkNumeric <- function(n){
  if (!is.numeric(n)){
    errorStr <- paste("Input value of ", toString(n), "is not of type numeric.")
    stop(errorStr)
  }
  else{
    return(n)
  }
}

#' Check validity of n
#' 
#' n should always be a positive whole number, check the user's input
#' 
#' @param n the input n from the user
#' 
#' @return n, if n is a positive whole number

checkN <- function(n) {
  checkNumeric(n)
  if ((n <= 0) || (n%%1 != 0)) {
    stop("n must be a positive whole number.")
  } else if (length(n) > 1) {
    stop("n must be of length 1.")
  }
  else {
    return(n)
  }
}

#' One-Dimensional Range Parameter Check
#' 
#' Helper function for checkRange. Checks if a supplied range is an ordered pair.
#' Coerces any vector of length two or greater into an ordered pair, and issues 
#' an error for shorter vectors.
#' 
#' @param rng Range. A numeric vector that ought to be an ordered pair.
#' @param varType String describing the variable type. One of 'logical' or 'numeric'.
#' 
checkRange1D <- function(rng, varType) {
  rngStr <- paste('c(',toString(rng),')')
  if (varType == 'logical') {
    rng <- c(0,1)
  } else if (length(rng) < 2) {
    errorStr <- paste('Error in range argument provided,', rngStr, ': requires upper and lower values as vector of length 2.')
    stop(errorStr)
  } else if (length(rng) > 2) {
    warningStr <- paste('Range argument of', rngStr, 'has more than two values.  Will proceed using min and max values as range.')
    warning(warningStr)
    rng <- c(min(rng), max(rng))
  } else {
    rng <- sort(rng)
  }
  return(rng)
}

#' Range Parameter Check
#' 
#' Checks if a supplied range is an ordered pair.
#'    
#' @param rng Range. A numeric vector of one of two formats:
#'   1. Of length two, that ought to be an ordered pair, representing 
#'   the maximum and minimum bounds on the data
#'   2. A sequence of ordered pairs in columns, where each row represents
#'    the maximum and minimum bounds on some subsets of the data (e.g. of different data columns)
#' @param varType The variable type; e.g. 'logical', 'integer', 'numeric', 'character'.
#' @return An ordered pair.
#' @examples
#'
#' checkRange(c(1,3))
#' checkRange(1:3)s
#' \dontrun{checkRange(1)}
#' @rdname checkRange
checkRange <- function(rng, varType) {
  if (NCOL(rng) > 1){
    newRng <- matrix(NA, nrow(rng), 2)
    for (i in 1:nrow(rng)) {
      newRng[i,] <- checkRange1D(rng[i,], varType)
    }
    return(newRng)
  }
  else{
    rng <- checkRange1D(rng, varType)
    return(rng)
  }
}

#' Epsilon Parameter Check
#' 
#' Utility function for checking that epsilon is acceptably defined.
#'
#' @param epsilon A vector, that ought to be positive.
#' @param multipleEps A boolean value. If TRUE, multiple epsilon paramaters may be defined. Default to FALSE.
#' @param expectedLength Integer value. Specifies how long the output ought to be.
#' 
#' @return The supplied epsilon if acceptable, otherwise an error 
#'    message interupts.
#'    
checkEpsilon <- function(epsilon, multipleEps=FALSE, expectedLength=1) {
  checkNumeric(epsilon)
  if (length(epsilon) > 1 && !multipleEps) {
    stop(paste("Privacy parameter epsilon must be a single value, but is currently a vector of length", length(epsilon)))
  }
  for (eps in epsilon){
    if (eps <= 0) {
      stop("Privacy parameter epsilon must be a value greater than zero.")
    }
    if (eps >= 3){
      strEps <- toString(eps)
      warningStr <- paste('Epsilon value of ', strEps, ' is in use. This is a higher global value than recommended for most cases.')
      warning(warningStr)
    }
  }
  actualLength <- length(epsilon)
  if (actualLength != expectedLength) {
    errorStr = paste("Epsilon parameter has improper length. Actual length is ", toString(actualLength),
                     " while expected length is ", toString(expectedLength),".")
    stop(errorStr)
  }
  return(epsilon)
}

#' Utility function for checking that accuracy is acceptably defined.
#' 
#' Verifies accuracy is greater than 0 and is a single value.
#'
#' @param accuracy 
#'
#' @return accuracy or errors.
checkAccuracy <- function(accuracy){
  checkNumeric(accuracy)
  if (!all(accuracy > 0)){
    stop("Accuracy must be greater than 0.")
  }
  return(accuracy)
}

#' Checking variable types
#' 
#' Verifies that the variable is an element in the set of acceptable types.
#' 
#' @param type A character specifying the type of the variable.
#' @param inTypes A vector of acceptable types of variables.
#' 
#' @return The original character string indicating the variable type.
#' @examples 
#' 
#' checkVariableType(type='Numeric', inTypes=c('Numeric', 'Factor'))
#' @rdname checkVariableType
checkVariableType <- function(type, inTypes) {
  type <- tolower(type)
  inTypes <- tolower(inTypes)
  if (!(type %in% inTypes)) {
    stop(paste('Variable type', type, 'should be one of', paste(inTypes, collapse = ', ')))
  } 
  return(type)
}

#' Error check imputation range for numeric or integer variables
#' 
#' Check that the entered imputation range is within the entered data range. If yes, return
#' the entered imputation range, which will be used as the imputation range for the call
#' to the utility function `fillMissing()`. If not, return the data range. 
#' If the imputation range is NULL, default to the data range.
#' 
#' We check if the imputation range is within the data range because it is a privacy concern.
#' If the imputation range is outside of the data range, NA values will be replaced with values 
#' outside of the data range, which will show that there are NA values in the data or skew the 
#' result when the differentially private estimate is released.
#' 
#' @param imputationRange The imputation range entered by the user
#' @param rng The data range entered by the user
#' @param varType The variable type for the histogram data
#' 
#' @return the imputation range that will be used for `fillMissing()`.

checkImputationRange <- function(imputationRange, rng, varType) {
  # if no imputation range was entered, return the data range.
  if (is.null(imputationRange)) {
    return(rng)
  }
  
  # for numeric and integer variables, the imputation range should be a 2-tuple
  # with the minimum and maximum of the imputation range.
  # if an imputation range was entered, check that it is
  # within the data range. If it is not, clip it to be within the data range
  if (varType %in% c('numeric', 'integer')) {
    
    checkNumeric(imputationRange)
    
    if (length(imputationRange)!=2){
      stop("Imputation range must have length 2.")
    }
    
    lowerBound <- NULL
    upperBound <- NULL
    
    # if the imputation range lower bound is below the data range lower bound,
    # clip the lower bound to the data range
    if (imputationRange[1] < rng[1]) {
      warning('Lower bound of imputation range is outside of the data range. Setting lower bound of the imputation range to the lower bound of the data range.')
      lowerBound <- rng[1]
    } else {
      lowerBound <- imputationRange[1]
    }
    
    # if the imputation range upper bound is above the data range upper bound,
    # clip the upper bound to the data range
    if (imputationRange[2] > rng[2]) {
      warning('Upper bound of imputation range is outside of the data range. Setting upper bound of the imputation range to the upper bound of the data range.')
      upperBound <- rng[2]
    } else {
      upperBound <- imputationRange[2]
    }
    
    for (entry in imputationRange) {
      if (!is.numeric(entry)) {
        warning('Imputation range for a numeric variable must be numeric. Setting imputation range to data range.')
        lowerBound <- rng[1]
        upperBound <- rng[2]
      }
    }
    
    # return the (potentially clipped) imputation range
    return(c(lowerBound,upperBound))
    
  } else {
    # if the variable type is something other than numeric or integer,
    # default to the data range
    warning('Imputation range entered for variable that is not of numeric or integer type. Setting imputation range to data range.')
    return(rng)
  }
}