checkEmpty <- function(n, emptyOkay=FALSE){
  isEmpty <- is.null(n) || is.na(n)
  if (isEmpty && !emptyOkay){
    stop("Input may not be NA or NULL.")
  }
  else{
    return(n)
  }
}

checkNumeric1D <- function(n, emptyOkay=FALSE){
  checkEmpty(n, emptyOkay)
  if (is.numeric(n) || is.null(n) || is.na(n)){
    return(n)
  } else{
    errorStr <- paste("Input value of ", toString(n), "is not of type numeric.")
    stop(errorStr)
  }
}

#' Helper function that generates error message if non-numeric typed value \code{n}
#' is passed as input, and otherwise returns \code{n}.
#' 
#' If \code{emptyOkay=TRUE}, NA or NULL values are allowed in n. Otherwise, they will
#' raise an error. 
#'
#' @param n Input value of arbitrary type.
#' @param A boolean. True if NA or NULL values are allowed, FALSE otherwise.
#'
#' @return n or errors.
checkNumeric <- function(n, emptyOkay=FALSE){
  for(i in 1:length(n)){
    checkNumeric1D(n[i], emptyOkay)
  }
  return(n)
}

#' Check if input xs has length n
#' 
#' @param xs A vector, factor, or object for which length is defined.
#' @param n A numeric value of length 1. The expected length of xs.
#'
#' @return xs or error.
checkLength <- function(xs, n){
  len <- length(xs)
  if (len == n){
    return(xs)
  }
  else {
    errorStr <- paste("Input was expected to be of length ", toString(n), " but is instead of length ", toString(len))
    stop(errorStr)
  }
}

#' Helper function for helpN
#' 
#' Checks that n is always a positive whole number; n may be NA or NULL only if emptyOkay=TRUE.
#' 
#' @param n the input n (often number of data points in data set) from the user
#' @param emptyOkay A boolean. True if NA or NULL values are allowed, FALSE otherwise.
#' 
#' @return n, if n is a series of positive integers with expected length and only containing NAs or NULL values if allowed.
checkN1D <- function(n, emptyOkay=FALSE){
  isEmpty <- is.null(n) || is.na(n)
  if (emptyOkay && isEmpty){
    return(n)
  } else if (isEmpty){
    stop("Input n may not be NA or NULL.")
  }
  
  checkNumeric(n)
  if (!all(n > 0) || (n%%1 != 0)) {
    stop("n must be a positive whole number.")
  } else {
    return(n)
  }
}

#' Check validity of n
#' 
#' n(s) should always be a positive whole number, check the user's input
#' 
#' @param n the input n (often number of data points in data set) from the user, or an array of n's from the user
#' @param expectedLength Positive integer. The expected length of n.
#' @param emptyOkay A boolean. True if NA or NULL values are allowed, FALSE otherwise.
#' 
#' @return n, if n is a series of positive integers with expected length and only containing NAs or NULL values if allowed.

checkN <- function(n, expectedLength=1, emptyOkay=FALSE) {
  checkLength(n, expectedLength)
  for (i in 1:length(n)) {
    checkN1D(n[i], emptyOkay) 
  }
  return(n)
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

#' One-Dimensional Range Parameter Check
#' 
#' Helper function for checkRange. Checks if a supplied range is an ordered pair.
#' Coerces any vector of length two or greater into an ordered pair, and issues 
#' an error for shorter vectors.
#' 
#' If varType is 'logical', range is coerced to c(0,1).
#' 
#' If emptyOkay, no error will be raised if the input range is null or NA, and an NA value will be returned.
#' 
#' @param rng Range. A numeric vector that ought to be an ordered pair.
#' @param varType String describing the variable type. One of 'logical' or 'numeric'.
#' @param emptyOkay Boolean. TRUE if a null or NA range is allowed, default to FALSE.
#' @return Vector of length 2 or NULL
checkRange1D <- function(rng, varType, emptyOkay=FALSE) {
  
  #Special case for logical variables.
  if (varType == 'logical') {
    rng <- c(0,1)
    return(rng)
  } 
  
  # Check for null or NA range
  lengthRng <- length(rng)
  emptyFlag <- is.null(rng) || (is.na(rng) && lengthRng <=1)
  if (emptyFlag && !emptyOkay){
    stop("Input range may not be empty.")
  } else if (emptyFlag){
    return(NULL)
  }
  
  rngStr <- paste('c(',toString(rng),')')
  
  # Check for NA values within range
  naFlag <- NA %in% rng
  if (naFlag && emptyOkay){
    warningStr <- paste('Range argument provided', rngStr,'has NA value. Setting range to NULL.')
    warning(warningStr)
    return(NULL)
  }else if (naFlag){
    stop("Input range may not contain NA values.")
  }
  
  # Range validation
  if (lengthRng < 2) {
    errorStr <- paste('Error in range argument provided,', rngStr, ': requires upper and lower values as vector of length 2.')
    stop(errorStr)
  } else if (lengthRng > 2) {
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
#' Checks if a supplied range(s) is(are) an ordered pair(s).
#' 
#' If emptyOkay, no error will be raised if the input range(s) is(are) null or NA. Ranges that
#' were input as NULL will be output as NULL.
#' 
#' In order to handle potential of some rows with different lengths when some rows contain NULL
#' or NA values, a rng input of 
#' 
#' @param rng Range. A numeric vector of one of two formats:
#'   1. A vector of length two, that ought to be an atomic ordered pair, representing 
#'   the maximum and minimum bounds on the data.
#'   2. A sequence of ordered pairs as a matrix or as a list, where each row represents
#'    the maximum and minimum bounds on some subsets of the data (e.g. of different data columns)
#'    Matrix and list types are supported. Internally, matrices are coerced to lists to allow 
#'    varying dimensions across rows.
#' @param varType The variable type; e.g. 'logical', 'integer', 'numeric', 'character'.
#' @param formatType One of 'vector' or 'list', describing which of the two the input range should be, where 'vector' returns to an atomic vector.
#'    Since matrices are coerced to lists within the function, when using a matrix as range input, `format type = 'list'` should be specified. 
#' @param expectedLength Integer value. Specifies how long the output ought to be. Defaults to NULL and only used on list or matrix inputs.
#' @param emptyOkay Boolean. TRUE if a null or NA range is allowed. Defaults to FALSE.
#' @return An ordered pair, a list of ordered pairs, or NULL.
#' 
#' Note that you can input a single ordered pair as a first element of a list, e.g. \code{rng = list(c(1,2))},
#' but performance will be slightly worse.
#' 
#' @examples
#'
#' checkRange(c(1,3))
#' checkRange(1:3)s
#' \dontrun{checkRange(1)}
#' @rdname checkRange
checkRange <- function(rng, varType, formatType, expectedLength=NULL, emptyOkay=FALSE) {
  # If input is matrix, coerce to a list. 
  if (is.matrix(rng)) {
    if (is.matrix(rng)){
      rng <- matrixToList(rng)
      formatType <- "list"
    }
  }
  
  # If input is just an atomic vector, just run validation checks on that.
  if ((isVector(rng) && formatType == "vector") || is.null(rng)) {
    rng <- checkRange1D(rng, varType, emptyOkay)
    return(rng)
  }
  
  else if ((is.list(rng)) && formatType == "list"){
    
    checkLength(rng, expectedLength)
    
    # Run validation on each tuple input in list
    newRngs <- rep(list(NA), length(rng)) #initialization
    
    for (i in 1:length(rng)) {
      newRngVal <- checkRange1D(rng[[i]], varType, emptyOkay)
      # If newRngVal is NULL, have to do some weird massaging to add it to the list.
      if (is.null(newRngVal)){
        newRngs[i] <- list(NULL)
      } else{
        newRngs[[i]] <- newRngVal
      }
    }
    return(newRngs)
  }
  else{
    stop("Input rng is not of expected format.")
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
checkEpsilon <- function(epsilon, expectedLength=1) {
  checkNumeric(epsilon)
  
  if (length(epsilon) > 1 && expectedLength<=1) {
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
checkAccuracy <- function(accuracy, expectedLength=1){
  checkNumeric(accuracy)
  checkLength(accuracy, expectedLength)
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
#' @export
checkVariableType <- function(type, inTypes) { 
  type <- tolower(type)
  inTypes <- tolower(inTypes)
  if (!(type %in% inTypes)) {
    stop(paste('Variable type', type, 'should be one of', paste(inTypes, collapse = ', ')))
  } 
  return(type)
}

checkMechanism <- function(inMech, allowedMechs) {
  checkLength(inMech, 1)
  checkEmpty(inMech)
  checkVariableType(typeof('inMech'), 'character')
  
  mech <- tolower(inMech)
  indexMech <- match(mech, tolower(allowedMechs)) #get index of mech in allowedMechs.
  # match returns NA if mech not in tolower(allowedMechs).
  if (is.na(indexMech)){
    stop(paste('Input mechanism', mech, 'should be one of', paste(allowedMechs, collapse = ', ')))
  }
  return(allowedMechs[indexMech])
}

#' Helper function that converts a matrix to a list, s.t. each element of the list is a row of the matrix.
#'
#' @param m An arbitrary matrix
#'
#' @return List form of m, where ith element of list is the ith row of the matrix m.
matrixToList <- function(m){
  getRow <- function(i){
    return(m[i,])
  }
  outList <- lapply(1:nrow(m), getRow)
  return(outList)
}

#' Checks if input value is an atomic vector (rather than a list).
#'
#' @param x arbitrary input
#'
#' @return TRUE if x is an atomic vector, FALSE otherwise.
isVector <- function(x){
  if (is.vector(x) && !is.list(x)){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}
