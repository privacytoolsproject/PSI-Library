#' Censoring data
#' 
#' For numeric types, checks if x is in rng = (min, max) and censors values to 
#'    either min or max if it is out of the range. For categorical types, 
#'    values not in `levels` are coded NA.
#'
#' @param x A vector of numeric or categorial values to censor.
#' @param varType Character indicating the variable type of \code{x}.
#'    Possible values include: numeric, logical, ...
#' @param rng For numeric vectors, a vector (min, max) of the bounds of the 
#'    range. For numeric matrices with nrow N and ncol P, a Px2 matrix of 
#'    (min, max) bounds.
#' @param levels For categorical types, a vector containing the levels to 
#'    be returned.
#' 
#' @return Original vector with values outside the bounds censored to the bounds.
#' @examples
#' 
#' censorData(x=1:10, varType='integer', rng=c(2.5, 7))
#' censorData(x=c('a', 'b', 'c', 'd'), varType='character', levels=c('a', 'b', 'c'))
#' @rdname censorData
#' @export
censorData <- function(x, varType, rng=NULL, levels=NULL, rngFormat=NULL) {
    if (varType %in% c('character', 'factor')) {
        if (is.null(levels)) {
            x <- factor(x, exclude=NULL)
        } else {
            x <- factor(x, levels=levels, exclude=NULL)
        }
    } else {
        if (NCOL(x) > 1) {
            checkRange(rng, varType, rngFormat, expectedLength=ncol(x)) 
            for(i in 1:NCOL(x)){
                x[,i] <- censorData1D(x[,i], rng[[i]])
            }
        } else if (rngFormat=="vector"){
            checkRange(rng, varType, rngFormat)
            x <- censorData1D(x,rng)
        } else {
          stop("range Format (rngFormat) must be either 'list' or 'vector'.")
        }
    }
    return(x)
}


censorData1D <- function(x, rng){
  x[x < rng[1]] <- rng[1]
  x[x > rng[2]] <- rng[2]
  return(x)
}


#' Logical variable check
#' 
#' Utility function to verify that a variable is dichotomous.
#'
#' This function effectively allows the user to ask for any variable containing
#' at most two unique values to treat the variable as logical. If the variable
#' contains numeric values, the highest value is recoded 1 and the the lower
#' value is recoded 0. If the variable is categorical and contains only two unique
#' values, the least frequently observed is recoded 1.
#'
#' @param x Vector containing two unique types of values.
#' 
#' @return Logical form of \code{x} coded 0-1.
#' @examples
#' 
#' makeLogical(sample(c('cat', 'dog'), size=8, replace=TRUE))
#' makeLogical(sample(c(0, 1), size=8, replace=TRUE))
#' makeLogical(sample(c(-6.87, 3.23), size=8, replace=TRUE))
#' @rdname makeLogical
#' @export
makeLogical <- function(x) {
    if (!length(unique(x)) <= 2) { # how to handle if contains 1 value only?
        stop('Variable has more than two values')
    }
    if (class(x) %in% c('character', 'factor')) {
        tab <- data.frame(table(x))
        min_level <- tab[which.min(tab[, 2]), 1]
        x <- ifelse(x == min_level, 1, 0)
    } else {
        x <- ifelse(x == max(x), 1, 0)
    }
    return(x)
}

#' Scaling helper function for fillMissing
#' 
#' Takes input array and scales to upper and lower bounds, which are either defined by lower and upper or calculated depending on
#' the type of variable. (Note that input array will always be numeric; varType refers to the variable type of the input array in
#' the fillMissing function.)
#' 
#' @param vals Input array of values to scale Type: numeric array
#' @param varType Variable type of input array to fillMissing function; affects how scaling occurs. 
#'   Type: one of following strings: 'character', 'factor', 'logical', 'integer', 'numeric'.
#' @param lower Lower bound for scaling. Type: numeric
#' @param upper Upper bound for scaling. Type: numeric
#' @param categories List of categories. Type: factor
#'
#' @return Array of values, either characters, integers, logicals, numerics depending on varType, scaled according to either the 
#' number of categories if varType='factor' or 'character', or based on lower and upper when varType='logical','numeric', or 'integer'.
#'
scaleValues = function(vals, varType, lower=NULL, upper=NULL, categories=NULL) {
    if (varType %in% c('character', 'factor')) { 
        lower <- 1
        upper <- length(categories)
    }
    
    if (varType == 'logical') {       # logical values can only be 0 or 1 so set bounds accordingly
        lower <- 0
        upper <- 2                       # upper bound of 2 not 1 because as.integer always rounds down.
    }
    
    out <- vals * (upper - lower) + lower  # scale uniform random numbers based on upper and lower bounds
    
    if (varType %in% c('logical', 'integer')) { # if logical or integer, trim output to integer values
        out <- as.integer(out)
    } else if(varType == 'logical') {
        
    } else if (varType %in% c('character', 'factor')) { # if character or factor, assign output to categories.
        out <- categories[as.integer(out)]
    }
    
    return(out)
}


#' Helper function for fillMissing; fills missing values in one-dimensional array
#'
#' Imputes uniformly in the range of the provided variable.
#'
#' @param x Input array of missing values.
#' @param varType Variable type of input array to fillMissing function; affects how scaling occurs. 
#'   Type: one of following strings: 'character', 'factor', 'logical', 'integer', 'numeric'.
#' @param lower Lower bound for scaling. Type: numeric. Default NULL.
#' @param upper Upper bound for scaling. Type: numeric. Default NULL.
#' @param categories List of categories. Type: factor. Default NULL.
#' @return Vector \code{x} with missing values imputed
fillMissing1D <- function(x, varType, lower=NULL, upper=NULL, categories=NULL) {
    naIndices <- is.na(x)         # indices of NA values in x
    nMissing <- sum(naIndices)    # number of missing values
    
    if (nMissing == 0) { 
        return(x) 
    }
    
    u <- dpUnif(nMissing) # array of uniform random numbers of length nMissing
    scaledVals <- scaleValues(u, varType, lower, upper, categories) # scale uniform vals
    x[naIndices] <- scaledVals #assign to NAs in input array
    return(x)
}

#' Helper function for fillMissing. Fills missing values column-wise for matrix.
#'
#' Impute uniformly in the range of the provided variable
#'
#' @param x Numeric matrix with missing values
#' @param varType Variable type of input array.
#'   Type: one of following strings: 'character', 'factor', 'logical', 'integer', 'numeric'.
#' @param imputeRng Px2 matrix where the pth row contains the range
#'      within which the pth variable in x is imputed.
#' @return Matrix \code{x} with missing values imputed
#'
#' @seealso \code{\link{fillMissing}}
fillMissing2D <- function(x, varType, imputeRng=NULL) {
    for (j in 1:ncol(x)) {
        x[, j] <- fillMissing1D(x[, j], varType, lower=imputeRng[j, 1], upper=imputeRng[j, 2])
    }
    return(x)
}


#' Fill missing values
#'
#' Impute uniformly in the range of the provided variable
#'
#' @param x Vector with missing values to impute
#' @param varType Character specifying the variable type
#' @param lower Numeric lower bound, default NULL
#' @param upper Numeric upper bound, default NULL
#' @param categories Set of possible categories from which to choose,
#'      default NULL
#' @return Vector \code{x} with missing values imputed
#' @examples
#'
#' # numeric example
#' y <- rnorm(100)
#' y[sample(1:100, size=10)] <- NA
#' y_imputed <- fillMissing(x=y, varType='numeric', lower=-1, upper=1)
#'
#' # categorical example
#' cats <- as.factor(c('dog', 'zebra', 'bird', 'hippo'))
#' s <- sample(cats, size=100, replace=TRUE)
#' s[sample(1:100, size=10)] <- NA
#' s_imputed <- fillMissing(x=s, varType='factor', categories=cats)
#'
#' @seealso \code{\link{dpUnif}}
#' @rdname fillMissing
#' @export
fillMissing = function(x, varType, imputeRng=NULL, categories=NULL) {
    if (varType %in% c('numeric', 'integer', 'logical')) {
        if (NCOL(x) > 1) {
            x <- fillMissing2D(x, varType, imputeRng)
        } else {
            x <- fillMissing1D(x, varType, imputeRng[1], imputeRng[2])
        }
    } else {
        x <- fillMissing1D(x, varType, categories=categories)
    }
    return(x)
}


#' Extract function arguments
#' 
#' Utility function to match arguments of a function with input list and/or attributes of input object
#'
#' @param targetFunc A character vector containing the name of the function to find the matching inputs of
#'    with arguments that need to be filled by output.
#' @param inputList A named list of values
#' @param inputObject An object with a predefined getFields() function that takes all attributes of object and returns them as a list
#' @return List of arguments and values needed for specification of \code{target.func}
#' @rdname getFuncArgs
getFuncArgs <- function(targetFunc, inputList=NULL, inputObject=NULL){
    funcArgs <- list()                      # Initialize list of function arguments
    argNames <- names(formals(targetFunc))  # Names of all arguments of targetFunc
    inputListNames <- names(inputList)      # Names of all items in inputList
    
    if (!is.null(inputObject)){             
        inputObjectNames <- names(inputObject$getFields()) # Names of all fields of inputObject
    }
    else{
        inputObjectNames <- NULL
    }
    
    for (argument in argNames) {
        # If argument of targetFunc has an associated named attribute in the input list, save that attribute to funcArgs
        if (argument %in% inputListNames) {
            funcArgs[[argument]] <- inputList[[argument]] 
        }
        # If argument of targetFunc has an associated field of inputObject, save that field value to funcArgs 
        if (argument %in% inputObjectNames) {          
            funcArgs[[argument]] <- inputObject[[argument]]
        }
    }
    return(funcArgs)
}
