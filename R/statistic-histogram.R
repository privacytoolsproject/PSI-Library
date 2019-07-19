#' Histogram accuracy
#' 
#' Determine accuracy of histogram release, given epsilon and delta, for the differentially 
#' private histogram release.
#' 
#' In differential privacy, "accuracy" is defined as the threshold value above which a given value 
#' is "significantly different" from the expected value. Mathematically, this is written as:
#' \deqn{\alpha = Pr[Y > a]}
#' Where \eqn{\alpha} is the statistical significance level, \eqn{a} is the accuracy,
#' and \eqn{Y} is a random  variable indicating the difference between the differentially 
#' private noisy output and the true value. This equation is saying that with probability 
#' \eqn{1-\alpha}, the count of a hisotgram bin will be within \eqn{a} of the true count. 
#' 
#' The equation for \eqn{Y} is:
#' \deqn{Y = |X - \mu|}
#' Where \eqn{\mu} is the true value of a bin and \eqn{X} is the noisy count. \eqn{X}
#' follows a Laplace distribution centered at \eqn{\mu}. Subtracting \eqn{mu} centers
#' \eqn{Y} at \eqn{0}, and taking the absolute value "folds" the Lapalce distribution.
#' The absolute value is taken because the difference between noisy and true outputs 
#' is measured in magnitude.
#' 
#' Deriving the accuracy formula:
#' 
#' \enumerate{
#'     \item The probability density function (PDF) \eqn{f(x)} of the Laplace distribution is:\cr
#'           \eqn{f(x) = {1 / 2\lambda} * e^{-|x-\mu| / \lambda}}
#'     \item Using the definition of \eqn{Y} above, we can consider the differentially private PDF \eqn{g(Y)} to be:\cr
#'     \eqn{g(y) = {1 / \lambda} * e^{-y / \lambda}}
#'     \item Using \eqn{\alpha = Pr[Y > a]} and the PDF, we can solve for \eqn{a} and plug in \eqn{\lambda = 2 / \epsilon}, and end up with the accuracy formula: 
#'          \deqn{a = {2 / \epsilon} * ln(1 / \alpha)}
#'     \item The accuracy formula for the stability mechanism is derived by adding the accuracy formula above to the accuracy threshold (which is the worst-case potentially added noise in the stability mechanism): \eqn{{2 / \epsilon} * ln(2 / \delta)+1}
#' }
#' 
#' @references S. Vadhan The Complexity of Differential Privacy, Section 3.3 Releasing Stable Values p.23-24. March 2017.
#'
#' @param mechanism A string indicating the mechanism that will be used to construct the histogram
#' @param n.bins A numeric vector of length one specifying the number of cells 
#'    in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#'    Default to 10^-6.
#' @param alpha A numeric vector of length one specifying the numeric 
#'    statistical significance level. Default to 0.05.
#' 
#' @export histogram.getAccuracy
#' @return Accuracy guarantee for histogram release, given epsilon.
#' @rdname histogram.getAccuracy

histogram.getAccuracy <- function(mechanism, n.bins, n, epsilon, delta=10^-6, alpha=0.05) {
 	acc <- NULL
	if(mechanism == 'mechanismStability'){
		acc <- (2*log(2/(alpha*delta)) /epsilon) + 1
	}
	else{
		acc <- 2*log(1/alpha) /epsilon
	}
	return(acc)
}


#' Histogram epsilon
#' 
#' Function to find the epsilon value necessary to meet a desired level of accuracy for the
#' differentially private histogram release.
#' 
#' This calculation is the inverse of the calculation for `histogram.getAccuracy`.
#' 
#' @seealso \code{\link{histogram.getAccuracy}} for accuracy derivation
#' 
#' @param mechanism A string indicating the mechanism that will be used to construct the histogram
#' @param n.bins A numeric vector of length one specifying the number of cells 
#'    in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data.
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#'    Default to 10^-6.
#' @param alpha A numeric vector of length one specifying the numeric 
#'    statistical significance level. Default to 0.05.
#' 
#' @export histogram.getEpsilon
#' @return Differential privacy parameter epsilon
#' @rdname histogram.getEpsilon

histogram.getEpsilon <- function(mechanism, n.bins, n, accuracy, delta=10^-6, alpha=0.05) {
	eps <- NULL
	if(mechanism == 'mechanismStability'){
		eps <- 2*log(2/(alpha*delta)) /accuracy
	}
	else{
		eps <- 2*log(1/alpha) /accuracy
	}
	return(eps)
}

#' Histogram confidence intervals
#' 
#' Return the confidence interval for each bins of the differentially private histogram
#' release, given the accuracy.
#' 
#' A confidence interval indicates the range in which we estimate the true value of a
#' statistic to be. In this case, the confidence interval indicates the range in which
#' the true count for each histogram bin could be. To give an example: say a 
#' differentially private release of a histogram bucket has a count of 5. Say the 
#' confidence interval for that bucket is [3,7]. We can say "we are 95% confident that 
#' the true count of this histogram bin is between 3 and 7."
#'
#' @param release A numeric vector with a noisy estimate of bin counts.
#' @param n.bins A numeric vector of length one specifying the number of cells 
#'    in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#' @param accuracy A numeric vector representing the accuracy guarantee.
#'    
#' @return Confidence interval for the noisy counts in each bin.
#' @rdname histogram.getCI

histogram.getCI <- function(release, n.bins, n, accuracy) {
    release_asNumeric <- as.numeric(release)
    out <- list()
    for (k in 1:n.bins) {
        bin.count <- release_asNumeric[k]
        if (bin.count == 0) {
            out[[k]] <- c(0, accuracy)
        } else {
            out[[k]] <- c(max(0, bin.count - accuracy), accuracy + bin.count)
        }
    }
    out <- data.frame(do.call(rbind, out))
    names(out) <- c('lower bound', 'upper bound')
    rownames(out) <- names(release)
    return(out)
}


#' Convert normalized histogram release to data frame
#' 
#' Convert the release of the private histogram from a table to a data frame.
#' While converting, the normalizeHistogram method is called, which normalizes
#' the private histogram so that the sum of the bin counts is equal to n.
#'
#' @param release Table, the result of \code{fun.hist}
#' @param n Sample size
#' @return Data frame

normalizeReleaseAndConvertToDataFrame <- function(release, n) {
    if (is(release, 'matrix')) {
        bin.names <- rownames(release)
        if (anyNA(bin.names)) { 
            bin.names[is.na(bin.names)] <- 'NA' 
        }
        release <- apply(release, 2, normalizeHistogram, n=n)
        release <- data.frame(t(release))
    } else {
        bin.names <- names(release)
        if (anyNA(bin.names)) { 
            bin.names[is.na(bin.names)] <- 'NA'
        }
        release <- normalizeHistogram(release, n)
        release <- data.frame(matrix(release, ncol=length(release)))
    }
    names(release) <- bin.names
    return(release)
}


#' Constrain the sum of histogram bins to sample size
#'
#' @param h Histogram
#' @param n Sample size
#' 
#' @return The noisy histogram with the sum of the bin counts normalized to the input data size

normalizeHistogram <- function(h, n) {
    # get just the bins sizes for each bin in the histogram (without bin labels)
    binSizes <- as.numeric(h)
    
    # round all negative bins to 0 and all bins greater than n to n
    binSizes[binSizes < 0] <- 0
    binSizes[binSizes > n] <- n
    
    # get the total size of the histogram, it is not necessarily equal to n
    sumOfBins <- sum(binSizes)
    
    # get a vector of the bin weights
    binWeights <- binSizes / sumOfBins
    
    # get a vector of the normalized bin counts
    normalizedBinCounts <- binWeights * n
    
    return(normalizedBinCounts)
}

#' Histogram Herfindahl Index
#' 
#' Produce differentially private Herfindahl index for categorical types of data.
#' 
#' Let \eqn{s_i} be the percentage of data points in category \eqn{i} in the result histogram.
#' Then, the Herfindahl index is defined to be the sum of the squares of \eqn{s_i} over all
#' categories. Since the percentage of data points in category \eqn{i} may be computed
#' directly from the bin counts in a histogram, a differentially private Herfindahl index
#' may be calculated from a noisy histogram with no additional privacy loss incurred.
#'
#' @param release A numeric vector with a noisy estimate of bin counts.
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#'    
#' @return Herfindahl index.
#' @rdname histogram.postHerfindahl
histogram.postHerfindahl <- function(release, n) {
    share <- release / n
    herfindahl <- sum(share^2)
    return(herfindahl)
}


#' JSON doc for histogram
#' 
#' Produce a JSON doc for differentially private histograms.
#'
#' @param output.json Should the output be converted to JSON format. Default
#' to \code{TRUE}.
#'
#' @return JSON doc for histogram function.
#' @rdname histogram.getJSON
histogram.getJSON <- function(output.json=TRUE) {
    out <- list()
    out$statistic <- 'Histogram'
    out$description <- 'Differentially Private Histogram'
    out$mechanisms <- c('noisy', 'random', 'stability')
    out$variableTypes <- list('numeric' = list(), 'categorical' = list())
    out$variableTypes$numeric$rTypes <- c('numeric', 'integer')
    out$variableTypes$numeric$fields <- list(
            'n' = 'Number of observations',
            'range' = 'Ordered pair indicating effective lower and upper bounds',
            'n_bins' = 'Number of cells in output (optional, default Sturges method)'
        )
    out$variableTypes$categorical$rTypes <- c('character', 'factor')
    out$variableTypes$categorical$fields <- list(
            'n' = 'Number of observations',
            'bins' = 'Vector indicating levels for which to produce frequencies'
        )
    if (output.json) {
        out <- jsonlite::toJSON(out, pretty=TRUE)
    }
    return(out)
}


#' Histogram
#'
#' Function to evaluate a histogram for numeric and categorical types. This function
#' is used internally by \code{dpHistogram} to evaluate the true histogram prior to 
#' perturbation.
#'
#' @param x Vector of numeric or categorical type.
#' @param var.type Character indicating the variable type.
#' @param bins Vector indicating the bins into which \code{x} is partitioned.
#' @return Histogram with counts for each level of \code{x}.

fun.hist <- function(x, var.type, bins=NULL) {
    if (var.type %in% c('numeric', 'integer')) {
        histogram <- table(cut(x, breaks=bins, include.lowest=TRUE, right=TRUE))
    } else {
        histogram <- table(x, useNA='ifany')
    }
    return(histogram)
}

#' Bootstrap replication for histogram
#'
#' This is a wrapper for the histogram function used internally by the 
#' bootstrap mechanism.
#'
#' @param xi Bootstrapped vector of numeric or categorical type.
#' @param var.type Character indicating the variable type.
#' @param bins Vector indicating the bins into which \code{xi} is partitioned.
#' @return Histogram with counts for each level of \code{xi}.

boot.hist <- function(xi, var.type, bins=NULL) {
    histogram <- fun.hist(xi, var.type, bins)
    return(histogram)
}

#' Determine Mechanism
#' 
#' This is a set of helper functions to determine which mechanism to use when
#' calculating the histogram (either the stability mechanism or the Laplace
#' mechanism).
#' 
#' @param var.type Character, the variable type.
#' @param rng Numeric, a priori estimate of the lower and upper bounds of a
#'    variable taking numeric values. Ignored for categorical types. Maybe be 
#'    null for numeric or integer types, in which case the stability mechanism is used.
#' @param bins Character or numeric, the available bins or levels of a variable. 
#'    Character for categorical variables, a vector of numbers for numeric variables.
#' @param n.bins Integer, the number of bins to release.
#' @param granularity Numeric, the width of each histogram bin.
#' 
#' @return a string indicating the mechanism to use when constructing the differentially private histogram

determineMechanism <- function(var.type, rng, bins, n.bins, granularity) {
    if (!is.null(bins)) {
        # if the bins are specified, then the user has previous knowledge
        # of the data, so the stability mechanism is not necessary
        return('mechanismLaplace')
    } else {
        # if the bins are not specified, then we need to look at the
        # variable type to determine which mechanism to use
        return(determineMechanismByVariableType(var.type, rng, bins, n.bins, granularity))
    }
}

# only called by determineMechanism()
determineMechanismByVariableType <- function(var.type, rng, bins, n.bins, granularity) {
    if (var.type == 'logical') {
        # if the variable type if logical, we will never use the
        # stability mechanism because the user already knows the
        # possible values of the data.
        return('mechanismLaplace')
    } else if (var.type == 'categorical') {
        # if we have reached this conditional statement, then we
        # already know that the bins have not been specified. If
        # the data is categorical and the bins are not specified,
        # then the mechanism must be the stability mechanism.
        return('mechanismStability')
    } else if (var.type %in% c('numeric', 'integer')) {
        # if the variable type is numeric or integer, then the user
        # must enter either the number of bins or the granularity
        # of the histogram bins. If neither is entered, then an
        # error must be thrown. Otherwise, we determine the mechanism
        # based on whether the user has prior knowledge of the range
        # of the data.
        if (is.null(n.bins) & is.null(granularity)) {
            stop('number of bins or granularity must be specified')
        } else {
            return(determineMechanismByRange(var.type, rng, bins, n.bins, granularity))
        }
    } else {
        # if the variable type is none of the above or it is null,
        # then use the stability mechanism by default
        return('mechanismStability')
    }
}

# only called by determineMechansim()
determineMechanismByRange <- function(var.type, rng, bins, n.bins, granularity) {
    if (is.null(rng)) {
        # if the range is null, then the user has no prior knowledge of the
        # range of the data. Thus the stability mechanism must be used,
        # becuase it gives preserves privacy by giving plausible denialbility
        # about the range of the data (by removing buckets that have too low
        # of a count)
        return('mechanismStability')
    } else {
        # if we have gotten to this conditional statement, then we know the
        # number of bins or the granularity of the histogram bins has been
        # specified, and we also know that the range has been specified.
        # That means that the user has prior knowledge of the data, and we
        # do not need the stability mechanism to preserve privacy.
        return('mechanismLaplace')
    }
}

#' Determine Bins
#' 
#' Determine the bins of the histogram based on the inputs from the user
#' 
#' @param var.type Character, the variable type.
#' @param rng Numeric, a priori estimate of the lower and upper bounds of a
#'    variable taking numeric values. Ignored for categorical types. Maybe be null for 
#'    numeric or integer types, in which case the stability mechanism is used.
#' @param bins Character or numeric, the available bins or levels of a variable. Character 
#'    for categorical variables, a vector of numbers for numeric variables.
#' @param n.bins Integer, the number of bins to release.
#' @param impute Boolean, if true then the mechanism should replace missing values with known 
#'    values from the data.If false, the mechanism should leave missing values as `NA`
#' @param granularity Numeric, the width of each histogram bin, or the number of observations in each bin
#' @param object Object, the dpHistogram object for the given variable (used it access and assign variable type)
#' 
#' @return a vector of histogram bins. Character vector for categorical variables. Numeric 
#'    vector for logical, numeric, and integer variables.

determineBins <- function(var.type, rng, bins, n, n.bins, impute, granularity, object) {
    if (!is.null(bins)) {
        # if the user passed in bins, then the passed bins are the histogram bins
        # check entered bins for errors. If there are not errors, entered bins will be assigned as histogram bins.
        # if there are errors, an error message will be returned to the user.
        errorCheckBins(var.type, rng, bins)
        return(bins)
    } else {
        # if we have reached this condition, then the data is not categorical,
        # because we only call this function if the mechanism is not the stability
        # mechanism and the bins are not passed in by the user. If the variable is
        # categorical and the bins are not passed in, then the mechanism is the
        # stability mechanism. So we only need to check logical, numeric, and integer.
        if (var.type == 'logical') {
            return(determineLogicalBins(impute, object))
        } else {
            # if we have reached this conditional statement, then the variable type can
            # only be numeric or integer, and the user must have entered the range and
            # either the number of bins or the granularity.
            return(determineNumericIntegerBins(rng, n, n.bins, granularity))
        }
    }
}

# only called by determineBins()
determineLogicalBins <- function(impute, object) {
    if (impute) {
        # if the variable is logical and the user wants to impute empty data points,
        # then the histogram bins are only true and false
        return(c(0,1))
    } else {
        # if the variable is logical but the user does not want to impute empty
        # data points, then there needs to be a histogram bin for empty datapoints.
        # Because there is now a bin for NA values, the variable type should be 'factor',
        # not logical, because there is a non-logical bin.
        setVariableTypeAsFactor(object)
        return(c(0,1,NA))
    }
}

# only called by determineBins()
determineNumericIntegerBins <- function(rng, n, n.bins, granularity) {
	# first check if n.bins is NULL, n.bins is considered the truth for the number
	# of bins if the user has entered both n.bins and granularity.
    if (is.null(n.bins)) {
        # if the n.bins is null, then the the granularity
        # must be entered
    	nBinsFromGranularity <- n / granularity # get the number of bins
    	return(seq(rng[1], rng[2], length.out=(nBinsFromGranularity + 1)))
    } else {
        # if n.bins is not null, then we can use it to calculate
        # the histogram bins
    	return(seq(rng[1], rng[2], length.out=(n.bins + 1)))
    }
}

# only called by determineBins()
# checks given bins (only if bins are not null) to confirm they are within range of given data
errorCheckBins <- function(var.type, rng, bins) {
    errorCheckBinVariableType(var.type, bins)
    errorCheckBinRange(var.type, rng, bins)
}

# only called by erroCheckBins()
errorCheckBinVariableType <- function(var.type, bins) {
    # if variable type if character (categorical), confirm that given bins are character
    if (var.type == 'character') {
        # loop through all bins entered
        for (enteredBin in bins) {
            # check that each bin is of type `character`
            # if it is NOT, send error message to user
            if (!is.character(enteredBin)) {
                stop('Bins must be of type `character` for a variable of type `character`')
            }
        }
    }
    
    # if variable is numeric or integer, confirm that the given bins are numeric
    if (var.type %in% c('numeric', 'integer')) {
        # loop through all bins entered
        for (enteredBin in bins) {
            # check that each bin is of type `numeric` (integers are also numeric)
            # if it is NOT, send error message to user
            if (!is.numeric(enteredBin)) {
                stop('Bins must be numeric for a numeric variable')
            }
        }
    }
    
    # if variable is logical, conform that bins are only 0, 1, or NA
    if (var.type == 'logical') {
        # loop through all bins entered
        for (enteredBin in bins) {
            # check that each bins is only 0, 1, or NA
            # if it is NOT, send error message to user
            if (!(enteredBin %in% c(0,1,NA))) {
                stop('Histogram bins for a logical variable may only be 0, 1, or NA')
            }
        }
    }
}

# only called by errorCheckBins()
errorCheckBinRange <- function(var.type, rng, bins) {
    # if the datatype is numeric or integer AND a range was entered,
    # check that each bin is within the range
    if ((var.type %in% c('numeric', 'integer')) & (!is.null(rng))) {
        # if the user user has both specified bins and entered a range,
        # show an error message, because we do not need both. Default to
        # the bins entered.
        warning('You have entered both bins and a data range, when you do not need both. Default is to use the bins that have been entered. If you would like to use the range, please enter the range and the desired number of bins and omit the bins.')
    }
}

#' Set Variable Type as Factor
#' 
#' Set the variable type of the input variable to type 'factor'. This will only be used
#' when the input variable type of 'logical' and impute == FALSE.
#' 
#' @param object The dpHistogram object for the histogram on the specific variable. The object
#'   has a 'var.type' member that will be changed from 'logical' to 'factor'
#' 
#' @return No return value

setVariableTypeAsFactor <- function(object) {
    object$var.type <- 'factor'
}

#' Error check imputation bins for logical, factor, or character variables
#' 
#' Check that the entered imputation bins are a subset of the histogram bins. If yes, return
#' the entered imputation bins, which will be used as the imputation bins for the call
#' to the utility function `fillMissing()`. If not, return the histogram bins. 
#' If the imputation bins are NULL, default to the histogram bins. The histogram bins may be
#' NULL, in which case the stability mechanism will be used to determine the histogram bins and
#' imputation bins.
#' 
#' @param imputationBins The imputation bins entered by the user
#' @param rng The histogram bins entered by the user
#' @param var.type The variable type for the histogram data
#' 
#' @return a list of the imputation bins that will be used for `fillMissing()`.

checkImputationBins <- function(imputationBins, bins, var.type) {
    # if no imputation bins were entered, return the bins.
    # (Note:  bins may be NULL, in which case stability mechanism will be used.)
    if (is.null(imputationBins)) {
        return(bins)
    }
    
    # if imputation bins were entered, check that they are
    # within the list of bins. If not, do not use them in imputation.
    clippedBins <- c()
    if (var.type %in% c('logical','factor','character')) {
        # Loop through each bin in the imputation bins given by the user and check 
        # that they are in the list of histogram bins
        for (b in imputationBins) {
            if (b %in% bins) {
                # if the bin is in the list of histogram bins, add it to the list of imputation bins
                clippedBins <- c(clippedBins, b)
            } else {
                warning('Imputation bin entered is not in list of histogram bins, removing bins from imputation bins.')
            }
        }
        
        # check if the list of (potentially clipped) imputation bins is null
        # (this would be the case if all entered bins are outside of the histogram bins)
        # if yes, return the histogram bins, otherwise return the imputation bins
        ifelse(is.null(clippedBins), return(bins), return(clippedBins))
        
    } else {
        # if the variable type is something other than logical, factor, or character,
        # default to the histogram bins
        warning('Imputation bins entered for variable that is not of logical or categorical type.
                Setting imputation bins to histogram bins')
        return(bins)
    }
}

#' Check if delta value was entered
#' 
#' This method is to send the user a warning message if they entered a delta value
#' that will not be used.
#' 
#' If the mechanism is NOT the stabiltiy mechanism and the user entered a delta value, 
#' send a warning message saying the delta value with not be used and set the delta 
#' value to NULL. If the stability mechanism is being used and the user entered a 
#' delta value, set it as the delta value (the value of delta will be checked in the 
#' stability mechanism). If the stability mechanism is being used and the user did 
#' not enter a delta value, set the delta to the default value (2^-30).
#' 
#' @param mechanism The mechanism chosen by determineMechanism
#' @param delta The delta value entered by the user, may be NULL
#' 
#' @return The value that will be used as the delta value for the statistic, may be NULL.

checkDelta <- function(mechanism, delta) {
    # throw an error if the stability mechanism is NOT being used and the
    # user entered a delta value (because a delta value is only used in 
    # the stability mechanism) and return NULL so the delta value is ignored
    if (mechanism != 'mechanismStability') {
        if (!is.null(delta)) {
            warning('A delta parameter has been entered, but the stability mechanism is not being used. A delta value is only necessary for the stability mechanism. Entered delta value ignored.')
        }
        return(NULL)
    } else {
        # if the stability mechanism is being used, return the delta value
        if (is.null(delta)) {
            # default delta value
            return(2^-30)
        } else {
            return(delta)
        }
    }
}

#' Set histogram range if bins are entered 
#' 
#' If the user enters a logical variable, they should not need to enter a
#' histogram range, but calculating the sensitivity requires a range. This
#' method sets the range for a logical variable as c(0,1).
#' 
#' If the user enters numeric bins, they should not need to enter a histogram
#' range, but calculating the sensitivity requires a range. This method sets
#' the range for a numeric variable with bins entered to the bins lower bound
#' and the bin upper bound.
#' 
#' If neither of the cases are true for the histogram, the range the user entered
#' if returned.
#' 
#' @param rng The rng entered by the user, may be NULL
#' @param delta The variable type of the data
#' @param bins The histogram bins entered by the user, may be NULL
#' 
#' @return The 2-tuple that will be used as the range for the histogram when 
#' calculating sensitivity, censoring data, and imputing values.

setHistogramRange <- function(rng, var.type, bins) {
    if (var.type == 'logical') {
        return(c(0,1))
    } else if (var.type %in% c('numeric','integer') & !(is.null(bins))) {
        bin_range = c(bins[1],bins[length(bins)])
        return(bin_range)
    } else {
        return(rng)
    }
}

#' Check the histogram variable type entered by the user 
#' 
#' The variable type for a histogram must be 'numeric', 'integer', 'logical', or 'character'.
#' If it is not one of these, send an error message to the user.
#' 
#' @param var.type The variable type of the data that was entered by the user
#' 
#' @return No return value, will only send an error message if variable tyep is invalid.
checkHistogramVariableType <- function(var.type) {
    if (!(var.type %in% c("numeric", "integer", "logical", "character"))) {
        stop("Please enter a data type of 'numeric', 'integer', 'logical', or 'character'")
    }
}

#' Set the number of histogram bins 
#' 
#' Set the number of bins for the histogram, given either the entered number of bins,
#' the entered granularity, or the vector of bins.
#' 
#' @param n.bins the number of bins entered by the user, may be null
#' @param granularity the number of items to be in each bin (i.e. the height of each bin), may be null
#' @param var.type The variable type of the data that was entered by the user
#' @param bins the bin vector either entered by the user or set by determineBins()
#' 
#' @return No return value, will only send an error message if variable tyep is invalid.
setNumHistogramBins <- function(n.bins, granularity, var.type, bins) {
    if (var.type %in% c('numeric', 'integer')) {
        # if the variable type is numeric or integer, then the length of the bins vector will
        # be 1 larger than the number of bins. So if the user did not enter a number of bins,
        # calculate the number of bins by subtracting one from the length of the bins vector.
        if (is.null(n.bins)) {
            return(length(bins) - 1)
        } else {
            return(n.bins)
        }
    } else {
        # if the variable type is not numeric or integer, then the nuber of bins can either be entered
        # by the user, calculated from the granularity, or calculate the number of bins from the length
        # of the bin vector.
        if (is.null(n.bins)) {
            # if there is no input for number of bins, get it from from granularity or the list of bins
            if (!is.null(granularity)) {
                .return(n / granularity)
            } else {
                return(length(bins))
            }
        } else {
            return(n.bins)
        }
    }
}

#' Differentially private histogram
#'
#' @param var.type Character, the variable type.
#' @param variable Character, the variable name in the data frame.
#' @param n Integer, the number of observations.
#' @param epsilon Numeric, the privacy loss parameter.
#' @param accuracy Numeric, the desired accuracy of the query.
#' @param rng Numeric, a priori estimate of the lower and upper bounds of a
#'    variable taking numeric values. Ignored for categorical types.
#' @param bins Character, the available bins or levels of a categorical variable.
#' @param n.bins Integer, the number of bins to release.
#' @param granularity Numeric, the width of each histogram bin (i.e. the inverse of `n.bins`). Used 
#'    to calculate histogram bins in comination with `rng`.
#' @param alpha Numeric, level of statistical significance, default 0.05.
#' @param delta Numeric, probability of privacy loss beyond \code{epsilon}.
#' @param impute.rng Numeric, a 2-tuple indicating the lower and upper bounds of the range from which NA
#'    values in numeric or integer-type variables should be imputed 
#' @param impute.bins Character (or numeric for logical variables), a list of bins from which NA values
#'    values in character or logical-type variables should be imputed
#' @param impute Boolean, a boolean value indicating if logical-type variables should have NA values
#'    imputed or not. If true, a logical variable histogram will have 2 bins, 0 and 1. If false, the
#'    histogram will have 3 bins: 0, 1, and NA.
#' @param n.boot Numeric, the number of bootstrap iterations to do for bootstrapping (not used for version 1 release)
#' 
#'
#' @import methods
#' @export dpHistogram
#' @exportClass dpHistogram
#'
#' @include mechanism.R
#' @include mechanism-laplace.R
#' @include mechanism-stability.R

dpHistogram <- setRefClass(
    Class = 'dpHistogram',
    contains = c('mechanismLaplace', 'mechanismStability')
)

dpHistogram$methods(
    initialize = function(var.type, variable, n, epsilon=NULL, accuracy=NULL, rng=NULL, 
                          bins=NULL, n.bins=NULL, granularity=NULL, alpha=0.05, delta=NULL,
                          impute.rng=NULL, impute.bins=NULL, impute=FALSE, n.boot=NULL, ...) {
        .self$name <- 'Differentially private histogram'
        
        # check variable type, can only continue initialization for certain variable type: numeric, integer, logical, character
        checkHistogramVariableType(var.type)
        
        # determine  which mechanism to use based on inputs
        .self$mechanism <- determineMechanism(var.type, rng, bins, n.bins, granularity)
        
        # set parameters of the histogram
        .self$var.type <- var.type
        .self$variable <- variable
        .self$n <- n
        .self$epsilon <- epsilon
        .self$accuracy <- accuracy
        .self$bins <- bins # may be null
        .self$n.bins <- n.bins # may be null
        .self$alpha <- alpha
        .self$impute.rng <- impute.rng
        .self$impute <- impute
        .self$n.boot <- n.boot
        .self$granularity <- granularity # may be null
        .self$boot.fun <- boot.hist
        
        # if the mechanism used is NOT the stability mechanism:
        # 1) determine the bins of the histogram. (If the mechanism is 
        #    the stability mechanism, then the bins will be determined in 
        #    the stability mechanism.)
        # 2) determine the number of bins from the input number of bins, the granularity, or the list of bins.
        if (.self$mechanism != 'mechanismStability') {
            .self$bins <- determineBins(.self$var.type, rng, bins, .self$n, n.bins, impute, granularity, .self)
            .self$n.bins <- setNumHistogramBins(n.bins, granularity, .self$var.type, .self$bins)
        }
        
        # check the data range
        # if numeric bins have been entered, set the range to the range of the bins 
        # if logical variable is entered, set the range to c(0,1)
        # (may be NULL)
        .self$rng <- setHistogramRange(rng, .self$var.type, bins)
        
        # get the epsilon and accuracy
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- histogram.getEpsilon(mechanism, n.bins, n, accuracy, delta, alpha)
        } else {
            .self$epsilon <- epsilon
            .self$accuracy <- histogram.getAccuracy(mechanism, .self$n.bins, n, epsilon, delta, alpha)
        }
        
        # get the delta parameter (will be NULL unless stability mechanism is being used)
        .self$delta <- checkDelta(.self$mechanism, delta)
        
        # set the range for data imputation (will be null if no range entered)
        .self$impute.rng <- checkImputationRange(impute.rng, rng, var.type)
        
        # set the bins for data imputation (will be null if no bins entered)
        .self$impute.bins <- checkImputationBins(impute.bins, bins, var.type)
})

dpHistogram$methods(
    release = function(data) {
        x <- data[, variable]
        noisy <- export(mechanism)$evaluate(fun.hist, x, 2, .self$postProcess)
        .self$result <- noisy
})

dpHistogram$methods(
    postProcess = function(out) {
        out$variable <- variable
        out$release <- normalizeReleaseAndConvertToDataFrame(out$release, n)
        out$accuracy <- accuracy
        out$epsilon <- epsilon
        out$mechanism <- mechanism
        if (mechanism == 'mechanismStability') out$delta <- delta
        if (length(out$release) > 0) {
            if (mechanism == 'mechanismLaplace') {
                out$intervals <- histogram.getCI(out$release, n.bins, n, out$accuracy)
            }
        }
        if (var.type %in% c('factor', 'character')) {
            out$herfindahl <- sum((out$release / n)^2)
        }
        if (var.type %in% c('logical', 'factor')) {
            temp.release <- out$release[na.omit(names(out$release))]
            out$mean <- as.numeric(temp.release[2] / sum(temp.release))
            out$median <- ifelse(out$mean < 0.5, 0, 1)
            out$variance <- out$mean * (1 - out$mean)
            out$std.dev <- sqrt(out$variance)
        }
        return(out)
})
