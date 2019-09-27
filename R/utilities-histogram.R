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
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#'    Default to 10^-6.
#' @param alpha A numeric vector of length one specifying the numeric 
#'    statistical significance level. Default to 0.05.
#' 
#' @export histogramGetAccuracy
#' @return Accuracy guarantee for histogram release, given epsilon.
#' @rdname histogramGetAccuracy

histogramGetAccuracy <- function(mechanism, epsilon, delta=2^-30, alpha=0.05, sensitivity) {
    acc <- NULL
    if(mechanism == 'mechanismStability'){
        acc <- stabilityGetAccuracy(sensitivity, epsilon, delta, alpha)
    }
    else{
        acc <- laplaceGetAccuracy(sensitivity, epsilon, alpha)
    }
    return(acc)
}


#' Histogram epsilon
#' 
#' Function to find the epsilon value necessary to meet a desired level of accuracy for the
#' differentially private histogram release.
#' 
#' This calculation is the inverse of the calculation for `histogramGetAccuracy`.
#' 
#' @seealso \code{\link{histogramGetAccuracy}} for accuracy derivation
#' 
#' @param mechanism A string indicating the mechanism that will be used to construct the histogram
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#'    Default to 10^-6.
#' @param alpha A numeric vector of length one specifying the numeric 
#'    statistical significance level. Default to 0.05.
#' 
#' @export histogramGetEpsilon
#' @return Differential privacy parameter epsilon
#' @rdname histogramGetEpsilon

histogramGetEpsilon <- function(mechanism, accuracy, delta=10^-6, alpha=0.05, sensitivity) {
    eps <- NULL
    if(mechanism == 'mechanismStability'){
        eps <- stabilityGetEpsilon(sensitivity, accuracy, delta, alpha)
    }
    else{
        eps <- laplaceGetEpsilon(sensitivity, accuracy, alpha)
    }
    return(eps)
}

#' Get accuracy for a stable statistic (only used for histogram statistic)
#' 
#' Function to find the accuracy guarantee of a stable statistic release at a given epsilon 
#' value.
#' 
#' @param sensitivity the sensitivity of the statistic
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#'    Default to 2^-30.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @return Accuracy guarantee for statistic release given epsilon.

stabilityGetAccuracy <- function(sensitivity, epsilon, delta = 2^-30, alpha=0.05) {
    accuracy <- (sensitivity/epsilon) * log(2/(alpha*delta)) + 1
    return(accuracy)
}


#' Get epsilon for a stable statistic (only used for histogram statistic)
#' 
#' Function to find the epsilon value necessary to meet a desired level of 
#' accuracy for a stable statistic release.
#' 
#' @param sensitivity the sensitivity of the statistic
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param delta The probability of an arbitrary leakage of information from 
#'    the data. Should be of length one and should be a very small value. 
#'    Default to 2^-30.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @return The scalar epsilon necessary to guarantee the needed accuracy.

stabilityGetEpsilon <- function(sensitivity, accuracy, delta = 2^-30, alpha=0.05) {
    epsilon <- sensitivity * log(2/(alpha*delta)) / (accuracy - 1)
    return(epsilon)
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
#' @param nBins A numeric vector of length one specifying the number of cells 
#'    in which to tabulate values.
#' @param accuracy A numeric vector representing the accuracy guarantee.
#'    
#' @return Confidence interval for the noisy counts in each bin.
#' @rdname histogramGetCI

histogramGetCI <- function(release, nBins, accuracy) {
    release_asNumeric <- as.numeric(release)
    out <- list()
    for (k in 1:nBins) {
        binCount <- release_asNumeric[k]
        if (binCount == 0) {
            out[[k]] <- c(0, accuracy)
        } else {
            out[[k]] <- c(max(0, binCount - accuracy), accuracy + binCount)
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
#' @param release Table, the result of \code{funHist}
#' @param n Sample size
#' @return Data frame

normalizeReleaseAndConvertToDataFrame <- function(release, n) {
    if (is(release, 'matrix')) {
        binNames <- rownames(release)
        if (anyNA(binNames)) { 
            binNames[is.na(binNames)] <- 'NA' 
        }
        release <- apply(release, 2, normalizeHistogram, n=n)
        release <- data.frame(t(release))
    } else {
        binNames <- names(release)
        if (anyNA(binNames)) { 
            binNames[is.na(binNames)] <- 'NA'
        }
        release <- normalizeHistogram(release, n)
        release <- data.frame(matrix(release, ncol=length(release)))
    }
    names(release) <- binNames
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
#' @rdname histogramPostHerfindahl
histogramPostHerfindahl <- function(release, n) {
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
#' @rdname histogramGetJSON
histogramGetJSON <- function(output.json=TRUE) {
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
#' @param varType Character indicating the variable type.
#' @param bins Vector indicating the bins into which \code{x} is partitioned.
#' @return Histogram with counts for each level of \code{x}.

funHist <- function(x, varType, bins=NULL) {
    if (varType %in% c('numeric', 'integer')) {
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
#' @param varType Character indicating the variable type.
#' @param bins Vector indicating the bins into which \code{xi} is partitioned.
#' @return Histogram with counts for each level of \code{xi}.

bootHist <- function(xi, varType, bins=NULL) {
    histogram <- funHist(xi, varType, bins)
    return(histogram)
}

#' Determine Mechanism
#' 
#' This is a set of helper functions to determine which mechanism to use when
#' calculating the histogram (either the stability mechanism or the Laplace
#' mechanism).
#' 
#' @param varType Character, the variable type.
#' @param rng Numeric, a priori estimate of the lower and upper bounds of a
#'    variable taking numeric values. Ignored for categorical types. Maybe be 
#'    null for numeric or integer types, in which case the stability mechanism is used.
#' @param bins Character or numeric, the available bins or levels of a variable. 
#'    Character for categorical variables, a vector of numbers for numeric variables.
#' @param nBins Integer, the number of bins to release.
#' @param granularity Numeric, the width of each histogram bin.
#' 
#' @return a string indicating the mechanism to use when constructing the differentially private histogram

determineMechanism <- function(varType, rng, bins, nBins, granularity) {
    if (!is.null(bins)) {
        # if the bins are specified, then the user has previous knowledge
        # of the data, so the stability mechanism is not necessary
        return('mechanismLaplace')
    } else {
        # if the bins are not specified, then we need to look at the
        # variable type to determine which mechanism to use
        return(determineMechanismByVariableType(varType, rng, bins, nBins, granularity))
    }
}

# only called by determineMechanism()
determineMechanismByVariableType <- function(varType, rng, bins, nBins, granularity) {
    if (varType == 'logical') {
        # if the variable type if logical, we will never use the
        # stability mechanism because the user already knows the
        # possible values of the data.
        return('mechanismLaplace')
    } else if (varType == 'categorical') {
        # if we have reached this conditional statement, then we
        # already know that the bins have not been specified. If
        # the data is categorical and the bins are not specified,
        # then the mechanism must be the stability mechanism.
        return('mechanismStability')
    } else if (varType %in% c('numeric', 'integer')) {
        # if the variable type is numeric or integer, then the user
        # must enter either the number of bins or the granularity
        # of the histogram bins. If neither is entered, then an
        # error must be thrown. Otherwise, we determine the mechanism
        # based on whether the user has prior knowledge of the range
        # of the data.
        if (is.null(nBins) & is.null(granularity)) {
            stop('number of bins or granularity must be specified')
        } else {
            return(determineMechanismByRange(varType, rng, bins, nBins, granularity))
        }
    } else {
        # if the variable type is none of the above or it is null,
        # then use the stability mechanism by default
        return('mechanismStability')
    }
}

# only called by determineMechansim()
determineMechanismByRange <- function(varType, rng, bins, nBins, granularity) {
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
#' If the user inputs a list of bins, the input bins will override the data and will
#' be released as the histogram bins. If a given bin does not exist in the data, it
#' will still be released in the result. It is possible that this non-existent bin
#' will still have a count, because it will be an option during data imputation in 
#' the call to `fillmissing()`. If the input list of bins does not include a value 
#' that exists in the data, the existing value will be changed to `NA` in the call 
#' to `censorData()` and will then be imputed as one of the input bins in `fillMissing()`.
#' 
#' @param varType Character, the variable type.
#' @param rng Numeric, a priori estimate of the lower and upper bounds of a
#'    variable taking numeric values. Ignored for categorical types. Maybe be null for 
#'    numeric or integer types, in which case the stability mechanism is used.
#' @param bins Character or numeric, the available bins or levels of a variable. Character 
#'    for categorical variables, a vector of numbers for numeric variables.
#' @param nBins Integer, the number of bins to release.
#' @param impute Boolean, if true then the mechanism should replace missing values with known 
#'    values from the data.If false, the mechanism should leave missing values as `NA`
#' @param granularity Numeric, the width of each histogram bin, or the number of observations in each bin
#' @param object Object, the dpHistogram object for the given variable (used it access and assign variable type)
#' 
#' @return a vector of histogram bins. Character vector for categorical variables. Numeric 
#'    vector for logical, numeric, and integer variables.

determineBins <- function(varType, rng, bins, n, nBins, impute, granularity, object) {
    if (!is.null(bins)) {
        # if the user passed in bins, then the passed bins are the histogram bins
        # check entered bins for errors. If there are not errors, entered bins will be assigned as histogram bins.
        # if there are errors, an error message will be returned to the user.
        errorCheckBins(varType, rng, bins)
        return(bins)
    } else {
        # if we have reached this condition, then the data is not categorical,
        # because we only call this function if the mechanism is not the stability
        # mechanism and the bins are not passed in by the user. If the variable is
        # categorical and the bins are not passed in, then the mechanism is the
        # stability mechanism. So we only need to check logical, numeric, and integer.
        if (varType == 'logical') {
            return(determineLogicalBins(impute, object))
        } else {
            # if we have reached this conditional statement, then the variable type can
            # only be numeric or integer, and the user must have entered the range and
            # either the number of bins or the granularity.
            return(determineNumericIntegerBins(rng, n, nBins, granularity))
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
determineNumericIntegerBins <- function(rng, n, nBins, granularity) {
    # first check if nBins is NULL, nBins is considered the truth for the number
    # of bins if the user has entered both nBins and granularity.
    if (is.null(nBins)) {
        # if the nBins is null, then the the granularity
        # must be entered
        nBinsFromGranularity <- n / granularity # get the number of bins
        return(seq(rng[1], rng[2], length.out=(nBinsFromGranularity + 1)))
    } else {
        # if nBins is not null, then we can use it to calculate
        # the histogram bins
        return(seq(rng[1], rng[2], length.out=(nBins + 1)))
    }
}

# only called by determineBins()
# checks given bins (only if bins are not null) to confirm they are within range of given data
errorCheckBins <- function(varType, rng, bins) {
    errorCheckBinVariableType(varType, bins)
    errorCheckBinRange(varType, rng, bins)
}

# only called by erroCheckBins()
errorCheckBinVariableType <- function(varType, bins) {
    # if variable type if character (categorical), confirm that given bins are character
    if (varType == 'character') {
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
    if (varType %in% c('numeric', 'integer')) {
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
    if (varType == 'logical') {
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
errorCheckBinRange <- function(varType, rng, bins) {
    # if the datatype is numeric or integer AND a range was entered,
    # check that each bin is within the range
    if ((varType %in% c('numeric', 'integer')) & (!is.null(rng))) {
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
#'   has a 'varType' member that will be changed from 'logical' to 'factor'
#' 
#' @return No return value

setVariableTypeAsFactor <- function(object) {
    object$varType <- 'factor'
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
#' @param varType The variable type for the histogram data
#' 
#' @return a list of the imputation bins that will be used for `fillMissing()`.

checkImputationBins <- function(imputationBins, bins, varType) {
    # if no imputation bins were entered, return the bins.
    # (Note:  bins may be NULL, in which case stability mechanism will be used.)
    if (is.null(imputationBins)) {
        return(bins)
    }
    
    # if imputation bins were entered, check that they are
    # within the list of bins. If not, do not use them in imputation.
    clippedBins <- c()
    if (varType %in% c('logical','factor','character')) {
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

setHistogramRange <- function(rng, varType, bins) {
    if (varType == 'logical') {
        return(c(0,1))
    } else if (varType %in% c('numeric','integer') & !(is.null(bins))) {
        binRange = c(bins[1],bins[length(bins)])
        return(binRange)
    } else {
        return(rng)
    }
}

#' Check the histogram variable type entered by the user 
#' 
#' The variable type for a histogram must be 'numeric', 'integer', 'logical', or 'character'.
#' If it is not one of these, send an error message to the user.
#' 
#' @param varType The variable type of the data that was entered by the user
#' 
#' @return No return value, will only send an error message if variable type is invalid.
checkHistogramVariableType <- function(varType) {
    if (!(varType %in% c("numeric", "integer", "logical", "character"))) {
        stop("Please enter a data type of 'numeric', 'integer', 'logical', or 'character'")
    }
}

#' Set the number of histogram bins 
#' 
#' Set the number of bins for the histogram, given either the entered number of bins,
#' the entered granularity, or the vector of bins.
#' 
#' @param nBins the number of bins entered by the user, may be null
#' @param granularity the number of items to be in each bin (i.e. the height of each bin), may be null
#' @param varType The variable type of the data that was entered by the user
#' @param bins the bin vector either entered by the user or set by determineBins()
#' 
#' @return No return value, will only send an error message if variable tyep is invalid.
setNumHistogramBins <- function(nBins, granularity, varType, bins) {
    if (varType %in% c('numeric', 'integer')) {
        # if the variable type is numeric or integer, then the length of the bins vector will
        # be 1 larger than the number of bins. So if the user did not enter a number of bins,
        # calculate the number of bins by subtracting one from the length of the bins vector.
        if (is.null(nBins)) {
            return(length(bins) - 1)
        } else {
            return(nBins)
        }
    } else {
        # if the variable type is not numeric or integer, then the nuber of bins can either be entered
        # by the user, calculated from the granularity, or calculate the number of bins from the length
        # of the bin vector.
        if (is.null(nBins)) {
            # if there is no input for number of bins, get it from from granularity or the list of bins
            if (!is.null(granularity)) {
                return(n / granularity)
            } else {
                return(length(bins))
            }
        } else {
            return(nBins)
        }
    }
}


#' Utility function to include NA level for categorical types when vector of bins
#' does not include all observed levels in the data vector.
#'
#' @param x Vector, categorical type
#' @param bins Vector, depositor-provided list of levels for which to count values

histogramCategoricalBins <- function(x, bins) {
    x <- factor(x, levels=bins, exclude=NULL)
    return(x)
}


#' Check histogram bins argument
#' 
#' Utility function to check bins argument to histogram. If number of bins 
#'    is not provided, the Sturges method is used.
#' 
#' @param nBins The number of cells in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#'
#' @return Number of bins
#' @rdname checkHistogramNBins
checkHistogramNBins <- function(nBins, n) {
    if (!is.null(nBins)) { # nBins may be null, in which case we do not want to change it
        if (nBins < 2) {
            stop('number of bins must be at least 2')
        } else if (as.logical(nBins %% 1)) {
            warning('non-integer value for number of bins converted to next highest integer value')
            return(ceiling(nBins))
        } else {
            return(nBins) # if no issues with input value, return input value
        }
    } else {
        return(nBins) # return the input value if the input nBins is null
    }
}
