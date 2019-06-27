#' Histogram accuracy
#' 
#' Determine accuracy of histogram release, given epsilon and delta, for the differentially 
#' private histogram release.
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
#' @param error The error term of the statistical significance level. Default
#'    to 1e-9. 
#' 
#' @export histogram.getAccuracy
#' @return Accuracy guarantee for histogram release, given epsilon.
#' @rdname histogram.getAccuracy

#JM replaced below with getaccuracy function from dpmodules/Jack/Histogramnew.R

histogram.getAccuracy <- function(mechanism, n.bins, n, epsilon, delta=10^-6, alpha=0.05, error=1e-10) {
 	acc <- NULL
	if(mechanism == 'mechanismStability'){
		acc <- 2*log(2/(alpha*delta)) /epsilon
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
#' @param error The error term of the statistical significance level. Default
#'    to 1e-9. 
#' 
#' @export histogram.getParameters
#' @return Differential privacy parameter epsilon
#' @rdname histogram.getParameters

histogram.getParameters <- function(mechanism, n.bins, n, accuracy, delta=10^-6, alpha=0.05, error=1e-10) {
	eps <- NULL
	if(mechanism == 'mechanismStability'){
		eps <- 2*log(2/(alpha*delta)) /accuracy
	}
	else{
		eps <- 2*log(1/alpha) /accuracy
	}
	return(eps)
}

#' Histogram confidence interval
#' 
#' Return the confidence interval for the differentially private histogram release given the
#' accuracy.
#'
#' @param release A numeric vector with a noisy estimate of bin counts.
#' @param n.bins A numeric vector of length one specifying the number of cells 
#'    in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#'    
#' @return Confidence interval for the noisy counts in each bin.
#' @rdname histogram.getCI

histogram.getCI <- function(release, n.bins, n, accuracy) {
    release <- as.numeric(release)
    accxn <- accuracy * n
    out <- list()
    for (k in 1:n.bins) {
        bin.count <- release[k]
        if (bin.count == 0) {
            out[[k]] <- c(0, accxn)
        } else {
            out[[k]] <- c(max(0, bin.count - accxn), accxn + bin.count)
        }
    }
    out <- data.frame(do.call(rbind, out))
    #names(out) <- c('lower', 'upper')
    rownames(out) <- names(release)
    return(out)
}


#' Format the release of private histogram
#'
#' Convert the release from a table to a data frame
#'
#' @param release Table, the result of \code{fun.hist}
#' @param n Sample size
#' @return Data frame

histogram.formatRelease <- function(release, n) {
    if (is(release, 'matrix')) {
        bin.names <- rownames(release)
        if (anyNA(bin.names)) { bin.names[is.na(bin.names)] <- 'NA' }
        release <- apply(release, 2, normalizeHistogram, n=n)
        release <- data.frame(t(release))
    } else {
        bin.names <- names(release)
        if (anyNA(bin.names)) { bin.names[is.na(bin.names)] <- 'NA' }
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
#' @param granularity Numeric, the width of each histogram bin.
#' @param object Object, the dpHistogram object for the given variable (used it access and assign variable type)
#' 
#' @return a vector of histogram bins. Character vector for categorical variables. Numeric 
#'    vector for logical, numeric, and integer variables.

determineBins <- function(var.type, rng, bins, n.bins, impute, granularity, object) {
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
            return(determineNumericIntegerBins(rng, n.bins, granularity))
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
determineNumericIntegerBins <- function(rng, n.bins, granularity) {
	# first check if n.bins is NULL, n.bins is considered the truth for the number
	# of bins if the user has entered both n.bins and granularity.
    if (is.null(n.bins)) {
        # if the n.bins is null, then the the granularity
        # must be entered
    	r <- rng[2] - rng[1] # get the width of the range
    	nBinsFromGranularity <- r / granularity # get the number of bins
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
        warning('You have entered both bins and a data range, when you do not need both.
            Default is to use the bins that have been entered.
            If you would like to use the range, please enter the range and the desired number of bins and omit the bins.')
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
#' @param alpha Numeric, level of statistical significance, default 0.05.
#' @param delta Numeric, probability of privacy loss beyond \code{epsilon}.
#' @param error Numeric, error.
#' @param granularity Numeric, the width of each histogram bin (i.e. the inverse of `n.bins`). Used 
#'    to calculate histogram bins in comination with `rng`.
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
                          bins=NULL, n.bins=NULL, alpha=0.05, delta=2^-30, error=1e-9,
                          impute.rng=NULL, impute=FALSE, n.boot=NULL, granularity=NULL, ...) {
        .self$name <- 'Differentially private histogram'
        
        # determine  which mechanism to use based on inputs
        .self$mechanism <- determineMechanism(var.type, rng, bins, n.bins, granularity)
        
        # set parameters of the histogram
        .self$var.type.orig <- .self$var.type <- var.type
        .self$variable <- variable
        .self$n <- n
        .self$epsilon <- epsilon
        .self$accuracy <- accuracy
        .self$rng <- rng # may be null
        .self$bins <- bins
        .self$n.bins <- n.bins # may be null
        .self$alpha <- alpha
        .self$delta <- delta
        .self$error <- error
        .self$impute.rng <- impute.rng
        .self$impute <- impute
        .self$n.boot <- n.boot
        .self$granularity <- granularity # may be null
        .self$boot.fun <- boot.hist
        
        # if the mechanism used is NOT the stability mechanism, determine the
        # bins of the histogram. If the mechanism is the stability mechanism,
        # then the bins will be determined in the stability mechanism.
        # Once the bins are determined, get the number of bins.
        if (.self$mechanism != 'mechanismStability') {
            .self$bins <- determineBins(var.type, rng, bins, n.bins, impute, granularity, .self)
            .self$n.bins <- ifelse(is.null(.self$n.bins), length(.self$bins), .self$n.bins)
        }
        
        # get the epsilon and accuracy
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- histogram.getParameters(mechanism, n.bins, n, accuracy, delta, alpha, error)
        } else {
            .self$epsilon <- epsilon
            .self$accuracy <- histogram.getAccuracy(mechanism, .self$n.bins, n, epsilon, delta, alpha, error)
        }
        
        # set the range for data imputation
        if (is.null(impute.rng)) {
            .self$impute.rng <- rng
        } else {
            .self$impute.rng <- impute.rng
        }
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
        out$release <- histogram.formatRelease(out$release, n)
        out$accuracy <- accuracy
        out$epsilon <- epsilon
        if (length(out$release) > 0) {
            if (mechanism == 'mechanismLaplace') {
                out$interval <- histogram.getCI(out$release, n.bins, n, out$accuracy)
            }
        }
        if (var.type %in% c('factor', 'character')) {
            out$herfindahl <- sum((out$release / n)^2)
        }
        if (var.type.orig == 'logical') {
            temp.release <- out$release[na.omit(names(out$release))]
            out$mean <- as.numeric(temp.release[2] / sum(temp.release))
            out$median <- ifelse(out$mean < 0.5, 0, 1)
            out$variance <- out$mean * (1 - out$mean)
            out$std.dev <- sqrt(out$variance)
        }
        return(out)
})
