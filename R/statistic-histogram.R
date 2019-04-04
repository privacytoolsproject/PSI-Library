#' Histogram accuracy
#' 
#' Determine accuracy of histogram release, given epsilon and delta, for the differentially 
#' private histogram release.
#'
#' @param n.bins A numeric vector of length one specifying the number of cells 
#'    in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data.
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param stability A logical vector indicating whether the stability 
#'    mechanism is used.
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

histogram.getAccuracy <- function(n.bins, n, epsilon, stability, delta=10^-6, alpha=0.05, error=1e-10) {
 	acc <- NULL
	if(stability){
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
#' @param n.bins A numeric vector of length one specifying the number of cells 
#'    in which to tabulate values.
#' @param n A numeric vector of length one specifying the number of
#'    observations in the data.
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param stability A logical vector indicating whether the stability 
#'    mechanism is used.
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

histogram.getParameters <- function(n.bins, n, accuracy, stability, delta=10^-6, alpha=0.05, error=1e-10) {
	eps <- NULL
	if(stability){
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
        release <- apply(release, 2, histogram.compose, n=n)
        release <- data.frame(t(release))
    } else {
        bin.names <- names(release)
        if (anyNA(bin.names)) { bin.names[is.na(bin.names)] <- 'NA' }
        release <- histogram.compose(release, n)
        release <- data.frame(matrix(release, ncol=length(release)))
    }
    names(release) <- bin.names
    return(release)
}


#' Constrain the sum of histogram bins to sample size
#'
#' @param h Histogram
#' @param n Sample size

histogram.compose <- function(h, n) {
    h <- as.numeric(h)
    j <- length(h)
    suppress.index <- sample(1:j, size=1)
    h[suppress.index] <- 0
    h[h < 0] <- 0
    h <- round(h)
    h[suppress.index] <- n - sum(h)
    return(h)
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


#' Differentially private histogram
#'
#' @param mechanism Character, the mechanism used to perturb histogram bins.
#' @param var.type Character, the variable type.
#' @param variable Character, the variable name in the data frame.
#' @param n Integer, the number of observations.
#' @param epsilon Numeric, the privacy loss parameter.
#' @param accuracy Numeric, the desired accuracy of the query.
#' @param rng Numeric, a priori estimate of the lower and upper bounds of a
#'    variable taking numeric values. Ignored for categorical types.
#' @param bins Character, the available bins or levels of a categorical variable.
#'    Ignored for numeric types.
#' @param n.bins Integer, the number of bins to release.
#' @param alpha Numeric, level of statistical significance, default 0.05.
#' @param delta Numeric, probability of privacy loss beyond \code{epsilon}.
#' @param error Numeric, error.
#'
#' @import methods
#' @export dpHistogram
#' @exportClass dpHistogram
#'
#' @include mechanism.R
#' @include mechanism-laplace.R
#' @include mechanism-bootstrap.R

dpHistogram <- setRefClass(
    Class = 'dpHistogram',
    contains = c('mechanismLaplace', 'mechanismBootstrap')
)

dpHistogram$methods(
    initialize = function(mechanism, var.type, variable, n, epsilon=NULL, accuracy=NULL, rng=NULL, 
                          bins=NULL, n.bins=NULL, alpha=0.05, delta=2^-30, error=1e-9,
                          impute.rng=NULL, impute=FALSE, n.boot=NULL, ...) {
        .self$name <- 'Differentially private histogram'
        .self$mechanism <- mechanism
        .self$variable <- variable
        .self$var.type.orig <- .self$var.type <- var.type
        .self$n <- n
        .self$rng <- rng
        .self$alpha <- alpha
        .self$delta <- delta
        .self$error <- error
        .self$impute <- impute
        if (var.type %in% c('numeric', 'integer')) {
            if (is.null(n.bins)) {
                stop('number of bins must be specified')
            }
            .self$n.bins <- n.bins
            .self$bins <- seq(rng[1], rng[2], length.out=(n.bins + 1))
            .self$stability <- FALSE
        } else if (var.type == 'logical') {
            .self$n.bins = ifelse(impute, 2, 3)
            if (!impute) {
                .self$bins = c(0, 1, NA)
                .self$var.type = 'factor'
            } else {
                .self$bins = c(0, 1)
            }
            .self$stability = FALSE
        } else {
            .self$bins <- bins
            .self$n.bins <- length(bins)
            .self$stability <- ifelse(is.null(bins), TRUE, FALSE)
        }
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- histogram.getParameters(n.bins, n, accuracy, stability, delta, alpha, error)
        } else {
            .self$epsilon <- epsilon
            .self$accuracy <- histogram.getAccuracy(.self$n.bins, n, epsilon, stability, delta, alpha, error)
        }
        if (is.null(impute.rng)) {
            .self$impute.rng <- rng
        } else {
            .self$impute.rng <- impute.rng
        }
        .self$boot.fun <- boot.hist
        .self$n.boot <- n.boot
})

dpHistogram$methods(
    release = function(data) {
        x <- data[, variable]
        noisy <- export(mechanism)$evaluate(fun.hist, x, 2, .self$postProcess)
        ### NEED TO ADD PAPER CITATION FOR STABILITY MECHANISM!
        if (stability) {
            if (check_histogram_n(noisy$accuracy, n, n.bins, epsilon, delta, alpha)) {
               # a <- accuracy * n / 2  
               # noisy$release <- noisy$release[noisy$release >= a]
               #JM changed to below after conversation with Victor
               a <- 1+2*log(2/delta)/epsilon 
               noisy$release <- noisy$release[noisy$release > a]
            }
        } #else {
          #  noisy$release <- round(noisy$release)
          #  noisy$release[noisy$release < 0] <- 0
        #}
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
