#' LaPlace mechanism for releasing differentially private queries
#'
#' @param fun A user supplied function, or string naming a function.  Must accept only one argument, named \code{x}.
#' @param x A vector of data to run the function on.
#' @param var_type Data type of vector x
#' @param n The number of samples
#' @param range An a priori estimate of the range
#' @param sensitivity numeric
#' @param epsilon numeric
#' @return Differentially private release of function \code{fun} on data \code{x}
#' @examples
#' n <- 1000
#' range <- c(0,1)
#' x <- runif(n, min=min(range), max=max(range))
#' sensitivity <- diff(range)/n
#' mechanism.laplace(mean,x=x,sensitivity=sensitivity, epsilon=.1)

#mechanism.laplace = function(fun, x, var_type, range, sensitivity, epsilon, levels=NULL, ...){
mechanism.laplace = function(fun, x, var_type, range, sensitivity, epsilon, ...) {

    epsilon <- checkepsilon(epsilon)

    if (hasArg(var_levels)) { 
        var_levels <- list(...)$var_levels
    } 

    if (var_type == 'logical') { # allow any dichotomous variable to be treated as logical type
        x <- make_logical(x)
    }

    if (var_type %in% c('numeric', 'integer', 'logical')) {
        range <- checkrange(range)
        x <- censordata(x, var_type, range=range)
    } else {
        x <- censordata(x, var_type, levels=var_levels)
    }

    truevalue <- fun(x, ...)
    #noise <- rlaplace(n=length(truevalue), sensitivity=sensitivity, epsilon=epsilon)
    noise <- rlap(mu=0, b=(sensitivity / epsilon), size=length(truevalue))
    release <- truevalue + noise
    return(release)
}


mechanism.histogram.random <- function(x, var_type, epsilon, levels, bins, n_bins) {

    if (var_type %in% c('factor', 'character')) {

        x <- censordata(x, var_type=var_type, levels=bins)
        bins <- unique(x)
        n_bins <- length(bins)

    } else {

        n_bins <- check_histogram_bins(n_bins=n_bins, n=n)
        x <- censordata(x, var_type=var_type, range=range)
        levels <- seq(range[1], range[2], length.out=(n_bins + 1))
    }

    if (n_bins < 2) {
        stop('mechanism `random` requires at least 2 bins for a numeric type')
    }

    if (epsilon > log((n_bins + 1) / (n_bins - 1))) {
        stop('`epsilon` is too large to guarantee differential privacy')
    }

    levels_mapper <- function(data, map, pr_same) {
        idx <- map[which(map[, 'level'] == data), 'index']
        if (runif(1) >= pr_same) {
            sampled_idx <- sample(1:(dim(map)[1] - 1), size=1)
            if (sampled_idx >= idx) {
                idx <- idx + 1
            }
        }
        return(map[idx, 'level'])
    }

    pr_diff <- 1 / n_bins * exp(-epsilon)
    pr_same <- 1 - (n_bins - 1) * pr_diff

    if (var_type %in% c('numeric', 'integer')) {

        levels_bins <- cut(x, breaks=levels, include.lowest=TRUE, right=TRUE)
        levels_uniq <- data.frame('level' = sort(unique(levels_bins)), 'index' = 1:n_bins)
        levels_out <- unlist(lapply(levels_bins, levels_mapper, levels_uniq, pr_same))

    } else {

        levels_uniq <- data.frame('level' = sort(unique(x), na.last=TRUE), 'index' = 1:n_bins)
        levels_out <- unlist(lapply(x, levels_mapper, levels_uniq, pr_same))
    }

    m <- 1 / (pr_same - pr_diff)
    b <- pr_diff * n * m
    release <- m * table(levels_out) - b
    return(release)
}


## NOTES ##
# [1] The user sometimes needs public information, like n, to determine the sensitivity.  We might check n is correct?  For example:
#	if(length(x) != n)  stop("Supplied value of n is not the length of the dataset.")  
# Not important at present, but there might become cases where argument x is not the data, but a reference to find the data.
#
# [2] Could possibly include an args argument with a NULL default as a way to override the args construction and supply a more general function.
#	if(is.null(args)){   
#		args = list(x=x)
#	}
# 
# [3] On line 26, should `n` be the length of `truevalue`? Number of draws equal to the number 
#   of values in the statistic. Cell counts for histogram, for example.
# 
# [4] Using ellipsis to pass flexible arguments to function in mechanism. Unsure how to do that 
#   in context of do.call(). Motivating case is to allow users to specify number of bins in histogram, 
#   also unsure if that is needed or appropriate functionality, or if we may want to allow this 
#   flexibility for other functions. 
# 
# [5] Including additional argument for data type. I still want to use this function for categorical types, but assume we do not 
#   need a numeric range or censor for these types. I am doing a censor of categorical types in the histogram.release function based
#   on the bins supplied to the function. Not sure if that belongs here.
