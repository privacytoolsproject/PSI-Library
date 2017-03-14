#' Function to evaluate the mean
#' 
#' @param x Numeric vector
#' @return List with fields `name` specifying the statistic and `stat` with the value of the statistic

dp.mean <- function(x) { 
    value <- sum(x) / length(x)
    out <- list('name' = 'mean', 'stat' = value)
    return(out)
} 


#' Function to evaluate a histogram for a numeric variable
#' 
#' @param x Vector of numeric values
#' @param var.levels Vector specifying the bins
#' @return List with fields `name` specifying the statistic and `stat` with the value of the statistic

dp.histogram.numeric <- function(x, hist.type, var.levels) { 
    values <- table(cut(x, breaks=var_levels, include.lowest=TRUE, right=TRUE))
    out <- list('name' = 'histogram', 'stat' = values)
    return(out)
}


#' Function to evaluate a histogram for a categorical variable
#' 
#' @param x Vector of categorical values
#' @param var.levels Vector specifying the bins
#' @return List with fields `name` specifying the statistic and `stat` with the value of the statistic

dp.histogram.categorical <- function(x, var.levels) { 
    n <- length(var.levels)  # avoid unused argument error in mechanism
    values <- table(x, useNA='ifany')
    out <- list('name' = 'histogram', 'stat' = values)
    return(out)
}
