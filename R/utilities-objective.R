#' Function to evaluate the p-norm of vectors in a matrix
#'
#' @param X matrix of numeric values
#' @param p The p in p-norm
#' @param margin The subscripts over which the norm is applied, where 1 gives row norms
#'      and 2 gives the column norms
#' @return A vector of norms

vectorNorm <- function(X, p=1, margin=1) {
    fun <- function(vec, p) {
        vec.norm <- sum(abs(vec)^p)^(1 / p)
        return(vec.norm)
    }
    if (!is.null(dim(X))) {
        p.norm <- apply(X, margin, fun, p)
    } else {
        p.norm <- fun(X, p)
    }
    return(p.norm)
}


#' Function to map rows of a numeric matrix to the unit ball
#'
#' This function will ensure that the row having the largest p-norm will be located on the unit ball
#'
#' @param X Numeric matrix
#' @param p The p in p-norm
#' @return List that includes the transformed matrix and the value corresponding to the maximum 
#'      observed p-norm

mapMatrixUnit <- function(X, p=1) {
    maxNorm <- max(vectorNorm(X, p=p))
    normed.matrix <- X / maxNorm
    return(list('matrix' = normed.matrix, 'maxNorm' = maxNorm))
}


#' Scale coefficient estimates
#'
#' This function puts coefficient estimates on the scale of the original features
#'
#' @param estimates Numeric, coefficient estimates vector.
#' @param xScaler Numeric, maximum norm from \code{mapMatrixUnit} fit on features.
#' @param yScaler Numeric, maximum norm from \code{mapMatrixUnit} fit on response,
#'    default NULL.
#' @return Transformed coefficients vector

scaleRelease <- function(estimates, xScaler, yScaler=NULL) {
    if (!is.null(yScaler)) {
        estimates <- estimates * yScaler
    }
    p <- length(estimates)
    estimates[2:p] <- estimates[2:p] / xScaler
    return(estimates)
}
