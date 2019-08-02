#' Linear regression on covariance matrix
#' 
#' Function to extract regression coefficients from the covariance matrix via 
#'    the sweep operator.
#'
#' @param formula An object of class 'formula' containing the desired 
#'    regression formula.
#' @param release A numeric private release of the covariance matrix.
#' @param n A numeric vector of length one specifying the number of
#'    observations in in the data.
#' @param intercept A logical vector indicating whether the intercept is 
#'    included in \code{formula}.
#'    
#' @return A numeric vector of regression coefficients corresponding 
#'    to \code{formula}.
#' @rdname linear.reg
linear.reg <- function(formula, release, n, intercept) {
    if (!all(eigen(release)$values > 0)) {  # could do is.positive.definite() but that requires matrixcalc package
        coefs <- "The input matrix is not invertible"
        return(coefs)
    } else {
        xy.locs <- extract.indices(as.formula(formula), release, intercept)
        x.loc <- xy.locs$x.loc
        y.loc <- xy.locs$y.loc
        loc.vec <- rep(TRUE, (length(x.loc) + 1))
        loc.vec[y.loc] <- FALSE
        sweep <- amsweep((as.matrix(release) / n), loc.vec)
        coefs <- sweep[y.loc, x.loc]
        se <- sqrt(sweep[y.loc, y.loc] * diag(solve(release[x.loc, x.loc])))
        coefs <- data.frame(cbind(coefs, se))
        coefs <- format(round(coefs, 5), nsmall=5)
        rownames(coefs) <- xy.locs$x.label
        names(coefs) <- c('Estimate', 'Std. Error')
        return(coefs)
    }
}


#' Moore Penrose Inverse Function
#' 
#' @references Gill, Jeff, and Gary King. "What to do when your Hessian is not invertible: Alternatives to model respecification in nonlinear estimation." Sociological methods & research 33, no. 1 (2004): 54-87.
#' 
#' Generate the Moore-Penrose pseudoinverse matrix of \code{X}.
#' 
#' @param X A numeric, symmetric covariance matrix.
#' @param tol Convergence requirement. 

mpinv <- function(X, tol = sqrt(.Machine$double.eps)) {
    ## Moore-Penrose Inverse function (aka Generalized Inverse)
    ##   X:    symmetric matrix
    ##   tol:  convergence requirement
    s <- svd(X)
    e <- s$d
    e[e > tol] <- 1/e[e > tol]
    s$v %*% diag(e,nrow=length(e)) %*% t(s$u)
}


#' Sweep operator
#' 
#' General sweep operator citation:
#' @references Goodnight, James H. "A tutorial on the SWEEP operator." The American Statistician 33, no. 3 (1979): 149-158.
#' This implementation is from pseudocode from:
#' @references Schafer, Joseph L. Analysis of incomplete multivariate data. Chapman and Hall/CRC, 1997.
#' Code ported from:
#' @references Honaker, James, Gary King, and Matthew Blackwell. "Amelia II: A program for missing data." Journal of statistical software 45, no. 7 (2011): 1-47.
#' 
#' Sweeps a covariance matrix to extract regression coefficients.
#' 
#' @param g Each unit of a numeric, symmetric covariance matrix divided by
#'    the number of observations in the data the covariance matrix was 
#'    derived from.
#' @param m A logical vector of length equal to the number of rows in \code{g}
#'    in which the \code{TRUE} values correspond to the x values in the matrix
#'    and the \code{FALSE} values correspond to the y value in the matrix.
#' @param reverse Logical vector specifying whether the sign of the matrix 
#'    should be flipped. Default to \code{FALSE}.
#' 
#' @return The coefficients from \code{g}.
#' @rdname amsweep
amsweep <- function(g, m, reverse=FALSE) {
    if (identical(m, vector(mode='logical', length=length(m)))) {
        return(g)
    } else {
        p <- nrow(g)
        rowsm <- sum(m)
        if (rowsm == p) {
            h <- solve(g)
            h <- -h
        } else {
            kseq <- 1:p
            k <- kseq[m]
            kcompl <- kseq[-k]
            g11 <- g[k, k, drop=FALSE]
            g12 <- g[k, kcompl, drop=FALSE]
            g21 <- t(g12)
            g22 <- g[kcompl, kcompl, drop=FALSE]
            h11a <- try(solve(g11), silent=TRUE)
            if (inherits(h11a, "try-error")) {
                h11a <- mpinv(g11)
            }
            h11 <- as.matrix((-h11a))
            if (reverse) {sgn2 <- -1} else {sgn2 <- 1}
            h12 <- as.matrix(sgn2 * (h11a %*% g12))
            h21 <- as.matrix(t(h12))
            h22 <- g22 - (g21 %*% h11a %*% g12)
            hwo <- rbind(cbind(h11, h12), cbind(h21, h22))
            xordering <- c(k, kcompl)
            h <- matrix(0, p, p)
            h[xordering, xordering] <- hwo
        }
        return(h)
    }
}
