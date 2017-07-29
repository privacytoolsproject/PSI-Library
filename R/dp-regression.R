#' Function to perform regression using the differentially private covariance matrix via the sweep operator
#' 
#' @param formula Formula, regression formula used on data
#' @param release Numeric, private release of covariance matrix
#' @param n Integer, indicating number of observations
#' @param method Character vector, the regression technique used
#' @export

regression.release <- function(formula, release, n, method='linear') {
  intercept <- ifelse('intercept' %in% names(release), T, F)
  if (method=='linear') {
    coefficients <- linear.reg(formula, release, n, intercept)
  }else{
    coefficients <-NULL
    warning("the method you gave me is craycray")
  }
  release <- list(name='regression', 
                  n=n, 
                  formula=formula, 
                  coefficients=coefficients)
  return(release)
}
