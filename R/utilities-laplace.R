#' Get accuracy for Laplace statistics
#' 
#' Function to find the accuracy guarantee of a statistic release at a given epsilon 
#' value.
#' 
#' @param sensitivity the sensitivity of the statistic
#' @param epsilon A numeric vector representing the epsilon privacy parameter.
#'    Should be of length one and should be between zero and one.
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @return Accuracy guarantee for statistic release given epsilon.

laplaceGetAccuracy <- function(sensitivity, epsilon, alpha=0.05) {
  accuracy <- log(1 / alpha) * (sensitivity / epsilon)
  return(accuracy)
}


#' Get epsilon for Laplace statistics
#' 
#' Function to find the epsilon value necessary to meet a desired level of 
#' accuracy for a statistic release.
#' 
#' @param sensitivity the sensitivity of the statistic
#' @param accuracy A numeric vector representing the accuracy needed to 
#'    guarantee (percent).
#' @param alpha A numeric vector specifying the statistical significance level.
#' 
#' @return The scalar epsilon necessary to guarantee the needed accuracy.

laplaceGetEpsilon <- function(sensitivity, accuracy, alpha=0.05) {
  epsilon <- log(1 / alpha) * (sensitivity / accuracy)
  return(epsilon)
}

