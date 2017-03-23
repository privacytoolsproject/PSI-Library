#' Release differentially private most common values.
#'
#' @param x A vector of the data
#' @param epsilon Epsilon value for differential privacy
#' @param delta Delta value for differential privacy
#' @param k top-k elements to return
#' @param symb Symbol to return on failure (default : NA)
#' @return A vector whose values are k most common values or a failure return symbol.
#' @author Victor Balcer
#'
#' Contains functions for privately releasing the top-k heavy hitters, that is,
#'  the k most common unique values.  For example, k=1 releases the mode.
#'
#' @examples
#' n <- 10000
#' range <- c(0,20)
#' x <- rbinom(n, size=max(range), prob=0.7)
#' heavyhitters.release(x=x, epsilon=.1, delta=.0000001)

heavyhitters.release <- function(x, epsilon, delta, k=1, symb=NA){
  
  hist <- table(factor(x, exclude=c()), useNA="ifany") 

  # TODO: possible violation of privacy
  if(k > length(hist) - 1)
  {
    return(symb)
  }

  # TODO: should be replaced with a partial sorting algorithm
  hist <- sort(hist, decreasing=TRUE)[1:(k+1)]

  gap <- hist[k] - hist[k+1] + rlaplace(n=1, sensitivity=2, epsilon=epsilon)

  testfailure <- gap < -2 / epsilon * log(delta)

  if(testfailure){
    return(symb)
  }else{
    return(names(hist[1:k]))
  }

}

  
#' Get the accuracy of heavyhitters statistic for a given value of epsilon
#'
#' @param epsilon Epsilon value for differential privacy
#' @param delta Delta value for differential privacy
#' @param beta The true value is within the accuracy range with
#'    probability 1-beta
#' @param gap 
#' @return The accuracy guaranteed by the given epsilon
#' @author Victor Balcer
#'
#' The accuracy is interpreted as follows: The alpha value returned means that with 
#'   probability 1 - beta, simultaneously for all t with min <= t <= max, the algorithm's 
#'   estimate of the count in [min, t] is within alpha of the true value.
  
heavyhitters.getAccuracy <- function(epsilon, delta, beta=0.05, gap) {
  accuracy <- exp(-eps * gap / 2) / del
  return(accuracy)
}


#' Get the epsilon value necessary to guarantee a desired level of accuracy of a heavyhitters release
#'
#' @param delta Delta value for DP
#' @param beta the true value is within the accuracy range (alpha)
#    with probability 1-beta
#' @param gap something here
#' @return The epsilon value necessary to gaurantee the given accuracy
#' @author Victor Balcer

heavyhitters.getParameters <- function(delta, beta=0.05, gap) {
  epsilon <- -2 * log(beta * del) / gap
  return(epsilon)
}
