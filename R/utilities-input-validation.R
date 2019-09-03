#' Check validity of n
#' 
#' n should always be a positive whole number, check the user's input
#' 
#' @param n the input n from te user
#' 
#' @return n, if n is a positive whole number

checkN <- function(n) {
  if ((n > 0) & (n%%1 == 0)) {
    return(n)
  } else {
    stop("n must be a positive whole number")
  }
}

