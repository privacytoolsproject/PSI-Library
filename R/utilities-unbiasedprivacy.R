# Utility functions for the udp package


#' Generate Data function
#'
#'Creates a dataset with specified intercept and coefficient. Outputs a continuous outcome \code{Y1} and a binary outcome \code{Y2}.
#'
#'Generates data as follows:
#' \deqn{X \sim N(0, \sqrt{7})}
#' \deqn{Y1 = \alpha + \beta X + N(0, \sqrt{\sigma^2_y})}
#' \deqn{Y2 = \frac{exp(\alpha + \beta X)}{exp(\alpha + \beta X) + 1}}
#' @param N Number of rows desired in the output dataset
#' @param alpha Intercept
#' @param beta Coefficient value for relationship between X and Y
#' @param var Variance of Y1
#' @param n_sims Number of times to replicate X (useful for creating sets of simulated data with same X)
#' @param seed Random seed
#'
#' @return 
#' \item{data.frame(Y1, Y2, X)}
#' 
#' @examples
#' generateData(N = 10000, alpha = 10, beta = 0.8, var = 10)
#' 
#' @export
#' 
generateData <- function(N, alpha, beta, var, n_sims = 1, seed = 1234){
  set.seed(1)
  X <- rnorm(N, 0, 7)
  set.seed(seed)
  
  # continuous outcome
  Y1 <- alpha + beta*X + rnorm(N, 0, sd = sqrt(var))
  
  # stack X if we want more than dataset returned
  if(n_sims > 1){
    X <- rep(X, n_sims)
    Y1 <- alpha + beta*X + rnorm(N*n_sims, 0, sd = sqrt(var))
  }

  
  # binary outcome
  pi <-exp(alpha + beta*X)/(exp(alpha + beta*X) + 1) 
  Y2 <- rbinom(N, size = 1, prob = pi)
  
  # exp outcome 
  Y3 <- exp(alpha + beta*X + rnorm(N*n_sims, 0, sd = sqrt(var)))
  
  return(data.frame(Y1, Y2, Y3, X))
}


#' Weighted OLS coefficient function
#'
#'Runs a specified weighted regression and returns the coefficient of interest
#'
#' @param data Dataset
#' @param w Weights to use in the regression
#' @param form Formula of desired regression model
#' @param coef Coefficient of interest
#'
#' @return Weighted coefficient for specified coefficient
#' 
#' @examples
#' \dontrun{coefFn(dat, w = w, form = as.formula(Y ~ X), coef = 'X')}
#' 
#' @export
#' 
coefFn<- function(data, w, form, coef){
  data[, 'w'] <- w
  coef(lm(form, data = data, w = w))[coef]
}

#' Weighted logistic regression coefficient function
#'
#'Runs a specified weighted logistic regression and returns the coefficient of interest
#'
#' @param data Dataset
#' @param w Weights to use in the regression
#' @param form Formula of desired regression model
#' @param coef Coefficient of interest
#'
#' @return Weighted coefficient for specified coefficient
#' 
#' @examples
#' \dontrun{coefFn(dat, w = w, form = as.formula(Y ~ X), coef = 'X')}
#' 
#' @export
#' 
logitCoefFn <- function(data, w, form, coef){
  data[, 'w'] <- w
  coef(glm(form, data = data, weights = w, family = 'binomial'))[coef]
}


#' Calculate lambda value for any amount of censoring
#'
#' Calculates the SE of theta to return the lambda that will provide the specified amount of censoring in expectation
#'
#' @param true_theta True QOI 
#' @param X X value of dataset
#' @param P number of partitions
#' @param y_var variance of error term in Y
#' @param alpha Desired amount of censoring
#' @param upper Whether we have upper or lower censoring
#'
#' @return Lambda
#' 
#' @export
#' 
trueLambdaCalc <- function(true_theta, X, P, y_var, alpha, upper = T){
  se <- sqrt(y_var * solve(t(X)%*%X) * P)
  
  if(upper){
    lambda <- abs(qnorm(1 - alpha, mean = true_theta, sd = se))
  }else{
    lambda <- abs(qnorm(alpha, mean = true_theta, sd = se))
  }
  
  return(lambda)
}


#' Generate data and run UDP algorithm for any parameter set
#'
#' @export
#'
udpSim <- function(param_row, save_path){
 # TODO remove browser comments
 # browser()
  pr <- param_row
  dat <- generateData(pr$N, pr$intercept, pr$beta, pr$y_var, seed = pr$seed) # generate new dataset
  
  # Calculate lambda based on specified value of alpha (prop. censoring) and OLS SE
  l <- trueLambdaCalc(pr$beta, dat$X, pr$P, pr$y_var, pr$alpha)
  #true_sigma <- se * sqrt(pr$P)
  
  # QOI
  form <- as.formula(Y1~X)
  coef <- 'X'
  
  sim <- algorithmUDP(data = dat, statistic = coefFn, B = pr$R, n = pr$b, P = pr$P, lambda = l, lambda_var = 0.025, delta = 0.01,
                      epsilon = pr$e, epsilon_alpha = pr$e_alpha, parallelize = F, censoring_cutoff = 0.9,
                      bias_cutoff = 0.1, form = form, coef = coef)
  sim$params <- param_row
  
  fname <- paste0(save_path, '/sim_', pr$seed, '.Rdata')
  save(sim, file = fname)
  return(sim)
}
