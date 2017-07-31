#' Differentially private objective function for Logistic regression
#' 
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Formula for the logistic regression model
#' @param intercept Logical indicating whether the intercept should be added to the model

dp.logit <- function(n, epsilon, formula, intercept) {
    objective.logit <- function(theta, X, y, b, n) {
        p <- as.numeric(1 / (1 + exp(-1 * as.matrix(X) %*% as.matrix(theta))))
        llik <- ((b %*% as.matrix(theta)) / n) + (sum(y * log(p) + ((1 - y) * log(1 - p))) / n)
        return(-llik)
    }
    return(list('name' = 'logit',
                'objective' = objective.logit,
                'n' = n,
                'epsilon' = epsilon,
                'formula' = formula,
                'intercept' = intercept))
}


#' Differentially private objective function for Poisson regression
#' 
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Formula for the logistic regression model
#' @param intercept Logical indicating whether the intercept should be added to the model

dp.poisson <- function(n, epsilon, formula, intercept) {
    objective.poisson <- function(theta, X, y, b, n) {
        lp <- X %*% as.matrix(theta)
        llik <- ((b %*% as.matrix(theta)) / n) + (sum((y * lp) - exp(lp)))
        return(-llik)
    }
    return(list('name' = 'logit',
                'objective' = objective.poisson,
                'n' = n,
                'epsilon' = epsilon,
                'formula' = formula,
                'intercept' = intercept))
}


#' Differentially private objective function for ordinary least squares regression
#' 
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Formula for the logistic regression model
#' @param intercept Logical indicating whether the intercept should be added to the model

dp.ols <- function(n, epsilon, formula, intercept) {
  objective.ols <- function(theta, X, y, b, n) {
    s <- exp( theta[length(theta)] )  # Constrain variance to be positive
    beta <- theta[1:(length(theta)-1)]      # Separate coefficients on covariates from variance
    xb <- as.matrix(X) %*% as.matrix(beta)  
    llik <- ((b %*% as.matrix(beta)) / n) + ((-n/2)*log(2*pi)-n*log(s)-(0.5/s^2)*sum((y-xb)^2))
    return(-llik)
  }
  return(list('name' = 'ols',
              'objective' = objective.ols,
              'n' = n,
              'epsilon' = epsilon,
              'formula' = formula,
              'intercept' = intercept))
}


#' Release a differentially private vector of parameter estimates for a logistic regression model
#'
#' @param x Dataframe with data needed to estimate the model
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Formula for the logistic regression model
#' @param intercept Logical indicating whether the intercept should be added to the model, default TRUE
#' 
#' @examples
#' set.seed(847523)
#' n <- 1000
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- sample(c('a', 'b', 'c', 'd', 'e'), replace=TRUE, size=n)
#' z <- -0.3 - 1.6 * x1 + 0.3 * x2
#' p <- 1 / (1 + exp(-z))
#' y <- rbinom(n, 1, p)
#' data <- data.frame(cbind(y, x1, x2))
#' form <- as.formula('y ~ x1 + x2')
#' logit.true <- glm(form, data=data, family='binomial')$coef
#' logit.private <- glm.release(data, nrow(data), epsilon=0.5, formula=form, objective=dp.logit)$release
#' form2 <- as.formula('y ~ x1 + x2 + x3')  # add a factor variable
#' logit.private2 <- glm.releases(data, nrow(data), epsilon=0.5, formula=form2, objective=dp.logit)$release

glm.release <- function(x, n, epsilon, formula, objective, n.boot=NULL, intercept=TRUE) {
    release <- mechanism.objective(fun=objective, x=x, n=n, epsilon=epsilon, formula=formula, n.boot=n.boot, intercept=intercept)
    return(release)
}
