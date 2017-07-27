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


#' Release a differentially private vector of parameter estimates for a Poisson regression model
#'
#' @param x Dataframe with data needed to estimate the model
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Formula for the logistic regression model
#' @param intercept Logical indicating whether the intercept should be added to the model, default TRUE
#' 
#' @examples
#' set.seed(84753)
#' n <- 1000
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rbinom(n, 1, prob=0.3)
#' X <- cbind(1, x1, x2, x3, x4)
#' z <- -0.3 - 1.6 * x1 + 0.3 * x2 + 0.01 * x3 - 0.73 * x4
#' lambda <- exp(z)
#' y <- rpois(n, lambda=lambda)
#' data <- data.frame(y, x1, x2, x3, x4)
#' form <- as.formula('y ~ x1 + x2 + x3 + x4')
#' poisson.true <- glm(form, data=data, family=poisson(link=log))$coef
#' poisson.private <- poisson.release(data, nrow(data), epsilon=0.5, formula=form)$release

poisson.release <- function(x, n, epsilon, formula, n.boot=NULL, intercept=TRUE) {
    release <- mechanism.objective(fun=dp.poisson, x=x, n=n, epsilon=epsilon, formula=formula, n.boot=n.boot, intercept=intercept)
    return(release)
}
