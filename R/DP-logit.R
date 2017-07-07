#' Differentially private objective function for logistic regression
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
#' z <- -0.3 - 1.6 * x1 + 0.3 * x2
#' p <- 1 / (1 + exp(-z))
#' y <- rbinom(n, 1, p)
#' data <- data.frame(cbind(y, x1, x2))
#' form <- as.formula('y ~ x1 + x2')
#' logit.true <- glm(form, data=data, family='binomial')$coef
#' logit.private <- logit.release(data, nrow(data), epsilon=0.5, formula=form)$release

logit.release <- function(x, n, epsilon, formula, intercept=TRUE) {
    release <- mechanism.objective(fun=dp.logit, x=x, n=n, epsilon=epsilon, formula=formula, intercept=intercept)
    return(release)
}
