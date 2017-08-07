#' Differentially private objective function for Logistic regression
#' 
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Formula for the Logistic regression model
#' @param intercept Logical indicating whether the intercept should be added to the model
#' @export
dp.logit <- function(n, epsilon, formula, intercept) {
    objective.logit <- function(theta, X, y, b, n) {
        xb <- as.matrix(X) %*% as.matrix(theta)
        p <- as.numeric(1 / (1 + exp(-1 * xb)))
        noise <- (b %*% as.matrix(theta)) / n
        llik <- sum(y * log(p) + ((1 - y) * log(1 - p))) / n
        llik.noisy <- noise + llik
        return(-llik.noisy)
    }
    return(list('name' = 'glm',
                'name' = 'logit',
                'objective' = objective.logit,
                'n' = n,
                'epsilon' = epsilon,
                'formula' = formula,
                'intercept' = intercept))
}


#' Differentially private objective function for Probit regression
#' 
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Formula for the Logistic regression model
#' @param intercept Logical indicating whether the intercept should be added to the model
#' @export
dp.probit <- function(n, epsilon, formula, intercept) {
    objective.probit <- function(theta, X, y, b, n) {
        xb <- as.matrix(X) %*% as.matrix(theta)
        p <- pnorm(xb)
        noise <- (b %*% as.matrix(theta)) / n
        llik <- sum(y * log(p) + ((1 - y) * log(1 - p))) / n
        llik.noisy <- noise + llik
        return(-llik.noisy)
    }
    return(list('name' = 'glm',
                'model' = 'probit',
                'objective' = objective.probit,
                'n' = n,
                'epsilon' = epsilon,
                'formula' = formula,
                'intercept' = intercept))
}


#' Differentially private objective function for Poisson regression
#' 
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Formula for the Poisson regression model
#' @param intercept Logical indicating whether the intercept should be added to the model
#' @export
dp.poisson <- function(n, epsilon, formula, intercept) {
    objective.poisson <- function(theta, X, y, b, n) {
        lp <- X %*% as.matrix(theta)
        noise <- (b %*% as.matrix(theta)) / n
        llik <- sum((y * lp) - exp(lp)) / n
        llik.noisy <- noise + llik
        return(-llik.noisy)
    }
    return(list('name' = 'ols',
                'model' = 'poisson',
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
#' @param formula Formula for the linear regression model
#' @param intercept Logical indicating whether the intercept should be added to the model
#' @export
dp.ols <- function(n, epsilon, formula, intercept) {
    objective.ols <- function(theta, X, y, b, n) {
        s <- exp(theta[length(theta)])
        beta <- theta[1:(length(theta) - 1)]
        xb <- as.matrix(X) %*% as.matrix(beta) 
        noise <- (b %*% as.matrix(theta)) / n
        llik <- ((-n / 2) * log(2 * pi) - n * log(s) - (0.5 / s^2) * sum((y - xb)^2)) / n
        llik.noisy <- noise + llik
        return(-llik.noisy)
    }
    return(list('name' = 'glm',
                'model' = 'ols',
                'objective' = objective.ols,
                'n' = n,
                'epsilon' = epsilon,
                'formula' = formula,
                'intercept' = intercept))
}


#' Release a differentially private vector of parameter estimates for a generalized linear model
#'
#' @param x Dataframe with data needed to estimate the model
#' @param n Integer indicating number of observations
#' @param epsilon Numeric epsilon parameter for differential privacy
#' @param formula Model formula
#' @param objective Perturbed objective function
#' @param n.boot Number of bootstrap samples on which to evaluate the 
#'    private parameter estimates, each at \code{epsilon / n.boot} privacy cost
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
#' logit.private <- glm.release(data, nrow(data), epsilon=0.5, 
#'    formula=form, objective=dp.logit)$release
#' form2 <- as.formula('y ~ x1 + x2 + x3')  # add a factor variable
#' logit.private2 <- glm.release(data, nrow(data), epsilon=0.5, 
#'    formula=form2, objective=dp.logit)$release
#' @export
glm.release <- function(x, n, epsilon, formula, objective, n.boot=NULL, intercept=TRUE) {
    postlist <- NULL
    if (!is.null(n.boot)) {
        postlist <- list('summary' = 'postSummary')
    }
    release <- mechanism.objective(fun=objective, x=x, n=n, epsilon=epsilon, formula=formula, 
                                   n.boot=n.boot, intercept=intercept, postlist=postlist)
    return(release)
}


#' Summary statistics for differentially private GLM via the bootstrap
#'
#' @param release Numeric matrix with differentially private estimates for each bootstrap sample
#' @param alpha Numeric proportion of vector to be trimmed, specifically the 
#'      least and greatest \code{alpha / 2} are trimmed
#' @return Data frame summary statistics, including estimates and standard errors

glm.postSummary <- function(release, n, model, alpha=0.10) {
    trimmed.release <- apply(release, 2, trimVector, alpha=alpha)
    estimate <- apply(trimmed.release, 2, mean)
    std.error <- apply(trimmed.release, 2, sd)
    t.values <- estimate / std.error
    p.values <- 2 * pt(abs(t.values), df=(n - length(t.values)), lower.tail=FALSE)
    dp.summary <- data.frame(estimate, std.error, t.values, p.values)
    names(dp.summary) <- c('Estimate', 'Std. Error', 'Statistic', 'p-value')
    rownames(dp.summary) <- names(release)
    out.summary <- list('coefficients' = dp.summary)
    if (model == 'ols') {
        variance <- as.numeric(dp.summary[nrow(dp.summary), 'Estimate'])
        out.summary$coefficients <- dp.summary[1:(nrow(dp.summary) - 1), ]
        out.summary$variance <- variance
    }
    return(out.summary)
}
