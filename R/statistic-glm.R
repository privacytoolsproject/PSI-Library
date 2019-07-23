#' Accuracy of the differentially private GLM
#'
#' Function to fin dthe accuracy guarantee of a GLM release at a given epsilon
#'
#' @param epsilon Numeric representing the epsilon privacy parameter. Should be 
#'    of length one and should be between zero and one.
#' @param n Integer specifying the number of observations.
#' @param alpha Numeric specifying the statistical significance level.
#'
#' @return Accuracy guarantee for GLM release
#' @rdname glm.getAccuracy
#' @export

glm.getAccuracy <- function(epsilon, n, alpha) {
    accuracy <- log(1 / alpha) / (n * epsilon)
    return(accuracy)
}


#' Privacy parameters for GLM
#'
#' Function to find the epsilon value necessary to meet a desired level of 
#' accuracy for a GLM release.
#'
#' @param accuracy Numeric giving the accuracy needed to guarantee.
#' @param n Integer specifying the number of observations.
#' @param alpha Numeric specifying the statistical significance level.
#' 
#' @return The scalar epsilon necessary to guarantee the needed accuracy.
#' @rdname glm.getParameters
#' @export

glm.getParameters <- function(accuracy, n, alpha) {
    epsilon <- log(1 / alpha) / (n * accuracy)
    return(epsilon)
}

#' Summary statistics for differentially private GLM via the bootstrap
#'
#' @param release Numeric matrix with differentially private estimates for each bootstrap sample
#' @param n Integer indicating number of observations
#' @param model Character indicating model form
#' @param alpha Numeric proportion of vector to be trimmed, specifically the 
#'      least and greatest \code{alpha / 2} are trimmed
#' @return Data frame summary statistics, including estimates and standard errors

glm.postSummary <- function(release, n, model, alpha) {
    trimmed.release <- apply(release, 2, trimVector, alpha=alpha)
    estimate <- apply(trimmed.release, 2, mean)
    std.error <- apply(trimmed.release, 2, sd)
    lower <- apply(trimmed.release, 2, quantile, 0.025)
    upper <- apply(trimmed.release, 2, quantile, 0.975)
    dp.summary <- data.frame(estimate, std.error, upper, lower)
    names(dp.summary) <- c('Estimate', 'Std. Error', 'CI95 Lower', 'CI95 Upper')
    rownames(dp.summary) <- names(release)
    if (model == 'ols') {
        variance <- as.numeric(dp.summary[nrow(dp.summary), 'Estimate'])
        coefficients <- dp.summary[1:(nrow(dp.summary) - 1), ]
        variance <- variance
        return(list('coefficients' = coefficients, 'variance' = variance))
    } else {
        return(dp.summary)
    }
}


#' Differentially private objective function for Logistic regression
#'
#' @return List with the name and objective function
# old dp Logit code blocked out below. 
# dpLogit <- function() {
    # objective.logit <- function(theta, X, y, b, n) {
        # xb <- as.matrix(X) %*% as.matrix(theta)
        # p <- as.numeric(1 / (1 + exp(-1 * xb)))
        # noise <- (b %*% as.matrix(theta)) / n
        # llik <- sum(y * log(p) + ((1 - y) * log(1 - p))) / n
        # llik.noisy <- noise + llik
        # return(-llik.noisy)
    # }
    # return(list('name' = 'logit', 'objective' = objective.logit))
# }

#' Differentially private objective function for Probit regression
#'
#' @return List with the name and objective function

#' Differentially private objective function for Logistic regression
#'
#' @return List with the name and objective function
# new dp logit including regularization below
dpLogit <- function() {
    objective.logit <- function(theta, X, y, b, n, lambda) {
    	theta <- as.matrix(theta)
        xb <- as.matrix(X) %*% theta
        p <- as.numeric(1 / (1 + exp(-1 * xb)))
        noise <- (b %*% theta) / n
        regularizer <- (lambda/2)* t(theta)%*%theta
        llik <- sum(y * log(p) + ((1 - y) * log(1 - p))) / n
        llik.noisy <- noise + llik - regularizer # double check that subtracting regularizer here is correct (as opposed to adding)
        return(-llik.noisy)
    }
    return(list('name' = 'logit', 'objective' = objective.logit))
}

# ' Differentially private objective function for Probit regression
# '
# ' @return List with the name and objective function






dpProbit <- function() {
    objective.probit <- function(theta, X, y, b, n) {
        xb <- as.matrix(X) %*% as.matrix(theta)
        p <- pnorm(xb)
        noise <- (b %*% as.matrix(theta)) / n
        llik <- sum(y * log(p) + ((1 - y) * log(1 - p))) / n
        llik.noisy <- noise + llik
        return(-llik.noisy)
    }
    return(list('name' = 'probit', 'objective' = objective.probit))
}

#' Differentially private objective function for Poisson regression
#'
#' @return List with the name and objective function

dpPoisson <- function() {
    objective.poisson <- function(theta, X, y, b, n) {
        lp <- as.matrix(X) %*% as.matrix(theta)
        noise <- (b %*% as.matrix(theta)) / n
        llik <- sum((y * lp) - exp(lp)) / n
        llik.noisy <- noise + llik
        return(-llik.noisy)
    }
    return(list('name' = 'poisson', 'objective' = objective.poisson))
}

#' Differentially private objective function for linear regression
#'
#' @return List with the name and objective function

dpOLS <- function() {
    objective.ols <- function(theta, X, y, b, n, lambda) {
        s <- exp(theta[length(theta)])
        beta <- theta[1:(length(theta) - 1)]
        xb <- as.matrix(X) %*% as.matrix(beta) 
        noise <- (b %*% as.matrix(theta)) / n
        llik <- ((-n / 2) * log(2 * pi) - n * log(s) - (0.5 / s^2) * sum((y - xb)^2)) / n
        llik.noisy <- noise + llik
        return(-llik.noisy)
    }
    return(list('name' = 'ols', 'objective' = objective.ols))
}

#' Objective functions
#'
#' List of objective functions support by \code{dpGLM}

glmObjectives = list(
    'logit' = dpLogit,
    'probit' = dpProbit,
    'poisson' = dpPoisson,
    'ols' = dpOLS
)


#' Differentially private generalized linear models
#'
#' @import methods
#' @export dpGLM
#' @exportClass dpGLM
#'
#' @include mechanism.R
#' @include mechanism-objective.R
#'
#' @examples
#'
#' data(PUMS5extract10000)
#' n <- 10000
#' epsilon <- 0.5
#' rng <- matrix(c(0, 150000, 1, 16, 0, 1), ncol=2, byrow=TRUE)
#' form <- 'income ~ educ + sex'
#'
#' model <- dpGLM$new(mechanism='mechanismObjective', var.type='numeric',
#'                    n=n, rng=rng, epsilon=epsilon, formula=form, objective='ols')
#' model$release(PUMS5extract10000)
#' print(model$result)

dpGLM <- setRefClass(
    Class = 'dpGLM',
    contains = 'mechanismObjective'
)

dpGLM$methods(
    initialize = function(mechanism, var.type, n, rng, formula, objective, epsilon=NULL,
                          accuracy=NULL, impute.rng=NULL, n.boot=NULL, intercept=TRUE, 
                          alpha=0.05) {
        fn <- glmObjectives[[objective]]()
        .self$name <- fn$name
        .self$objective <- fn$objective
        .self$mechanism <- mechanism
        .self$var.type <- var.type
        .self$n <- n
        .self$rng <- checkrange(rng, var.type)
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- glm.getParameters(accuracy, n, alpha)
        } else {
            .self$epsilon <- epsilon
            .self$accuracy <- glm.getAccuracy(epsilon, n, alpha)
        }
        if (is.null(impute.rng)) {
            .self$impute.rng <- rng
        } else {
            .self$impute.rng <- checkImputationRange(impute.rng)
        }
        .self$formula <- formula
        .self$n.boot <- n.boot
        .self$intercept <- intercept
        .self$alpha <- alpha
})

dpGLM$methods(
    release = function(x, ...) {
        .self$result <- export(mechanism)$evaluate(x, .self$postProcess, ...)
})

dpGLM$methods(
    postProcess = function(out) {
        out$variable <- all.vars(as.formula(formula))
        out$accuracy <- accuracy
        out$epsilon <- epsilon
        if (!is.null(n.boot)) {
            out$summary <- glm.postSummary(out$release, n, name, alpha)
        }
        return(out)
})
