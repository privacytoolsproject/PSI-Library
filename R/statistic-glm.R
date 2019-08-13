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
#' @rdname glmGetAccuracy
#' @export

glmGetAccuracy <- function(epsilon, n, alpha) {
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
#' @rdname glmGetParameters
#' @export

glmGetParameters <- function(accuracy, n, alpha) {
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

glmPostSummary <- function(release, n, model, alpha) {
    trimmedRelease <- apply(release, 2, trimVector, alpha=alpha)
    estimate <- apply(trimmedRelease, 2, mean)
    stdError <- apply(trimmedRelease, 2, sd)
    lower <- apply(trimmedRelease, 2, quantile, 0.025)
    upper <- apply(trimmedRelease, 2, quantile, 0.975)
    dpSummary <- data.frame(estimate, stdError, upper, lower)
    names(dpSummary) <- c('Estimate', 'Std. Error', 'CI95 Lower', 'CI95 Upper')
    rownames(dpSummary) <- names(release)
    if (model == 'ols') {
        variance <- as.numeric(dpSummary[nrow(dpSummary), 'Estimate'])
        coefficients <- dpSummary[1:(nrow(dpSummary) - 1), ]
        variance <- variance
        return(list('coefficients' = coefficients, 'variance' = variance))
    } else {
        return(dpSummary)
    }
}


#' Differentially private objective function for Logistic regression
#'
#' @return List with the name and objective function
# old dp Logit code blocked out below. 
# dpLogit <- function() {
    # objectiveLogit <- function(theta, X, y, b, n) {
        # xb <- as.matrix(X) %*% as.matrix(theta)
        # p <- as.numeric(1 / (1 + exp(-1 * xb)))
        # noise <- (b %*% as.matrix(theta)) / n
        # llik <- sum(y * log(p) + ((1 - y) * log(1 - p))) / n
        # llikNoisy <- noise + llik
        # return(-llikNoisy)
    # }
    # return(list('name' = 'logit', 'objective' = objectiveLogit))
# }

#' Differentially private objective function for Probit regression
#'
#' @return List with the name and objective function

#' Differentially private objective function for Logistic regression
#'
#' @return List with the name and objective function
# new dp logit including regularization below
dpLogit <- function() {
    objectiveLogit <- function(theta, X, y, b, n, lambda) {
    	theta <- as.matrix(theta)
        xb <- as.matrix(X) %*% theta
        p <- as.numeric(1 / (1 + exp(-1 * xb)))
        noise <- (b %*% theta) / n
        regularizer <- (lambda/2)* t(theta)%*%theta
        llik <- sum(y * log(p) + ((1 - y) * log(1 - p))) / n
        llikNoisy <- noise + llik - regularizer # double check that subtracting regularizer here is correct (as opposed to adding)
        return(-llikNoisy)
    }
    return(list('name' = 'logit', 'objective' = objectiveLogit))
}

# ' Differentially private objective function for Probit regression
# '
# ' @return List with the name and objective function






dpProbit <- function() {
    objectiveProbit <- function(theta, X, y, b, n) {
        xb <- as.matrix(X) %*% as.matrix(theta)
        p <- pnorm(xb)
        noise <- (b %*% as.matrix(theta)) / n
        llik <- sum(y * log(p) + ((1 - y) * log(1 - p))) / n
        llikNoisy <- noise + llik
        return(-llikNoisy)
    }
    return(list('name' = 'probit', 'objective' = objectiveProbit))
}

#' Differentially private objective function for Poisson regression
#'
#' @return List with the name and objective function

dpPoisson <- function() {
    objectivePoisson <- function(theta, X, y, b, n) {
        lp <- as.matrix(X) %*% as.matrix(theta)
        noise <- (b %*% as.matrix(theta)) / n
        llik <- sum((y * lp) - exp(lp)) / n
        llikNoisy <- noise + llik
        return(-llikNoisy)
    }
    return(list('name' = 'poisson', 'objective' = objectivePoisson))
}

#' Differentially private objective function for linear regression
#'
#' @return List with the name and objective function

dpOLS <- function() {
    objectiveOLS <- function(theta, X, y, b, n, lambda) {
        s <- exp(theta[length(theta)])
        beta <- theta[1:(length(theta) - 1)]
        xb <- as.matrix(X) %*% as.matrix(beta) 
        noise <- (b %*% as.matrix(theta)) / n
        llik <- ((-n / 2) * log(2 * pi) - n * log(s) - (0.5 / s^2) * sum((y - xb)^2)) / n
        llikNoisy <- noise + llik
        return(-llikNoisy)
    }
    return(list('name' = 'ols', 'objective' = objectiveOLS))
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
#' model <- dpGLM$new(mechanism='mechanismObjective', varType='numeric',
#'                    n=n, rng=rng, epsilon=epsilon, formula=form, objective='ols')
#' model$release(PUMS5extract10000)
#' print(model$result)

dpGLM <- setRefClass(
    Class = 'dpGLM',
    contains = 'mechanismObjective'
)

dpGLM$methods(
    initialize = function(mechanism, varType, n, rng=NULL, formula, objective, epsilon=NULL,
                          accuracy=NULL, imputeRng=NULL, nBoot=NULL, intercept=TRUE, 
                          alpha=0.05) {
        fn <- glmObjectives[[objective]]()
        .self$name <- fn$name
        .self$objective <- fn$objective
        .self$mechanism <- mechanism
        .self$varType <- varType
        .self$n <- checkNValidity(n)
        .self$rng <- checkRange(rng, varType)
        if (is.null(epsilon)) {
            .self$accuracy <- accuracy
            .self$epsilon <- glmGetParameters(accuracy, n, alpha)
        } else {
            .self$epsilon <- epsilon
            .self$accuracy <- glmGetAccuracy(epsilon, n, alpha)
        }
        if (is.null(imputeRng)) {
            .self$imputeRng <- rng
        } else {
            .self$imputeRng <- checkImputationRange(imputeRng)
        }
        .self$formula <- formula
        .self$nBoot <- nBoot
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
        if (!is.null(nBoot)) {
            out$summary <- glmPostSummary(out$release, n, name, alpha)
        }
        return(out)
})
