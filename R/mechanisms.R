#' Base mechanism class
#'
#' @import methods
#' @export mechanism
#' @exportClass mechanism
#'
#' @field mechanism Name of the mechanism
#' @field name Name of the statistic
#' @field var.type Variable type
#' @field var.type.orig Variable type at instantiation
#' @field n Number of observations
#' @field epsilon Differential privacy parameter
#' @field delta Differential privacy parameter
#' @field rng A priori estimate of the variable range
#' @field result List with statistical output
#' @field alpha Level of statistical signficance
#' @field accuracy Accuracy guarantee of the estimate
#' @field bins Bins
#' @field n.bins Number of bins
#' @field k Number of bins desired for the release
#' @field error Error
#' @field n.boot Number of bootstrap replications
#' @field boot.fun Function passed to the bootstrap mechanism
#' @field impute.rng The range from which to impute missing values
#' @field impute Logical, impute categorical types?
#' @field formula R formula for regression models
#' @field columns Vector of column names
#' @field intercept Logical, is the intercept included?
#' @field stability Logical, use stability histogram
#' @field objective Objective function for regression models
#' @field gran Granularity
#' @field percentiles Percentiles evaluated by binary tree
#' @field tree.data Binary tree attributes needed for efficient estimation

mechanism <- setRefClass(
    Class = 'mechanism',
    fields = list(
        mechanism = 'character',
        name = 'character',
        var.type = 'character',
        var.type.orig = 'character',
        n = 'numeric',
        epsilon = 'numeric',
        delta = 'numeric',
        rng = 'ANY',
        result = 'ANY',
        alpha = 'numeric',
        accuracy = 'numeric',
        bins = 'ANY',
        n.bins = 'ANY',
        k = 'numeric',
        error = 'numeric',
        n.boot = 'ANY',
        boot.fun = 'function',
        impute.rng = 'ANY',
        impute = 'logical',
        formula = 'ANY',
        columns = 'ANY',
        intercept = 'logical', 
        stability = 'logical',
        objective = 'function',
        gran = 'numeric',
        percentiles = 'ANY',
        tree.data = 'ANY'
))

mechanism$methods(
    getFields = function() {
        f <- names(getRefClass()$fields())
        out <- setNames(vector('list', length(f)), f)
        for (fd in f) {
            out[[fd]] <- .self[[fd]]
        }
        return(out)
})

mechanism$methods(
    getFunArgs = function(fun) {
        f <- .self$getFields()
        spec <- list()
        for (arg in names(f)) {
            if (arg %in% names(formals(fun))) {
                spec[[arg]] <- f[[arg]]
            }
        }
        return(spec)
})


#' Laplace mechanism
#'
#' @import methods
#' @export mechanismLaplace
#' @exportClass mechanismLaplace

mechanismLaplace <- setRefClass(
    Class = 'mechanismLaplace',
    contains = 'mechanism'
)

mechanismLaplace$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismLaplace$methods(
    evaluate = function(fun, x, sens, postFun, ...) {
        x <- censordata(x, .self$var.type, .self$rng, .self$bins)
        if (.self$var.type %in% c('numeric', 'integer', 'logical')) {
            if (NCOL(x) > 1) {
                x <- fillMissing2d(x, .self$var.type, .self$impute.rng)
            } else {
                x <- fillMissing(x, .self$var.type, .self$impute.rng[1], .self$impute.rng[2])
            }
        } else {
            x <- fillMissing(x, .self$var.type, categories=.self$bins)
        }
        field.vals <- .self$getFunArgs(fun)
        ellipsis.vals <- getFuncArgs(list(...), fun)
        true.val <- do.call(fun, c(list(x=x), field.vals, ellipsis.vals))
        scale <- sens / .self$epsilon
        release <- true.val + dpNoise(n=length(true.val), scale=scale, dist='laplace')
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})


#' Exponential mechanism
#'
#' @import methods
#' @export mechanismExponential
#' @exportClass mechanismExponential

mechanismExponential <- setRefClass(
    Class = 'mechanismExponential',
    contains = 'mechanism'
)

mechanismExponential$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismExponential$methods(
    evaluate = function(fun, x, sens, postFun, ...) {
        x <- censordata(x, .self$var.type, levels=.self$bins)
        x <- fillMissing(x, .self$var.type, categories=.self$bins)
        field.vals <- .self$getFunArgs(fun)
        true.val <- do.call(fun, c(list(x=x), field.vals))
        quality <- true.val - max(true.val)
        probs <- ifelse(true.val == 0, 0, exp((.self$epsilon * quality) / (2 * sens)))
        gap <- as.numeric(true.val[.self$k] - true.val[.self$k + 1])
        if (gap < (-2 / epsilon * log(delta))) {
            out <- list('release' = NULL)
        } else {
            release <- sample(names(true.val), size=.self$k, prob=probs)
            out <- list('release' = release)
            out <- postFun(out, gap)
        }
        return(out)
})


#' Gaussian mechanism
#'
#' @import methods
#' @export mechanismGaussian
#' @exportClass mechanismGaussian

mechanismGaussian <- setRefClass(
    Class = 'mechanismGaussian',
    contains = 'mechanism'
)

mechanismGaussian$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismGaussian$methods(
    evaluate = function(fun, x, sens, postFun) {
        x <- censordata(x, .self$var.type, .self$rng)
        if (.self$var.type %in% c('numeric', 'integer', 'logical')) {
            if (NCOL(x) > 1) {
                x <- fillMissing2d(x, .self$var.type, .self$impute.rng)
            } else {
                x <- fillMissing(x, .self$var.type, .self$impute.rng[1], .self$impute.rng[2])
            }
        } else {
            x <- fillMissing(x, .self$var.type, categories=.self$bins)
        }
        field.vals <- .self$getFunArgs(fun)
        true.val <- do.call(fun, c(list(x=x), field.vals))
        scale <- sens * sqrt(2 * log(1.25 / .self$delta)) / .self$epsilon
        release <- true.val + dpNoise(n=length(true.val), scale=scale, dist='gaussian')
        out <- list('release' = release)
        out <- postFun(out)
        return(out)
})


#' Bootstrap replication for a function
#'
#' @param x Vector
#' @param n Number of observations
#' @param sensitivity Sensitivity of the function
#' @param epsilon Numeric differential privacy parameter
#' @param fun Function to evaluate
#' @return Value of the function applied to one bootstrap sample
#' @import stats
#' @export

bootstrap.replication <- function(x, n, sensitivity, epsilon, fun) {
    partition <- rmultinom(n=1, size=n, prob=rep(1 / n, n))
    max.appearances <- max(partition)
    probs <- sapply(1:max.appearances, dbinom, size=n, prob=(1 / n))
    stat.partitions <- vector('list', max.appearances)
    for (i in 1:max.appearances) {
        variance.i <- (i * probs[i] * (sensitivity^2)) / (2 * epsilon)
        stat.i <- fun(x[partition == i])
        noise.i <- dpNoise(n=length(stat.i), scale=sqrt(variance.i), dist='gaussian')
        stat.partitions[[i]] <- i * stat.i + noise.i
    }
    stat.out <- do.call(rbind, stat.partitions)
    return(apply(stat.out, 2, sum))
}


#' Bootstrap mechanism
#'
#' @import methods
#' @export mechanismBootstrap
#' @exportClass mechanismBootstrap

mechanismBootstrap <- setRefClass(
    Class = 'mechanismBootstrap',
    contains = 'mechanism'
)

mechanismBootstrap$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismBootstrap$methods(
    bootStatEval = function(xi) {
        field.vals <- .self$getFunArgs(boot.fun)
        stat <- do.call(boot.fun, c(list(xi=xi), field.vals))
        return(stat)
})

mechanismBootstrap$methods(
    bootSE = function(release, n.boot, sens) {
        se <- sd(release)
        c.alpha <- qchisq(0.01, df=(n.boot - 1))
        conservative <- sqrt(max(c(se^2 - (c.alpha * sens^2 * n.boot) / (2 * epsilon * (n.boot - 1)), 0)))
        naive <- sqrt(max(c(se^2 - (sens^2 * n.boot) / (2 * epsilon), 0)))
        return(list('sd' = se,
                    'conservative' = conservative,
                    'naive' = naive))
})

mechanismBootstrap$methods(
    evaluate = function(fun, x, sens, postFun) {
        x <- censordata(x, .self$var.type, .self$rng)
        x <- fillMissing(x, .self$var.type, .self$impute.rng[0], .self$impute.rng[1])
        epsilon.part <- epsilon / .self$n.boot
        release <- replicate(.self$n.boot, bootstrap.replication(x, n, sens, epsilon.part, fun=.self$bootStatEval))
        std.error <- .self$bootSE(release, .self$n.boot, sens)
        out <- list('release' = release, 'std.error' = std.error)
        out <- postFun(out)
        return(out)
})


#' Objective perturbation mechanism
#'
#' @import methods
#' @export mechanismObjective
#' @exportClass mechanismObjective

mechanismObjective <- setRefClass(
    Class = 'mechanismObjective',
    contains = 'mechanism'
)

mechanismObjective$methods(
    getFunArgs = function(fun) {
        callSuper(fun)
})

mechanismObjective$methods(
    evaluate = function(x, postFun, ...) {

        # subset data from formula
        cols <- all.vars(as.formula(.self$formula))
        x <- x[, cols]

        # censor & impute missing values
        x <- censordata(x, .self$var.type, .self$rng, .self$bins)
        x <- fillMissing2d(x, .self$var.type, .self$impute.rng)

        # extract X and y
        y <- x[, cols[1]]
        X <- x[, cols[2:length(cols)], drop=FALSE]
        X.names <- names(X)

        # scale inputs s.t. max Euclidean norm <= 1
        scaler <- mapMatrixUnit(X, p=2)
        X <- scaler$matrix

        # add intercept
        if (.self$intercept) {
            X <- cbind(1, X)
            X.names <- c('intercept', X.names)
        }

        # set start params, adjust for ols
        if (.self$name == 'ols') {
            start.params <- rep(0, ncol(X) + 1)
            X.names <- c(X.names, 'variance')
            y.scaler <- mapMatrixUnit(y, p=2)
            y <- y.scaler$matrix
            y.max.norm <- y.scaler$max.norm
        } else {
            start.params <- rep(0, ncol(X))
            y.max.norm <- NULL
        }

        # fit
        if (is.null(.self$n.boot)) {
            b.norm <- dpNoise(n=1, scale=(2 / .self$epsilon), dist='gamma', shape=length(start.params))
            b <- dpNoise(n=length(start.params), scale=(-.self$epsilon * b.norm), dist='laplace')
            estimates <- optim(par=start.params, fn=.self$objective, X=X, y=y, b=b, n=n)$par
            release <- data.frame(scaleRelease(estimates, scaler$max.norm, y.max.norm))
            names(release) <- 'estimate'
            rownames(release) <- X.names
        } else {
            local.epsilon <- .self$epsilon / .self$n.boot
            release <- vector('list', .self$n.boot)
            for (i in 1:.self$n.boot) {
                index <- sample(1:.self$n, .self$n, replace=TRUE)
                X.star <- X[index, ]
                y.star <- y[index]
                b.norm <- dpNoise(n=1, scale=(2 / local.epsilon), dist='gamma', shape=length(start.params))
                b <- dpNoise(n=length(start.params), scale=(-local.epsilon * b.norm), dist='laplace')
                estimates <- optim(par=start.params, fn=.self$objective, X=X.star, y=y.star, b=b, n=n)$par
                release[[i]] <- scaleRelease(estimates, scaler$max.norm, y.max.norm)
            }
            release <- data.frame(do.call(rbind, release))
            names(release) <- X.names
        }

        # format output
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})
