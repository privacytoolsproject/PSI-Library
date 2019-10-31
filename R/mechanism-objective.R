#' Objective perturbation mechanism
#'
#' @import methods
#' @export mechanismObjective
#' @exportClass mechanismObjective
#'
#' @include mechanism.R

mechanismObjective <- setRefClass(
    Class = 'mechanismObjective',
    contains = 'mechanism'
)

mechanismObjective$methods(
    evaluate = function(x, postFun, ...) {

        # subset data from formula
        cols <- all.vars(as.formula(.self$formula))
        x <- x[, cols]

        # censor & impute missing values
        #x <- censorData(x, .self$varType, .self$rng, .self$bins) # NEEDS fixed to support checkRange
        x <- fillMissing(x, .self$varType, imputeRng=.self$imputeRng)

        # extract X and y
        y <- x[, cols[1]]
        X <- x[, cols[2:length(cols)], drop=FALSE]
        xNames <- names(X)

        # scale inputs s.t. max Euclidean norm <= 1
        scaler <- mapMatrixUnit(X, p=2)
        X <- scaler$matrix

        # add intercept
        if (.self$intercept) {
            X <- cbind(1, X)
            xNames <- c('intercept', xNames)
        }

        # set start params, adjust for ols
        if (.self$name == 'ols') {
            startParams <- rep(0, ncol(X) + 1)
            xNames <- c(xNames, 'variance')
            yScaler <- mapMatrixUnit(y, p=2)
            y <- yScaler$matrix
            yMaxNorm <- yScaler$maxNorm
        } else {
            startParams <- rep(0, ncol(X))
            yMaxNorm <- NULL
        }
		
		# Set scalar c from [CMS11]
		if(.self$name == 'logit'){
			c <- .25
		}
		else{
			c <- .25 # This is wrong! It is only a placeholder so that the code returns something. Must calculate c for each loss function we use.
			}
		#Set regularization parameter lambda
		lambda <- 1 #.self$n/20 # What should the default value be?
		#Ensure lambda satisfies condition in Algorithm 2 of [CMS11]
		compare <- c/(.self$n*(exp(.self$epsilon/2)-1))
		if(lambda <= compare){
			lambda <- compare + .001
		}
		term <- c/n*lambda
		ep <- .self$epsilon-log(1+2*term+term^2)

		beta <- ep/2
        # fit
        if (is.null(.self$nBoot)) {
           # old noise draw based on gamma and laplace
           # bNorm <- dpNoise(n=1, scale=(2 / .self$epsilon), dist='gamma', shape=length(startParams))
           # b <- dpNoise(n=length(startParams), scale=(-.self$epsilon * bNorm), dist='laplace')
           
           # new noise draw based on exponential and uniform
            bNorm <- dpNoise(n=1, scale=(1/beta), dist='gamma', shape=1)
            randomVec <- dpNoise(n=length(startParams), scale=1, dist='gaussian')
            randomVecNorm <- sqrt(sum(randomVec^2))
            b <- randomVec*(bNorm/randomVecNorm)
            
            estimates <- optim(par=startParams, fn=.self$objective, X=X, y=y, b=b, n=n, lambda=lambda)$par
            release <- data.frame(scaleRelease(estimates, scaler$maxNorm, yMaxNorm))
            names(release) <- 'estimate'
            rownames(release) <- xNames
        } else {
            localEpsilon <- .self$epsilon / .self$nBoot
            release <- vector('list', .self$nBoot)
            for (i in 1:.self$nBoot) {
                index <- sample(1:.self$n, .self$n, replace=TRUE)
                xStar <- X[index, ]
                yStar <- y[index]
                bNorm <- dpNoise(n=1, scale=(2 / localEpsilon), dist='gamma', shape=length(startParams))
                b <- dpNoise(n=length(startParams), scale=(-localEpsilon * bNorm), dist='laplace')
                estimates <- optim(par=startParams, fn=.self$objective, X=xStar, y=yStar, b=b, n=n)$par
                release[[i]] <- scaleRelease(estimates, scaler$maxNorm, yMaxNorm)
            }
            release <- data.frame(do.call(rbind, release))
            names(release) <- xNames
        }

        # format output
        out <- list('release' = release)
        out <- postFun(out, ...)
        return(out)
})
