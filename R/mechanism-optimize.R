#' Optimization mechanism
#'
#' @import methods
#' @export mechanismOptimize
#' @exportClass mechanismOptimize
#'
#' @include mechanism.R
source('mechanism.R')

# Compute the multi-class cross entropy loss for each observation with mutually exclusive class membership
# expected: character vector of [batchSize]
# predicted: numeric matrix of [batchSize, numClasses], with column names the levels of expected
# returns: numeric matrix of [batchSize]
crossEntropy <- function(expected, predicted) {
  expected <- as.character(expected)

  # rows must sum to one
  predicted <- predicted / rowSums(predicted)

  # compute loss for each observation
  -sapply(1:length(expected), function(rowIdx)
    log(predicted[rowIdx, expected[rowIdx]]))
}


# Compute the binary log-likelihood loss for each observation with mutually exclusive level membership
# expected: numeric vector or factor of [batchSize], positive label is first in sort(expected)
# predicted: numeric vector of [batchSize]
# returns: numeric matrix of [batchSize]
binaryCrossEntropy <- function(expected, predicted) {
  predicted <- cbind(predicted, 1 - predicted)
  colnames(predicted) <- as.character(sort(unique(expected)))
  crossEntropy(expected, predicted)
}


# Compute the sum squared error
# expected: numeric vector of [batchSize]
# predicted: numeric vector of [batchSize]
sumSquaredError <- function(expected, predicted)
  sum((expected - predicted)^2)


# Compute the mean squared error
# expected: numeric vector of [batchSize]
# predicted: numeric vector of [batchSize]
meanSquaredError <- function(expected, predicted)
  sumSquaredError(expected, predicted) / length(expected)


costFunctions <- list(
  crossEntropy=crossEntropy,
  binaryCrossEntropy=binaryCrossEntropy,
  sumSquaredError=sumSquaredError,
  meanSquaredError=meanSquaredError
)

gaussianReleaseNoise <- function(size, sensitivity, epsilon, delta){
  # Source from 3.1 Differentially Private SGD Algorithm, Moments Accountant
  # https://arxiv.org/pdf/1607.00133.pdf
  sigma <- sqrt(2*log(1.25/delta)) / epsilon
  rnorm(n=size, mean=0, sd=sensitivity * sigma)
}


clip <- function(x, lower, upper){
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)
}


mechanismOptimize <- setRefClass(
    Class = 'mechanismOptimize',
    contains = 'mechanism',
    fields = list(
    predictors = 'character',
    forwardFun = 'function',
    theta = 'matrix',
    batchSize = 'ANY',
    stepSize = 'numeric',
    gradientFun = 'function',
    costFun = 'function',
    gradientNumericOffset = 'numeric',
    clippingInterval = 'numeric'
  ),
  methods=list(
    initialize = function(..., forwardFun, predictors, theta, clippingInterval, costFun,
    batchSize=NULL, stepSize=.01, gradientFun=NULL, gradientNumericOffset=.0001) {
      callSuper(...)
      .self$forwardFun <- forwardFun
      .self$predictors <- predictors
      .self$theta <- theta
      .self$clippingInterval <- clippingInterval
      .self$batchSize <- batchSize
      .self$stepSize <- stepSize

      if (is.character(costFun)) .self$costFun <- costFunctions[[costFun]]
      else {
        warning('Privacy is not guaranteed for a custom cost function.')
        .self$costFun <- costFun
      }

      if (is.null(gradientFun)) .self$gradientFun <- function(lossFun, stimulus, expected, theta) {
        original <- lossFun(expected, stimulus, theta)

        return(sapply(1:ncol(theta), function(y) {
          sapply(1:nrow(theta), function(x) {
            offset <- do.call(matrix, list(data=0, nrow=nrow(theta), ncol=ncol(theta)))
            offset[x, y] <- gradientNumericOffset
            print(offset)
            print(dim(offset))
            perturbed <- lossFun(expected, stimulus, theta + offset)
            (perturbed - original) / gradientNumericOffset
          })
        }, simplify="array"))
      }
      else {
        warning('Privacy is not guaranteed for a custom gradient function.')
        .self$gradientFun <- gradientFun
      }
      .self
    },

    # x: dataset
    evaluate = function(x, ...) {
      # x <- fillMissing(x, .self$var.type, rng=.self$rng, categories=.self$bins)
      # fun.args <- getFuncArgs(.self$forwardFun, inputList=list(...), inputObject=.self)
      # predictFun <- function(x, theta) do.call(.self$forwardFun, c(list(x), fun.args))

      # loss is the composition of cost and predict functions
      lossFun <- function(expected, stimulus, theta) .self$costFun(expected, .self$forwardFun(stimulus, theta))

      .self$batchSize <- if (is.null(.self$batchSize)) nrow(x) else .self$batchSize
      sensitivity <- 2 * .self$clippingInterval / .self$batchSize

      samplingRatio <- .self$batchSize / nrow(x)

      print('iterations')
      print(floor(epsilon / samplingRatio)*20)

      for (iteration in 1:floor(epsilon / samplingRatio)*20) {
        # print(iteration)
        batch <- if (.self$batchSize == length(x)) x else x[sample(1:nrow(x), .self$batchSize),]

        stimulus <- batch[,.self$predictors, drop=FALSE]
        expected <- batch[,.self$variable]

        # compute gradient
        grad <- .self$gradientFun(lossFun, stimulus, expected, .self$theta)

        # clip and collapse gradient
        grad <- clip(grad, -.self$clippingInterval, .self$clippingInterval)
        grad <- apply(grad, 2:3, mean)

        # compute update direction
        direction <- .self$iterate(grad, iteration)

        noise <- gaussianReleaseNoise(length(.self$theta), sensitivity, .self$epsilon * nrow(batch) / 2, .self$delta) / nrow(batch)

        .self$theta <- .self$theta + direction # + noise
      }
    },

    release = function(data) {
      x <- data[,c(.self$variable, .self$predictors)]
      .self$result <- .self$evaluate(x)
    },
    postProcess = function() {
      .self$results <- list(
        value=.self$theta
      )
    }
  )
)

dpOptimizerSGD <- setRefClass(
  Class = 'dpOptimizerSGD',
  contains = 'mechanismOptimize',
  methods = list(
    iterate = function(gradient, stepNumber) -.self$stepSize * gradient
  )
)

dpOptimizerAdam <- setRefClass(
  Class = 'dpOptimizerAdam',
  contains = 'mechanismOptimize',
  fields = list(
    gradientCache = 'numeric',
    gradientSquare = 'numeric',
    decayMoments = 'numeric'
  ),
  methods = list(
    initialize=function(..., gradCache=NULL, gradSquare=NULL, decayMoments=list(.9, .999), wedge=1e-8) {
      callSuper(...)
      .self$gradCache <- if(is.null(gradSquare)) matrix(0L, nrow=nrow(.self$theta), ncol=ncol(.self$theta)) else gradCache
      .self$gradSquare <- if(is.null(gradSquare)) matrix(0L, nrow=nrow(.self$theta), ncol=nrow(.self$theta)) else gradSquare

      .self$decay <- decay
      .self$wedge <- wedge
      .self
    },
    iterate = function(gradient, stepNumber) {
      .self$gradCache <- .self$decayMoments[[1]] * .self$gradCache +
        (1 - .self$decayMoments[[1]]) * gradient

      .self$gradSquare <- .self$decayMoments[[2]] * .self$gradSquare +
        (1 - .self$decayMoments[[2]]) * gradient %*% t(gradient)

      firstMoment <- .self$gradCache / (1 - .self$decayMoments[[1]] ^ stepNumber)
      secondMoment <- .self$gradSquare / (1 - .self$decayMoments[[2]] ^ stepNumber)

       -.self$stepSize * firstMoment / (sqrt(diag(secondMoment)) + .self$wedge)
    }
  )
)

dpOptimizerAdagrad <- setRefClass(
  Class = 'dpOptimizerAdagrad',
  contains = 'mechanismOptimize',
  fields = list(
    decay = 'numeric',
    gradSquare = 'numeric',
    wedge = 'numeric' # from (5), Adaptive Subgradient Methods for Online Learning and Stochastic Optimization, delta * I
  ),
  methods = list(
    initialize=function(..., gradSquare=NULL, decay=.9, wedge=1e-8) {
      callSuper(...)
      .self$gradSquare <- if(is.null(gradSquare)) matrix(0L, nrow=nrow(.self$theta), ncol=nrow(.self$theta)) else gradSquare

      .self$decay <- decay
      .self$wedge <- wedge
      .self
    },

    iterate = function(gradient, stepNumber) {
      # update stored hessian approximation
      .self$gradSquare <- .self$decay * .self$gradSquare + (1 - .self$decay) * gradient %*% t(gradient)

      # rescale each axis of the gradient
      -.self$stepSize * gradient / (sqrt(diag(.self$gradSquare)) + .self$wedge)
    }
  )
)


testLoss <- function() {
  expected <- c(1,2,1)

  # CASE 1: classifier-produced class probabilities
  predicted <- rbind(c(.6, .01, 0.1), c(.5, .8, .9))
  print(costFunctions$crossEntropy(expected, predicted))

  expected <- c(4, 4)
  predicted <- rbind(
    c(.25, .25, .25, .25),
    c(.01, .01, .01, .96)
  )
  print(costFunctions$crossEntropy(expected, predicted))

  # CASE 2: correct predictions
  # predicted <- matrix(c(1, 0, 1, 0, 1, 0), ncol=length(unique(expected)))
  # print(costFunctions$crossEntropy(expected, predicted))
}

# testLoss()

PUMSdata <- read.csv(file="/home/shoe/Desktop/MaPUMS5full.csv")
mydata<-PUMSdata[c("married","educ")]

theta <- matrix(0, 2, 1)

predict <- function(x, params)
  1/(1+exp(-as.matrix(cbind(1, x)) %*% params))

optimizer <- dpOptimizerSGD$new(
  epsilon=1, delta=1e-5,
  stepSize=c(1, .01),
  batchSize=nrow(mydata),
  predictors=c('educ'), variable='married',
  clippingInterval=10,
  forwardFun=predict, theta=theta, costFun='binaryCrossEntropy')
optimizer$release(mydata)

# results
print(optimizer$theta)

# baseline
print(coef(glm(married ~ educ, family="binomial", data=mydata)))
