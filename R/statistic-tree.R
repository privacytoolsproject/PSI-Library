# ToDo: add accuracy stuff, data imputation, percentiles, alpha
# ToDo: n.bins needs to be changed to nBins.
dpTreeStatistic <- setRefClass(
  Class = 'dpTreeStatistic',
  contains = 'mechanismLaplace',
  methods = list(
    initialize = function(variable, n.bins, rng, epsilon=NULL, impute.rng=NULL){
      .self$name <- 'differentially private binary tree'
      .self$mechanism <- 'mechanismLaplace'
      .self$variable <- variable
      .self$n.bins <- n.bins
      .self$rng <- rng
      .self$accuracy <- accuracy
    },
    release = function(data){
      x <- data[,variable]
    }
  )
)