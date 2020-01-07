context('sampleandaggregate')

data <- data(PUMS5extract10000, package = 'PSIlence')

test_that('errors are thrown for illegitimate function/mechanism calls', {
    # illegitimate innerFun
    expect_error(dpSampleAndAggregate(innerFun = 'max', aggregationFun = 'Mean', numSubsets = 10,
                             mechanism = 'mechanismLaplace', varType = 'numeric',
                             variable = 'income', n = nrow(PUMS5extract10000), epsilon = 0.1, rng = c(0, 750000)))

    # illegitimate aggregationFun
    expect_error(dpSampleAndAggregate(innerFun = 'mean', aggregationFun = 'Median', numSubsets = 10,
                             mechanism = 'mechanismLaplace', varType = 'numeric',
                             variable = 'income', n = nrow(PUMS5extract10000), epsilon = 0.1, rng = c(0, 750000)))
    # illegitimate mechanism
    expect_error(dpSampleAndAggregate(innerFun = 'mean', aggregationFun = 'Mean', numSubsets = 10,
                             mechanism = 'mechanismExponential', varType = 'numeric',
                             variable = 'income', n = nrow(PUMS5extract10000), epsilon = 0.1, rng = c(0, 750000)))

})

test_that('errors are thrown when numSubsets is invalid', {
    # numSubsets > n
    expect_error(dpSampleAndAggregate(innerFun = 'mean', aggregationFun = 'Mean', numSubsets = nrow(PUMS5extract10000)+1,
                             mechanism = 'mechanismLaplace', varType = 'numeric',
                             variable = 'income', n = nrow(PUMS5extract10000), epsilon = 0.1, rng = c(0, 750000)))

    # numSubsets = 0
    expect_error(dpSampleAndAggregate(innerFun = 'mean', aggregationFun = 'Mean', numSubsets = 0,
                            mechanism = 'mechanismLaplace', varType = 'numeric',
                            variable = 'income', n = nrow(PUMS5extract10000), epsilon = 0.1, rng = c(0, 750000)))
})

test_that('sample-and-aggregate is approximately equal to normal DP statistic when numSubsets = 1', {
    # instantiate sample and aggregate for Mean with one subset and normal dpMean calculation
    # NOTE: we do so with very large epsilon in order to
    sample_and_aggregate_test <- dpSampleAndAggregate(innerFun = 'mean', aggregationFun = 'Mean', numSubsets = 1,
                                                      mechanism = 'mechanismLaplace', varType = 'numeric',
                                                      variable = 'income', n = nrow(PUMS5extract10000), epsilon = 10^6, rng = c(0, 750000))
    mean_test <- dpMean(mechanism = 'mechanismLaplace', varType = 'numeric',
                                                      variable = 'income', n = nrow(PUMS5extract10000), epsilon = 10^6, rng = c(0, 750000))

    # release results
    sample_and_aggregate_test$release(data)
    mean_test$release(data)

    # ensure that the releases are approximately equal
    expect_equal(round(sample_and_aggregate_test$result$release, 0.1), round(mean_test$result$release, 0.1))
})

# TODO: once sample-and-aggregate contains aggregationFun/mechanism combinations that are not compatible,
#       test that it throws an informative error