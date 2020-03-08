context('unbiasedprivacy')

data <- data(PUMS5extract10000, package = 'PSIlence')

#TODO - add tests for validation

test_that('test current sim', {
   # param_combos <- read.csv('ellen_param_combos.csv')
   # param_row <- param_combos[5,]
    
#TODO - generate random seed?
 
    param_row <- list( N = 25000, 
                      sim=2,
                      alpha= 0.25,
                      e = 1,
                      P = 250,
                      seed = 61851,
                      rownum = 5,
                      intercept = 1,
                      beta = 3, 
                      y_var = 1000,
                      R = 100,
                      b = 100,
                      e_alpha = 1)
    
    # run UDP sim
    save_path = '.'
    returnValue = udpSim(param_row, save_path = save_path)
    
   # expect_equal(returnValue$theta_tilde, 3)
    #adjusted solution is in theta_tilde
   
})

#test_that('unbiased privacy sanity check', {
    # instantiate sample and aggregate for Mean with one subset and normal dpMean calculation
#    unbiased_privacy_test <- dpUnbiasedPrivacy(innerFun = 'mean', aggregationFun = 'Mean', numSubsets = 1,
#                                                      mechanism = 'mechanismLaplace', varType = 'numeric',
#                                                      variable = 'income', n = nrow(PUMS5extract10000), epsilon = 10^6, rng = c(0, 750000))
#    mean_test <- dpMean(mechanism = 'mechanismLaplace', varType = 'numeric',
#                                                      variable = 'income', n = nrow(PUMS5extract10000), epsilon = 10^6, rng = c(0, 750000))
    # release results
#    unbiased_privacy_test$release(data)
   
#})

# TODO: once sample-and-aggregate contains aggregationFun/mechanism combinations that are not compatible,
#       test that it throws an informative error