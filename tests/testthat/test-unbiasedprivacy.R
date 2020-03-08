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
    expect_equal(returnValue$theta_tilde, 2.998053, tolerance=1e-5)
    #adjusted solution is in theta_tilde
   
})

#test_that('unbiased privacy sanity check', {
    # instantiate sample and aggregate for Mean with one subset and normal dpMean calculation
#    pr <- param_row
#    dat <- generateData(pr$N, pr$intercept, pr$beta, pr$y_var, seed = pr$seed) # generate new dataset

    # Calculate lambda based on specified value of alpha (prop. censoring) and OLS SE
#    l <- trueLambdaCalc(pr$beta, dat$X, pr$P, pr$y_var, pr$alpha)
    #true_sigma <- se * sqrt(pr$P)

    # QOI
#    form <- as.formula(Y1~X)
#    coef <- 'X'

#    sim <- algorithmUDP(data = dat, statistic = coefFn, B = pr$R, n = pr$b, P = pr$P, lambda = l, lambda_var = 0.025, delta = 0.01,
#                    epsilon = pr$e, epsilon_alpha = pr$e_alpha, parallelize = F, censoring_cutoff = 0.9,
#                    bias_cutoff = 0.1, form = form, coef = coef)

   
#})

# TODO: once sample-and-aggregate contains aggregationFun/mechanism combinations that are not compatible,
#       test that it throws an informative error