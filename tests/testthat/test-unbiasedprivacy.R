context('unbiasedprivacy')

#data <- data(PUMS5extract10000, package = 'PSIlence')

test_that('test mean', {
  if (TRUE) {
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
    pr <- param_row
    dat <- generateData(pr$N, pr$intercept, pr$beta, pr$y_var, seed = pr$seed) # generate new dataset
    
    # Calculate lambda based on specified value of alpha (prop. censoring) and OLS SE
    l <- trueLambdaCalc(pr$beta, dat$X, pr$P, pr$y_var, pr$alpha)
    #true_sigma <- se * sqrt(pr$P)
    
    # QOI
    form <- as.formula(Y1~X)
    coefVal <- 'Y3'
    
    
    unbiased_privacy_test <- dpUnbiasedPrivacy(statistic = meanFn,
                                               B=param_row$R,
                                               n=param_row$b,
                                               P=param_row$P,
                                               lambda=l,
                                               lambda_var = 0.025,
                                               delta = 0.01,
                                               epsilon = param_row$e, 
                                               epsilon_alpha = param_row$e_alpha,
                                               censoring_cutoff = 0.9,
                                               bias_cutoff = 0.1 )
    
    returnValue <- unbiased_privacy_test$release(data = dat, form = form, coefVal = coefVal)       
    expect_equal(returnValue$theta_tilde, 0, tolerance=.5)
  }                                      
  
  
  

})

test_that('unbiased privacy sanity check', {
    
    #TODO: remove this TRUE/FALSE condition
    
    if (TRUE) {
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
    pr <- param_row
    dat <- generateData(pr$N, pr$intercept, pr$beta, pr$y_var, seed = pr$seed) # generate new dataset

    # Calculate lambda based on specified value of alpha (prop. censoring) and OLS SE
    l <- trueLambdaCalc(pr$beta, dat$X, pr$P, pr$y_var, pr$alpha)
    #true_sigma <- se * sqrt(pr$P)

    # QOI
    form <- as.formula(Y1~X)
    coefVal <- 'X'
  
   
    unbiased_privacy_test <- dpUnbiasedPrivacy(statistic = coefFn,
                                               B=param_row$R,
                                               n=param_row$b,
                                               P=param_row$P,
                                               lambda=l,
                                               lambda_var = 0.025,
                                               delta = 0.01,
                                               epsilon = param_row$e, 
                                               epsilon_alpha = param_row$e_alpha,
                                               censoring_cutoff = 0.9,
                                               bias_cutoff = 0.1 )
                                               
       returnValue <- unbiased_privacy_test$release(data = dat, form = form, coefVal = coefVal)   
       browser()
       expect_equal(returnValue$theta_tilde, 3, tolerance=.1)
    }                                      


   
})
#TODO - add tests for validation

test_that('test current sim', {
    # param_combos <- read.csv('ellen_param_combos.csv')
    # param_row <- param_combos[5,]
    
    #TODO - generate random seed?
    if (FALSE) {
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
    
    expect_equal(returnValue$theta_tilde, 3, tolerance=1e-2)
    }
    #adjusted solution is in theta_tilde
    
})
