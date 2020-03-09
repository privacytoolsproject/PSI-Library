#  This is being done in the test file now
#TODO remove if not needed
#source('functions-unbiasedprivacy.R')
#source('utilities-unbiasedprivacy.R')

#Commenting these because we are trying to do it in testthat
#param_combos <- read.csv('ellen_param_combos.csv')

#param_row <- param_combos[5,]

#param_row

# run UDP sim
#save_path = '.'
#udpSim(param_row, save_path = save_path)

#FROM UNIT TEST - 
# instantiate sample and aggregate for Mean with one subset and normal dpMean calculation
param_row <- list( N = 25000, 
                   sim=2,
                   alpha= 0.25,
                   e = 1,
                   P = 250,
                   seed = 61852,
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
coef <- 'X'
# initialize = function(statistic, B, n, P, lambda, lambda_var, delta, epsilon = 0.1, epsilon_alpha = 0.1, 
#                          censoring_cutoff = 0.6, bias_cutoff = 0.1, ...) {

#    sim <- algorithmUDP(data = dat, 
#                         statistic = coefFn, 
#                         B = pr$R, 
#                         n = pr$b, 
#                         P = pr$P, 
#                         lambda = l, 
#                         lambda_var = 0.025, 
#                         delta = 0.01,
#                         epsilon = pr$e, 
#                         epsilon_alpha = pr$e_alpha, 
#                         parallelize = F, 
#                         censoring_cutoff = 0.9,
#                         bias_cutoff = 0.1, 
#                         form = form, 
#                         coef = coef)
#TODO: remove this TRUE/FALSE condition
if (TRUE) {
  
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
 
 unbiased_sim <-  unbiased_privacy_test$release(data = dat, form = form, coef = coef)  
  save_path <- '.'
  fname <- paste0(save_path, '/sim_', pr$seed, '.Rdata')
  save(unbiased_sim, file = fname)
}                                      

