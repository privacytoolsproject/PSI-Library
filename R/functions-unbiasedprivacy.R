#' Bag of little bootstraps
#'
#' Implements the bag of little bootstraps algorithm (Kleiner et al., 2014). Splits the dataset and calculates the weighted QOI on each partition, also calculates alpha_1 and alpha_2
#
#' @param x Dataset partition
#' @param statistic Function that calculates quantity of interest
#' @param metric Summary metric (default is mean and variance)
#' @param P Number of bootstrap simulations
#' @param N Dataset n-size
#' @param lambda Bounding parameter for the QOI
#' @param ... Parameters required for \code{statistic}
#' @return 
#' Returns a list:
#' \item{metrics}{Summary metric across bootstraps: default is mean and variance}
#' \item{a_2}{Proportion of bootstrapped QOIs that are above lambda}
#' \item{a_1}{Proportion of bootstrapped QOIs that are above -lambda}
#' \item{res}{Vector of B bootstrapped QOIs}
#' 
#' @export
innerBLB <- function(x, statistic, metric, P, N, lambda, ...) {
  n <- nrow(x)
  resamples <- stats::rmultinom(P, N, rep(1 / n, n))
  res <- lapply(1:P, function(ii) {
    weights <- resamples[, ii]
    return(suppressWarnings(statistic(x, weights, ...)))
  })
  res <- unlist(res) #do.call(rbind, res)
  metrics <- unlist(lapply(1:length(metric), FUN = function(i) metric[[i]](res)))
  
  # alpha 
  a_2 <- mean(res > lambda)
  a_1 <- mean(res < -1*lambda)
  metrics <- c(metrics, a_1, a_2, res)
  return(metrics)
}

#' Split Data
#'
#' Splits the dataset into P partitions of size n
#' 
#' @param N Rows in dataset
#' @param n Number of rows desired in a single partition 
#' @param P Number of partitions
#' @param disjoint Whether samples should be allowed to overlap
#' 
#' @return 
#' \item{data_list}{List of row numbers}
#' 
#' @export
splitData <- function(N, n, P, disjoint = T) {
  #disjoint: T if we don't want samples to intersect
  if (disjoint == T) {
    data_list <- list()
    vec <- c(1:N)
    for (i in 1:P) {
      subsample <- sample(vec, n, replace = F)
      data_list[[i]] <- subsample
      vec <- vec[-which(vec %in% subsample)]
    }
  } else {
    data_list <- lapply(1:P, FUN = function(i) sample(1:N, n, replace = F))
  }
  
  return(data_list)
}


#' Censor the BLB parameters
#'
#' Censors \eqn{\hat \theta_p} at \eqn{\Lambda} and \eqn{-\Lambda} and adds DP noise.
#' 
#' @param param_vector Vector of quantities of interest calculated on each BLB partition
#' @param lambda Bounding parameter
#' @param delta privacy parameter
#' @param epsilon Privacy budget for QOI
#' @param epsilon_alpha Privacy budget for alpha
#' 
#' @return 
#' \item{param_hat_dp}{Censored QOI with added noise}
#' \item{param_noise}{SE of noise added to \code{param_hat_dp}}
#' \item{a_1}{Estimate of left-sided censoring \textbf{Not currently used}}
#' \item{a_2}{Estimate of right-sided censoring \textbf{Not currently used}}
#' \item{alpha_noise}{SE of noise to be added to \code{primary_alpha}}
#' 
#' @export
censorParam <- function(param_vector, lambda, delta, epsilon, epsilon_alpha){
  
  # calculate number of partitions that are censored + noise
  P  <- length(param_vector)
  S_1_a <- sqrt((log(delta))^2 - epsilon_alpha*log(delta)) + (epsilon_alpha/2 - log(delta))
  S_2_a <- 2*log(1.25/delta)
  D_a <- 2/P
  alpha_noise <- ((D_a/epsilon_alpha)* sqrt(min(S_1_a, S_2_a)))/2
  
  a_1 <- (1/P) * sum(param_vector < - lambda) + rnorm(1, 0, alpha_noise)
  a_2 <- (1/P) * sum(param_vector > lambda) + rnorm(1, 0, alpha_noise)
  
  # pulling alphas back to be between 0 and 1
  a_1 <- ifelse(a_1 > 1, 1, a_1)
  a_1 <- ifelse(a_1 < 0, 0, a_1)
  a_2 <- ifelse(a_2 > 1, 1, a_2)
  a_2 <- ifelse(a_2 < 0, 0, a_2)

  # do the censoring 
  param_censored <- param_vector
  param_censored[param_vector > lambda] <- lambda 
  param_censored[param_vector < -1*lambda] <- -1*lambda
  
  S_1 <- sqrt((log(delta))^2 - epsilon*log(delta)) + (epsilon/2 - log(delta))
  S_2 <- 2*log(1.25/delta)
  D <- 2*lambda/P
  param_noise <- (D/epsilon)* sqrt(min(S_1, S_2))
  param_hat_dp <- mean(param_censored) + rnorm(1, 0, param_noise)
  
  # return the mean of the censored theta and estimate of c
  #NOTE: a_1 and a_2 are returned here but are only used when we do not bias adjust. 
  #      when we bias adjust, we calculate a_1 and a_2 differently
  return(list(param_hat_dp, param_noise, a_1, a_2, alpha_noise))
}

#' Inverse \code{erf} function
#'
#' @import VGAM
#' @export
erfinv <- function(x){
  VGAM::erf(x, inverse = TRUE)
}


#' Upper bound bias adjustment equation
#'
#' Takes parameters and outputs function for upper bound bias adjustment equation.
#' 
#' @export
upperBound <- function(theta, sigma, alpha2, lambda){
  (theta + sigma*sqrt(2)*erfinv(2*(1-alpha2) - 1)) - lambda
}

#' Lower bound bias adjustment equation
#'
#' Takes parameters and outputs function for lower bound bias adjustment equation.
#' 
#' @export
lowerBound <- function(theta, sigma, alpha1, lambda){
  
  theta + sigma*sqrt(2)*erfinv(2*(alpha1) - 1) + lambda
}

#' Theta bias adjustment equation
#'
#' Takes parameters and outputs function for the third bias adjustment equation. 
#' 
#' @export
thetaHat <- function(theta, sigma, alpha1, alpha2, lambda, theta_dp){
  (1 - alpha2 - alpha1)*(theta + (sigma/sqrt(2*pi))*((exp(-0.5*((-lambda - theta)/sigma)^2) - exp(-0.5*((lambda - theta)/sigma)^2))/(1-alpha2 - alpha1))) + 
    lambda*(alpha2 - alpha1) - theta_dp
}

#' System of equations for right-sided censoring bias adjustment
#'
#' @export
boundFunUpper <- function(x, alpha2, theta_dp, lambda){
  # X is params you're solving for
  eqn1 <- upperBound(x[1], x[2], alpha2, lambda)
  eqn2 <- lowerBound(x[1], x[2], x[3], lambda)
  eqn3 <- thetaHat(x[1], x[2], x[3], alpha2, lambda, theta_dp)
  return(c(eqn1, eqn2, eqn3))
}

#' System of equations for left-sided censoring bias adjustment
#'
#' @export
boundFunLower <- function(x, alpha1, theta_dp, lambda){
  # X is params you're solving for
  eqn1 <- upperBound(x[1], x[2], x[3], lambda)
  eqn2 <- lowerBound(x[1], x[2], alpha1, lambda)
  eqn3 <- thetaHat(x[1], x[2], alpha1, x[3], lambda, theta_dp)
  return(c(eqn1, eqn2, eqn3))
}

#' Sigma estimate 
#' 
#' Calculates the standard deviation for the one-sided bias adjustment procedure
#'
#' @param theta Our DP censored estimate (\eqn{\hat \theta^{dp}})
#' @param c DP proportion of partitions censored (\eqn{\alpha^{dp}})
#' @param lambda Bounding parameter
#' @param upper indicator to calculate based on left or right-sided censoring
#' 
#' @export
#'
sigmaEstimate <- function(theta, c, lambda, upper){
  t <- sqrt(2)*erfinv(2*(1-c) - 1)
  if(upper){
    numerator <- 2*sqrt(pi)*(theta - lambda)
  }else{
    numerator <- 2*sqrt(pi)*(theta + 1*lambda)
  }
  denominator <- 2*sqrt(pi)*(c)*t - sqrt(2)*exp(-(t^2)/2) - 2*sqrt(pi)*t
  return(numerator/denominator)
}

#' Bias adjustment
#'
#'Calculates \eqn{\tilde \theta^{dp}} from \eqn{\hat \theta^{dp}}, which is biased due to censoring.
#' 
#' @param theta_dp Mean of censored QOI returned from BLB with dp noise added
#' @param alpha Primary alpha. Differentially private estimate of censoring on one side of the distribution of theta.
#' @param lambda Bounding parameter
#' @param upper Indicator for whether the alpha represents upper (right-sided) or lower (left-sided) censoring
#' @param two_sided Indicator for whether censoring is taking place on both sides or just one of the distribution of theta
#' 
#' @return 
#' \item{theta_tilde}{Bias corrected QOI}
#' \item{sigma}{Estimated sigma}
#' \item{alpha}{Estimated secondary alpha}
#' \item{two_sided_ba_ind}{Indicator for whether one or two-sided bias adjustment was used (just for debugging, not differentially private)}
#' 
#' @import nleqslv
#' 
#' @export
biasAdjustment <- function(theta_dp, alpha, lambda, upper, two_sided){
  # moving impossible thetas_dp to lambda
  theta_dp <- ifelse(theta_dp > lambda, lambda, theta_dp)
  theta_dp <- ifelse(theta_dp < -1*lambda, -1*lambda, theta_dp)
  
  if(two_sided){
    if(upper){
      soln <- tryCatch({nleqslv::nleqslv(c(theta_dp, 0.001, alpha/3), boundFunUpper, method = c('Broyden','Newton'), 
                                alpha2 = alpha, lambda = lambda, theta_dp = theta_dp, control = list(maxit = 1000))},
                       error = function(e){ret <- list('message' = paste0('Error', e))})
    }else{
      soln <- tryCatch({nleqslv::nleqslv(c(theta_dp, 0.001, alpha/3), boundFunLower, method = c('Broyden','Newton'), 
                                alpha1 = alpha, lambda = lambda, theta_dp = theta_dp, control = list(maxit = 1000))},
                       error = function(e){ret <- list('message' = paste0('Error', e))})
    }
    
    # return two sided results
    if(soln$message %in% c("Function criterion near zero", "x-values within tolerance 'xtol'")){
      two_sided_ba_ind <- 1
      ret <- list('theta_tilde' = soln$x[1], 'sigma' = soln$x[2], 'alpha' = soln$x[3], 
                  'two_sided_ba_ind' = two_sided_ba_ind)
    }else{
      # there was an error in two sided adj, proceed to the one sided
      two_sided <- FALSE
    }
  }

  # One sided bias adjustment procedure
  if(!two_sided){
    two_sided_ba_ind <- 0
    # Do one-sided bias adjustment
    if(upper){
      t <- sqrt(2)*erfinv(2*(1-alpha) - 1)
      B <- (1-alpha) + (sqrt(2)*exp(-(t^2/2)))/(2*t*sqrt(pi))
      theta_tilde <- theta_dp*(1/B) + lambda*((B-1)/B)
      sigma <- sigmaEstimate(theta_dp, alpha, lambda, upper)
      ret <- list('theta_tilde' = round(theta_tilde, 10), 'sigma' = ifelse(is.na(sigma), 0, sigma), 'alpha' = 0,
                  'two_sided_ba_ind' = two_sided_ba_ind)
    }else{
      K <- sqrt(2)*erfinv(2*alpha - 1)
      numerator <-  2*K*(sqrt(pi)*theta_dp - sqrt(pi)*-1*lambda)
      denominator <-  2*sqrt(pi)*alpha*K + sqrt(2)*exp(-(K^2)/2) - 2*sqrt(pi)*K
      theta_tilde <- -1*lambda - numerator/denominator
      sigma <- sigmaEstimate(theta_dp, alpha, lambda, upper)
      ret <- list('theta_tilde' = round(theta_tilde, 10), 'sigma' = ifelse(is.na(sigma), 0, sigma), 'alpha' = 0,
                  'two_sided_ba_find' = two_sided_ba_ind)
    }
  }
  return(ret)
}

phi <- function(x){
  (1/(sqrt(2*pi)) * exp(-0.5*(x^2)))
}

Phi <- function(x){
  0.5*(1 + VGAM::erf(x/sqrt(2)))
}

#' Censored variance estimate
#'
#' Estimates the variance of \eqn{\hat \theta^{dp}}.
#' 
#' @param theta Mean of censored QOI returned from BLB with dp noise added
#' @param a_1 DP proportion of partitions censored on the left
#' @param a_2 DP proportion of partitions censored on the right
#' @param lambda Bounding parameter 
#' @param sigma Estimated variance of BLB partitions
#' @param P Number of partitions
#' @param censoring_side Indicator for which side or both censoring occurred
#' 
#' @return 
#' \item{var}{Variance estimate of theta_hat_dp}
#' 
#' @export
censoredVarianceEstimate <- function(theta, theta_hat_dp, a_1, a_2, lambda, theta_noise, sigma, P, censoring_side){
  #theta: our dp censsored estimate
  #a_1, a_2: our dp proportion of partitions censored
  #lambda: our bound
  #theta_noise: dp-noise added to theta estimate
  #sigma: output of sigmaEstimate
  
  alpha <- (-1*lambda - theta)/sigma
  beta <- (lambda - theta)/sigma 
  
  part1 <- (alpha*phi(alpha) - beta*phi(beta))/(Phi(beta) - Phi(alpha))
  part2 <- ((phi(alpha) - phi(beta))/(Phi(beta) - Phi(alpha)))^2
  new_P <- P*(1 - a_1 - a_2)
  truncated_var <- (1/new_P)*(sigma^2)*(1 + part1 - part2)

  part0 <- theta + sigma * ((phi(alpha) - phi(beta))/(Phi(beta) - Phi(alpha)))
  part3 <- (lambda^2)*(a_2 + a_1) - (theta_hat_dp)^2
  
  var <- (1/P) * ((1 - a_2 - a_1)*(part0^2 + new_P*truncated_var) + part3)
  
  return(var)
}


#' Variance simulation
#'
#' Estimates the variance of \eqn{\tilde \theta^{dp}} via simulation.
#' 
#' @param theta_tilde Estimate of QOI
#' @param sigma_hat Estimate of variance  
#' @param theta_hat_dp_overall Differentially private censored QOI
#' @param a_1 DP estimate of proportion of partitions that were left-censored
#' @param a_2 DP estimate of proportion of partitions that were right-censored
#' @param P Number of partitions
#' @param lambda Bounding parameter
#' @param theta_noise SD of dp noise added to theta_hat_dp
#' @param alpha_noise SD of dp noise added to alpha_dp
#' @param nsims Number of simulations to boostrap variance
#' 
#' @return 
#' \item{var_est}{Variance estimate of theta_tilde}
#' \item{theta_tilde_sims}{Vector of simulated theta_tilde estimates}
#' \item{mvn_draws}{Matrix of theta_hat_dp and alpha draws from the multivariate normal}
#' \item{mat_not_pos_def}{Indicator for whether an adjustment was applied to fix numerical issue where the covariance matrix can appear to be not positive definite}
#' \item{sigma_matrix}{Covariance matrix used in the multivariate normal draws}
#' \item{var_theta_hat_nonoise}{Estimated variance of theta_hat_dp not including variance from dp noise added}
#' \item{orig_sigma_matrix}{Covariance matrix prior to adjustment to be positive definite}
#' 
#' @import MASS 
#' 
#' @export
varianceSimulation <- function(theta_tilde, sigma_hat, theta_hat_dp_overall, a_1, a_2, P, lambda, theta_noise, alpha_noise, nsims){
  # which side are we censoring on?
  censoring_side <- ifelse(theta_hat_dp_overall > 0, 'right', 'left')
  censoring_side <- ifelse(a_1 > 0.1 & a_2 > 0.1, 'both', censoring_side)

  if(theta_hat_dp_overall > 0){
    cov_hat <- (1/P) * (lambda - ((theta_hat_dp_overall - a_2 * lambda + a_1 * lambda)/(1 - a_2 - a_1))) * a_2*(1 - a_2)
  }else{
    cov_hat <- (1/P) * (-1*lambda - ((theta_hat_dp_overall - a_2 * lambda + a_1 * lambda)/(1 - a_2 - a_1))) * a_1*(1 - a_1)
  }

  var_theta_hat_nonoise <- censoredVarianceEstimate(theta_tilde, theta_hat_dp_overall, a_1, a_2, lambda, theta_noise, sigma_hat, P, censoring_side)

  var_theta_hat <- var_theta_hat_nonoise + (theta_noise^2)
  var_alpha <- ifelse(censoring_side %in% c('both','right'),  (1/P)*(1 - a_2)*a_2 + (alpha_noise^2),  (1/P)*(1 - a_1)*a_1 + (alpha_noise^2))

  # draw from the MVN and bias adjust
  sigma_mat <- matrix(c(var_theta_hat, cov_hat, cov_hat, var_alpha), ncol = 2)
  
  orig_sigma_mat <- sigma_mat
  # fixing non positive definite matrices
  # NOTE: this is a floating point issue (I think)
  if(det(sigma_mat) <= 0){
    sigma_mat <- lqmm::make.positive.definite(sigma_mat)
    mat_not_pos_def <- 1
  }else{
    mat_not_pos_def <- 0
  }
  
  draws <- MASS::mvrnorm(n = nsims, mu = c(theta_hat_dp_overall, a_2), Sigma = sigma_mat)
  draws[,2][draws[,2] < 0.1] <- 0.1 # fixing alpha draws that are too high or too low 
  draws[,2][draws[,2] > 0.9] <- 0.9
  fix_inds <- rep(0, nrow(draws))
  draws <- cbind(draws, fix_inds)
  while(sum(abs(draws[,1] >= lambda)) > 0){
    draws[,1][draws[,1] >= lambda] <- draws[,1][draws[,1] >= lambda] - (theta_noise * sqrt(2/pi))
    draws[,1][draws[,1] <= -1*lambda] <- draws[,1][draws[,1] <= -1*lambda] + (theta_noise * sqrt(2/pi))
    
    draws[,3][draws[,1] >= lambda] <-  draws[,3][draws[,1] >= lambda] + 1
    draws[,3][draws[,1] <= -1*lambda] <-  draws[,3][draws[,1] <= -1*lambda] + 1
  }
  theta_tilde_sims <- list()
  two_sided <- ifelse(a_2 > 0.1 & a_1 > 0.1, TRUE, FALSE)
  for(i in 1:nsims){
    bias_adj_sim <- biasAdjustment(draws[i, 1], draws[i, 2], lambda, upper = a_2 > a_1, two_sided = two_sided)
    theta_tilde_sims[[i]] <- bias_adj_sim$theta_tilde
  }
  
  theta_tilde_sims <- unlist(theta_tilde_sims)
  
  
  ret <- list('var_est' = var(theta_tilde_sims),
              'theta_tilde_sims' = theta_tilde_sims,
              'mvn_draws' = draws,
              'mat_not_pos_def' = mat_not_pos_def,
              'sigma_matrix' = sigma_mat,
              'var_theta_hat_nonoise' = var_theta_hat_nonoise,
              'orig_sigma_matrix' = orig_sigma_mat)
  return(ret)
}


#' Differential privacy algorithm
#'
#' Calculates a quantity of interest in a differentially-private way. Note that many returned items are not differentially-private and are simply used for debugging and illustrative purposes.
#
#' @param data Input data
#' @param statistic Function that calculates quantity of interest
#' @param B Number of bootstraps to run via BLB algorithm
#' @param n Split size
#' @param P Number of partitions 
#' @param lambda Bounding parameter for the QOI
#' @param lambda_var Bounding parameter for the variance
#' @param delta Privacy parameter
#' @param epsilon Privacy budget for the QOI
#' @param epsilon_alpha Privacy budget for estimating alpha^{dp}
#' @param censoring_cutoff Maximum amount of censoring to allow
#' @param bias_cutoff Maximum amount of censoring to allow withou doing bias correction
#' @param parallelize Whether to parallelize the BLB calculations
#' @param ... Parameters necessary for \code{statistic}
#' 
#'
#' @return 
#' \item{theta_tilde}{Differentially private estimate of quantity of interest}
#' \item{theta_hat}{Differentially private estimate of QOI, unadjusted for bias introduced by censoring}
#' \item{var_est}{Estimate of variance of theta_tilde}
#' \item{a_1}{Estimate of left-sided censoring}
#' \item{a_2}{Estimate of right-sided censoring}
#' \item{alpha_noise}{SD of differentially private noise added to alpha estimate}
#' \item{theta_noise}{SD of differentially private noise added to theta_hat estimate}
#' \item{blb_thetas}{Vector of QOIs calculated in each partition during bag of little bootstraps procedure}
#' \item{sigma_hat}{Estimated SD of true QOI}
#' \item{alpha_too_high_halt}{Indicator for whether alpha was greater than the censoring cutoff}
#' \item{bias_adj_no_converge}{Indicator when bias adjustment procedure has failed}
#' \item{theta_tilde_var_sims}{Simulated theta_tilde draws from the variance simulation}
#' \item{mvn_draws}{Matrix of draws from the multivariate normal from the variance simulation}
#' \item{var_sigma_mat_not_pos_def}{Indicator for whether the variance simulation covariance matrix needed to be adjusted}
#' \item{sigma_marix}{Covariance matrix used in variance simulation}
#' \item{var_theta_hat_dp_nonoise}{Variance of theta_hat_dp before accounting for variance introduced by dp noise}
#' \item{orig_sigma_mat}{Covariance matrix used in variance simulation before any adjustment to make it positive definite}
#' \item{fix_indicator}{Indicator for whether theta_hat_dp was outside the range of the bounding parameter and was brough back in}
#' \item{two_sided_ba_ind}{Indicator for whether the two sided bias adjustment procedure was used}
#'
#' @examples
#' \dontrun{algorithmUDP(dat, statistic = coefFn, B = 100, n = 100, P = 1000, lambda = 3.1,
#' lambda_var = 0.025, form = as.formula(Y1 ~ X), coef = 'X')}
#' 
#' @import parallel
#' 
#' @export

algorithmUDP <- function(data, statistic, B, n, P, lambda, lambda_var, delta, epsilon = 0.1, epsilon_alpha = 0.1, 
                         censoring_cutoff = 0.6, bias_cutoff = 0.1, parallelize = F, ...) {
  
  # STEP 0: Assign placeholder variables and check inputs
  var_sim <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA,NA)
  if(delta > 1){
    stop('ERROR: Delta must be between 0 and 1!')
  }
  
  # STEP 1: Split the sample into subsets
  N <- nrow(data)
  data_splits <- splitData(N, n, P, disjoint = T)
  
  # STEP 2: for each split, bootstrap statistic and return the QOI calculated on bootstrap distribution
  # define metric 
  metric <- list('mean' = mean, 'var' = var)
  
  if(parallelize){
    var_list <- append(list('innerBLB','data','statistic','B','n', 'metric', 'lambda'), names(list(...)))
    
    # we have to export the variables in our environment in order to use them with the parallel package
    # for ... arguments that are user specified, we need to assign them to their names in this environ
    for(i in 1:length(list(...))){
      assign(names(list(...)[i]), list(...)[[i]])
    }
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    parallel::clusterExport(cl, var_list, envir=environment())
    est <- parallel::parLapply(cl, data_splits, fun = function(i)
      innerBLB(data[i, ], statistic, metric, B, N, lambda, form, coef))
    parallel::stopCluster(cl)
  }
  else{
    est <- lapply(data_splits, FUN = function(i) innerBLB(data[i, ], statistic, metric, 
                                                          B, N, lambda, ...))
  }
  
  #store estimates of mean and variance of QOI and alpha_1 and alpha_2
  est <- do.call(rbind, est)
  s_est <- est[,1]
  v_est <- est[,2]
  a1_est <- est[,3]
  a2_est <- est[,4]
  
  
  # STEP 3: censor and calculate theta_hat_dp and alpha_hat_dp
  cen <- censorParam(s_est, lambda, delta, epsilon, epsilon_alpha)
  theta_hat_dp <- cen[[1]]
  theta_noise <- cen[[2]] # S^2
  alpha_noise <- cen[[5]] # S_{alpha}^2
  
  # Determine which alpha is our primary alpha (the side with more censoring) and calculate alpha^{dp}
  if(theta_hat_dp > 0){
    a_2 <- mean(a2_est) + rnorm(1, mean = 0, sd = alpha_noise)
    a_2  <- ifelse(a_2 < 0, 0.1, a_2)
    upper <- TRUE 
    primary_alpha <- a_2
  }else{
    a_1 <- mean(a1_est) + rnorm(1, mean = 0, sd = alpha_noise)
    a_1 <- ifelse(a_1 < 0, 0.1, a_2)
    upper <- FALSE 
    primary_alpha <- a_1
  }

  # STEP 3.1: Fail if censoring is above censoring cutoff
  if(primary_alpha > censoring_cutoff){
    stop('ERROR: Censoring was above set censoring cutoff.')
  }
  
 
  # STEP 4: Adjust for bias if necessary
  # Assign indicator whether we should use the one sided or two sided bias adjustment equations
  if(mean(a2_est) >= 0.1 & mean(a1_est) >= 0.1){
    two_sided <- T
  }else{
    two_sided <- F
  }
  
  # Employ "the fix" if necessary: 
  # reassign values of theta_hat_dp that are impossible (greater than lambda, less than -lambda)
  theta_hat_dp_before_fix <- theta_hat_dp
  fix_indicator <- 0
  while(abs(theta_hat_dp) >= lambda & fix_indicator < 20){
    fix_indicator <- fix_indicator + 1
    if(theta_hat_dp <= -1*lambda){
      theta_hat_dp <- theta_hat_dp + (theta_noise * sqrt(2/pi))
    }else{
      theta_hat_dp <- theta_hat_dp - (theta_noise * sqrt(2/pi))
    }
  }
  
  # pulling back thetas that the fix didn't fix
  # NOTE: the following lines shouldn't be necessary and are legacy code
  extra_fix_indicator <- 0
  if(abs(theta_hat_dp) > lambda){
    extra_fix_indicator  <- 1
  }
  theta_hat_dp <- ifelse(theta_hat_dp > lambda, lambda, theta_hat_dp)
  theta_hat_dp <- ifelse(theta_hat_dp < -1*lambda, -1*lambda, theta_hat_dp)
  
  fix_noise <- theta_noise * sqrt(2/pi)
  
  
  
  # Adjust for bias
  bias_adj <- biasAdjustment(theta_hat_dp, primary_alpha, lambda, upper = upper, two_sided = two_sided)

  # If the bias adjustment fails, fail
  if(bias_adj == 'Error'){
    stop('ERROR: bias adjustment procedure failed to converge.')
  }else{
    # If bias adjustment succeeeds, assign results!
    theta_tilde <- round(unlist(bias_adj[1]), 10)
    sigma_est <- unlist(bias_adj[2])
    
    # assign secondary alpha, estimated via bias adjustment
    if(theta_hat_dp > 0){
      a_1 <- unlist(bias_adj[3])
      secondary_alpha <- a_1
    }else{
      a_2 <- unlist(bias_adj[3])
      secondary_alpha <- a_2
    }
  }
  
  # STEP 5: Calculate variance of theta_tilde
  # Simulate the variance if we have used the bias adjustment procedure
  if(a_1 + a_2 >= bias_cutoff){
    bias_adj_ind <- 1
    #var_sim <- varianceSimulation(theta_tilde, sigma_est, theta_hat_dp, a_1, a_2, P, lambda, theta_noise, alpha_noise, nsims = 1000)
    var_sim <- varianceSimulation(theta_tilde, sigma_est, theta_hat_dp_before_fix, a_1, a_2, P, lambda, theta_noise, alpha_noise, nsims = 1000)
    var_est <- var_sim[[1]]

  }else{
    # Otherwise, calculate the variance from what we have estimated via BLB
    bias_adj_ind <- 0
    
    # reining in impossible values for theta_hat_dp
    theta_hat_dp <- ifelse(theta_hat_dp > lambda, lambda, theta_hat_dp)
    theta_hat_dp <- ifelse(theta_hat_dp < -1*lambda, -1*lambda, theta_hat_dp)
    
    # Censor and estimate variance
    var_censored <- censorParam(v_est, lambda = lambda_var, delta, epsilon, epsilon_alpha)
    
    var_hat_dp <- var_censored[[1]]
    var_noise <- var_censored[[2]]
    a_1_var <- var_censored[[3]]
    a_2_var <- var_censored[[4]]
    alpha_noise_var <- var_censored[[5]]
    
    # Adjust variance for bias if it was censored
    if(a_2_var > bias_cutoff){
      var_est <- biasAdjustment(var_hat_dp, a_2_var, lambda_var, two_sided = two_sided) + (theta_noise)^2
    }else{
      var_hat_dp <- ifelse(var_hat_dp > lambda_var, lambda_var, var_hat_dp)
      var_hat_dp <- ifelse(var_hat_dp < -1*lambda_var, -1*lambda_var, var_hat_dp)
      
      var_est <- var_hat_dp  + (theta_noise)^2
    }
    
    theta_tilde <- theta_hat_dp # we have not censored so theta_hat_dp = theta_tilde
  }
  


 
  # STEP 6: return relevant quantities
  ret <- list('theta_tilde' = ifelse(exists("theta_tilde"), theta_tilde, NA), 
              'theta_hat' = ifelse(exists("theta_hat_dp"), theta_hat_dp, NA),
              'theta_blb' = ifelse(exists("s_est"), mean(s_est), NA),
              'var_est' = ifelse(exists("var_est"), var_est, NA),
              'a_1' = ifelse(exists("a_1"), a_1, NA),
              'a_2' = ifelse(exists("a_2"), a_2, NA),
              'alpha_noise' = ifelse(exists("alpha_noise"), alpha_noise, NA),
              'theta_noise' = ifelse(exists("theta_noise"), theta_noise, NA),
              'blb_thetas' = ifelse(exists("s_est"), s_est, NA),
              'sigma_hat' = ifelse(exists("sigma_est"), sigma_est, NA),
              
              # params we specified to run the algorithm
              'P' = P,
              'lambda' = lambda,
              'B' = B,
              'n' = n,
              'lambda_var' = lambda_var,
              'delta' = delta,
              'epsilon' = epsilon,
              'epsilon_alpha' = epsilon_alpha,
              
              # error messages and branching info (mainly useful for debugging)
              'alpha_too_high_halt' = 0,
              'bias_adj_no_converge' = 0,
              'bias_adj_ind' = ifelse(exists("bias_adj_ind"), bias_adj_ind, NA),
              'theta_tilde_var_sims' = var_sim[[2]],
              'mvn_draws' = var_sim[[3]],
              'var_sigma_mat_not_pos_def' = var_sim[[4]],
              'sigma_matrix' = var_sim[[5]],
              'var_theta_hat_dp_nonoise' =  var_sim[[6]],
              'orig_sigma_mat' = var_sim[[7]],
              'fix_indicator' = fix_indicator,
              'two_sided_ba_ind' = unlist(bias_adj[4]),
              'fix_noise' = fix_noise,
              'extra_fix_indicator' = extra_fix_indicator,
              'theta_hat_dp_before_fix' = theta_hat_dp_before_fix
  )
  return(ret)
}

