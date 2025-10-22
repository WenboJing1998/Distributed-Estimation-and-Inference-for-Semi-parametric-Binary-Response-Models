library(foreach, warn.conflicts = F, quietly = T)
library(doParallel, warn.conflicts = F, quietly = T)
library(doRNG, warn.conflicts = F, quietly = T)
library(psych, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(optimization, warn.conflicts = F, quietly = T)
library(MASS, warn.conflicts = F, quietly = T)

#os <- "C://Users//92865//Dropbox//SMSE//p=10//"
#source(paste(os, "MSE_functions.R", sep=""))

source("MSE_functions.R")

noiseType = "norm"
m <- 1000 
p <- 1
logmn <- 1.6
n <- floor(m ^ {logmn})
L <- floor(n / m)
errorTolerence <- 0.001
maxIterations <- 10
confidenceLevel <- 0.05
z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
lambda <- 1 
alpha <- 2 # kernel order
beta_true <- rep(1/sqrt(p), p)
local_sizes <- rep(m, L) # local size on each machine

data = data_generate(n, p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
x <- data$x[1:sum(local_sizes), ]
y <- data$y[1:sum(local_sizes)]
minIterations = 4
v0 <- rep(1, p)
kernelOrder = 2
lambda_init = 1
delta=0.1
beta_init = "SMSE"
lambda_opt = F

 ## Tune c_init

results_cr <- rep(NA, 24)
results_mean <- rep(NA, 24)
results_median <- rep(NA, 24)
results_var <- rep(NA, 24)

for(logmn in seq(1.3, 1.9, .05)){
  print(logmn)
  m <- 1000
  n <- floor(m ^ {logmn})
  L <- floor(n / m)
  local_sizes <- rep(m, L)
  c_h_init_seq = c(1, 2, 3)
  cl <- makeCluster(6)
  registerDoParallel(cl, cores=6)
  results <- foreach(rep=1:100, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                     .errorhandling = 'remove') %dorng% {
                # for(rep in 1:100){
                       
                # print(rep)
                       
                # beta_temp <- rnorm(p, 0, 0.5)
                # beta_true <- beta_temp / sqrt(sum((beta_temp)^2))
                # sum(beta_true)  
                       
                  beta_true <- rep(1/sqrt(p), p)
                       
                  data <- data_generate(n, p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
                  x <- data$x[1:sum(local_sizes), ]
                  y <- data$y[1:sum(local_sizes)]
                       
                  beta_0 <- rep(0, p)

                  res_cr_each_rep <- c()
                  res_bias_each_rep <- c()
                       
                  for(c_h in seq(1, 1.5, 0.25)){
                    #for(c_h_init in seq(1, 2, 0.5)){
                    
                        res_DCMSE <- Ave_MSE(x, y, local_sizes = local_sizes)
                        beta_1_bias <- sum(res_DCMSE$beta_hat) - sum(beta_true)
                        coverage_1  <- abs(sum(res_DCMSE$beta_hat) - sum(beta_true)) < z * res_DCMSE$se
                        time_1 <- res_DCMSE$time
                           
                        res_DCSMSE <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                                              h_option="n", local_sizes=local_sizes)
                        beta_2_bias <- sum(res_DCSMSE$beta_hat) - sum(beta_true)
                        # time_2 <- res_DCSMSE$time
                        coverage_2  <- abs(sum(res_DCSMSE$beta_hat) - res_DCSMSE$bias - sum(beta_true)) < z * res_DCSMSE$se
                           
                        beta_pooled  <- Smoothed_MSE(x, y, beta_init=beta_0, inference = T, c_h=c_h)
                        beta_3_bias <- sum(beta_pooled$beta_hat) - sum(beta_true)
                        # time_3 <- beta_pooled$time
                        coverage_3 <- abs(sum(beta_pooled$beta_hat) - sum(beta_true) - beta_pooled$bias) < z * beta_pooled$se
                          
                        
                        results_mSMSE <- mSMSE(x, y, local_sizes, v0=rep(1, p), minIterations=5, maxIterations=7, 
                                               lambda_opt = F, c_h=c_h, c_h_init=c_h_init_seq,
                                               verbose = F)
                        # sum(results_mSMSE$beta_0 - beta_true)
                        beta_4_bias <- results_mSMSE$beta_each_iteration[1:5] - sum(beta_true)
                        coverage_4 <- apply(as.matrix(results_mSMSE$interval_each_iteration[ ,1:5]), 2, function(x) (sum(beta_true) < x[2]) & (sum(beta_true) > x[1]))
                           
                        cr <- c(coverage_1, coverage_2, coverage_3, coverage_4)
                        methods_names <- c("AvgMSE", "AvgSMSE", "pooled", paste("round", 1:5, sep=""))
                        names(cr) <- c(paste(methods_names, "_cr_", c_h, sep=""))
                        res_cr_each_rep <- c(res_cr_each_rep, cr)
                           
                        bias <- c(beta_1_bias, beta_2_bias, beta_3_bias, beta_4_bias)
                        names(bias) <- c(paste(methods_names, "_bias_", c_h, sep=""))
                        res_bias_each_rep <- c(res_bias_each_rep, bias)
                     # }
                    }
                       
                      return(c(res_cr_each_rep, res_bias_each_rep))
                       #res_cr_tune <- rbind(res_cr_tune, res_cr_each_rep)
                       #res_bias_tune <- rbind(res_bias_tune, res_bias_each_rep)
                }
  stopCluster(cl)
  results_cr <- rbind(results_cr, apply(results[ , 1:24], 2, mean))
  results_mean <- rbind(results_mean, apply(results[ , 25:48], 2, mean))
  results_median <- rbind(results_median, apply(results[ , 25:48], 2, median))
  results_var <- rbind(results_var, apply(results[ , 25:48], 2, var))
  print(rbind(results_cr, results_mean))
}

results_cr[-1, ]
results_mean[-1, ]
results_median[-1, ]
results_var[-1, ]


results[1:80, ]
results[1:199, ][results[1:199, 6] == 0, ]
results_record_1 <- results_all[-1, ]



results <- foreach(rep=1:102, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                            .errorhandling = 'remove') %dorng% {
# for(rep in 1:100){
  
  # print(rep)
                              
  #beta_temp <- rnorm(p, 0, 1)
  #beta_true <- beta_temp / sqrt(sum((beta_temp)^2))
  #sum(beta_true)  
                              
  beta_true <- rep(-1/sqrt(p), p)
  data <- data_generate(n, p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
  x <- data$x[1:sum(local_sizes), ]
  y <- data$y[1:sum(local_sizes)]
  
  beta_0 <- rep(0, p)
  # beta_true <- rep(-1/sqrt(p), p)

  
  res_cr_each_rep <- c()
  res_bias_each_rep <- c()
  
  for(c_h in seq(1.5, 1.5, 0.5)){
    for(c_h_init in seq(1.5, 2, 0.5)){
      
      res_DCSMSE <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                            h_option="n", local_sizes=local_sizes)
      
      beta_2_bias <- sum(res_DCSMSE$beta_hat) - sum(beta_true)
      
      # time_2 <- res_DCSMSE$time
      coverage_2  <- abs(sum(res_DCSMSE$beta_hat) - res_DCSMSE$bias - sum(beta_true)) < z * res_DCSMSE$se
      
      print(c(c_h, c_h_init))
      beta_pooled  <- Smoothed_MSE(x, y, beta_init=beta_0, inference = T, c_h=c_h)
      beta_3_bias <- sum(beta_pooled$beta_hat) - sum(beta_true)
      # time_3 <- beta_pooled$time
      coverage_3 <- abs(sum(beta_pooled$beta_hat) - sum(beta_true) - beta_pooled$bias) < z * beta_pooled$se
      
      results_mSMSE <- mSMSE(x, y, local_sizes, v0=rep(1, p), minIterations=5, 
                             lambda_opt = F, c_h=c_h, c_h_init=c_h_init,
                             verbose = F)
      # sum(results_mSMSE$beta_0 - beta_true)
      beta_4_bias <- apply(results_mSMSE$beta_each_iteration[ ,1:5], 2, sum) - sum(beta_true)
      coverage_4 <- apply(as.matrix(results_mSMSE$interval_each_iteration[ ,1:5]), 2, function(x) (sum(beta_true) < x[2]) & (sum(beta_true) > x[1]))
      
      cr <- c(coverage_2, coverage_3, coverage_4)
      methods_names <- c("AvgSMSE", "pooled", paste("round", 1:5, sep=""))
      names(cr) <- c(paste(methods_names, "_cr_", c_h, "_", c_h_init, sep=""))
      res_cr_each_rep <- c(res_cr_each_rep, cr)
      
      bias <- c(beta_2_bias, beta_3_bias, beta_4_bias)
      names(bias) <- c(paste(methods_names, "_bias_", c_h, "_", c_h_init, sep=""))
      res_bias_each_rep <- c(res_bias_each_rep, bias)
    }
  }
  
  return(c(res_cr_each_rep, res_bias_each_rep))
  #res_cr_tune <- rbind(res_cr_tune, res_cr_each_rep)
  #res_bias_tune <- rbind(res_bias_tune, res_bias_each_rep)
  
}

stopCluster(cl)



m <- 1000
n <- floor(m ^ (logmn))
L <- floor(n / m)
local_sizes <- rep(m, L)
data <- data_generate(n, p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
x <- data$x[1:sum(local_sizes), ]
y <- data$y[1:sum(local_sizes)]

results_mSMSE <- mSMSE(x, y, local_sizes, v0=rep(1, p), minIterations=5, maxIterations=6,
                       lambda_opt = F, c_h=1.5, c_h_init=c_h_init_seq,
                       verbose = T)
beta_4_bias <- apply(results_mSMSE$beta_each_iteration[ ,1:5], 2, sum) - sum(beta_true)
coverage_4 <- apply(as.matrix(results_mSMSE$interval_each_iteration[ ,1:5]), 2, function(x) (sum(beta_true) < x[2]) & (sum(beta_true) > x[1]))
print(coverage_4)

