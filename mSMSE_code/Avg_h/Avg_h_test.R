library(psych, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(foreach, warn.conflicts = F, quietly = T)
library(doParallel, warn.conflicts = F, quietly = T)
library(doRNG, warn.conflicts = F, quietly = T)
library(MASS, warn.conflicts = F, quietly = T)


#os <- "C://Users//WENBO JING//Documents//SMSE//p=10//"
#source(paste(os, "MSE_functions.R", sep=""))

source("MSE_functions.R")

p <- 10
noise_seq <- c("norm", "unif", "hetero")

taskId <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# taskId <- 2999
taskId <- as.numeric(taskId)
noiseID <- taskId %% length(noise_seq)
if(noiseID == 0) noiseID <- length(noise_seq) 

repID <- (taskId - noiseID) / length(noise_seq) + 1


logmn_seq <- seq(1.50, 1.9, length=5)
m <- 1000
noiseType <- "norm"
beta_true <- rep(1/sqrt(p), p)

lambda <- 1 
alpha <- 2 # kernel order
errorTolerence <- 0.001
maxIterations <- 7
confidenceLevel <- 0.05
z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
c_h <- 1


results_all <- rep(NA, 19)

set.seed(1)

#N <- floor(m ^ {1.6})
#data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
for(logmn_id in 1:length(logmn_seq)){
  
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L) # local size on each machine
  c_h_seq <- c(0.1, 1, 2, 3, 4, 4.5)
  #lambda_h_seq <- c(1, 10, 20, 30) 
  beta_true <- rep(1/sqrt(p), p)
  #beta_true[seq(2, p, 2)] <- -1/sqrt(p)
  noiseType <- "norm"
  confidenceLevel <- 0.05
  z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
  
  cl <- makeCluster(6)
  registerDoParallel(cl, cores=6)
  results <- foreach(rep=1:100, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                   .errorhandling = 'remove') %dorng% {

      data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
      x <- data_MSE$x[1:sum(local_sizes), ]
      y <- data_MSE$y[1:sum(local_sizes)]
      beta_0 <- "none"
      res_bias <- c()
      res_cr <- c()
      
      # method 2: AvgSMSE, h=h_n, h_m, \sqrt{h_nh_m}
      for(c_h in c_h_seq){
        
        #lambda_opt = ifelse(lambda_h==-1, T, F)
        #lambda_h = ifelse(lambda_h==-1, 1, lambda_h)
        
        res_DCSMSE_n <- DC_SMSE(x, y, beta_init=beta_0, 
                                c_h=c_h,
                                lambda_h=1, 
                                lambda_opt=F, 
                                h_option="n", 
                                local_sizes=local_sizes,
                                max_iter_GD=20,
                                stepsizeDecrease = T)
      #res_DCSMSE_m <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
      #                        h_option="m", local_sizes=local_sizes)
      #res_DCSMSE_ga <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
      #                        h_option="GA", local_sizes=local_sizes)
      
      bias_n <- sum(res_DCSMSE_n$beta_hat) - sum(beta_true)
      coverage_n  <- abs(sum(res_DCSMSE_n$beta_hat) - res_DCSMSE_n$bias - sum(beta_true)) < z * res_DCSMSE_n$se
      
      #bias_m <- sum(res_DCSMSE_m$beta_hat) - sum(beta_true)
      #coverage_m  <- abs(sum(res_DCSMSE_m$beta_hat) - res_DCSMSE_m$bias - sum(beta_true)) < z * res_DCSMSE_m$se
      
      #bias_ga <- sum(res_DCSMSE_ga$beta_hat) - sum(beta_true)
      #coverage_ga  <- abs(sum(res_DCSMSE_ga$beta_hat) - res_DCSMSE_ga$bias - sum(beta_true)) < z * res_DCSMSE_ga$se
      res_bias <- c(res_bias, bias_n)
      res_cr <- c(res_cr, coverage_n)
      
      }
      
      # return(c(logmn, bias_n, bias_m, bias_ga, coverage_n, coverage_m, coverage_ga))
      return(c(logmn, res_bias, res_cr))
   }
  # methods_names <- c("n", "m", "ga")
 
  results_all <- rbind(results_all, c(apply(results, 2, mean), apply(results[, 2:7], 2, var)))
  methods_names <- as.character(c_h_seq)
  colnames(results_all) <- c("logmn",
                      paste(methods_names, "_bias" , sep=""), 
                      paste(methods_names, "_coverage", sep=""),
                      paste(methods_names, "_var" , sep=""))
  print(results_all[-1, ])
}

save(results_all, file="results_all_stepdecrease.RData")


results_all_2 <- rep(NA, 13)

set.seed(1)

#N <- floor(m ^ {1.6})
#data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
for(logmn_id in 1:length(logmn_seq)){
  
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L) # local size on each machine
  c_h_seq <- c(0.1, 1.5, 3, 4.5)
  #lambda_h_seq <- c(1, 10, 20, 30) 
  beta_true <- rep(1/sqrt(p), p)
  #beta_true[seq(2, p, 2)] <- -1/sqrt(p)
  noiseType <- "norm"
  confidenceLevel <- 0.05
  z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
  
  cl <- makeCluster(6)
  registerDoParallel(cl, cores=6)
  results <- foreach(rep=1:100, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                     .errorhandling = 'remove') %dorng% {
                       
                       data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
                       x <- data_MSE$x[1:sum(local_sizes), ]
                       y <- data_MSE$y[1:sum(local_sizes)]
                       beta_0 <- "none"
                       res_bias <- c()
                       res_cr <- c()
                       
                       # method 2: AvgSMSE, h=h_n, h_m, \sqrt{h_nh_m}
                       for(c_h in c_h_seq){
                         
                         #lambda_opt = ifelse(lambda_h==-1, T, F)
                         #lambda_h = ifelse(lambda_h==-1, 1, lambda_h)
                         
                         res_DCSMSE_n <- DC_SMSE(x, y, beta_init=beta_0, 
                                                 c_h=c_h,
                                                 lambda_h=1, 
                                                 lambda_opt=F, 
                                                 h_option="n", 
                                                 local_sizes=local_sizes,
                                                 max_iter_GD=10,
                                                 stepsizeDecrease = T)
                         #res_DCSMSE_m <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                         #                        h_option="m", local_sizes=local_sizes)
                         #res_DCSMSE_ga <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                         #                        h_option="GA", local_sizes=local_sizes)
                         
                         bias_n <- sum(res_DCSMSE_n$beta_hat) - sum(beta_true)
                         coverage_n  <- abs(sum(res_DCSMSE_n$beta_hat) - res_DCSMSE_n$bias - sum(beta_true)) < z * res_DCSMSE_n$se
                         
                         #bias_m <- sum(res_DCSMSE_m$beta_hat) - sum(beta_true)
                         #coverage_m  <- abs(sum(res_DCSMSE_m$beta_hat) - res_DCSMSE_m$bias - sum(beta_true)) < z * res_DCSMSE_m$se
                         
                         #bias_ga <- sum(res_DCSMSE_ga$beta_hat) - sum(beta_true)
                         #coverage_ga  <- abs(sum(res_DCSMSE_ga$beta_hat) - res_DCSMSE_ga$bias - sum(beta_true)) < z * res_DCSMSE_ga$se
                         res_bias <- c(res_bias, bias_n)
                         res_cr <- c(res_cr, coverage_n)
                         
                       }
                       
                       # return(c(logmn, bias_n, bias_m, bias_ga, coverage_n, coverage_m, coverage_ga))
                       return(c(logmn, res_bias, res_cr))
                     }
  # methods_names <- c("n", "m", "ga")
  
  results_all_2 <- rbind(results_all_2, c(apply(results, 2, mean), apply(results[, 2:5], 2, var)))
  methods_names <- as.character(c_h_seq)
  colnames(results_all_2) <- c("logmn",
                             paste(methods_names, "_bias" , sep=""), 
                             paste(methods_names, "_coverage", sep=""),
                             paste(methods_names, "_var" , sep=""))
  print(results_all_2[-1, ])
}

save(results_all_2, file="results_all_2.RData")

results_all_3 <- rep(NA, 19)

set.seed(2)

#N <- floor(m ^ {1.6})
#data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
for(logmn_id in 1:length(logmn_seq)){
  
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L) # local size on each machine
  c_h_seq <- c(0.1, 0.2, 0.4, 0.8, 1.6, 3.2)
  #lambda_h_seq <- c(1, 10, 20, 30) 
  beta_true <- rep(1/sqrt(p), p)
  #beta_true[seq(2, p, 2)] <- -1/sqrt(p)
  noiseType <- "norm"
  confidenceLevel <- 0.05
  z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
  
  cl <- makeCluster(6)
  registerDoParallel(cl, cores=6)
  results <- foreach(rep=1:30, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                     .errorhandling = 'remove') %dorng% {
                       
                       data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
                       x <- data_MSE$x[1:sum(local_sizes), ]
                       y <- data_MSE$y[1:sum(local_sizes)]
                       beta_0 <- "none"
                       res_bias <- c()
                       res_cr <- c()
                       
                       # method 2: AvgSMSE, h=h_n, h_m, \sqrt{h_nh_m}
                       for(c_h in c_h_seq){
                         
                         #lambda_opt = ifelse(lambda_h==-1, T, F)
                         #lambda_h = ifelse(lambda_h==-1, 1, lambda_h)
                         
                         res_DCSMSE_n <- DC_SMSE(x, y, beta_init=beta_0, 
                                                 c_h=c_h,
                                                 lambda_h=1, 
                                                 lambda_opt=F, 
                                                 h_option="n", 
                                                 local_sizes=local_sizes,
                                                 max_iter_GD=10,
                                                 stepsizeDecrease=T)
                         #res_DCSMSE_m <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                         #                        h_option="m", local_sizes=local_sizes)
                         #res_DCSMSE_ga <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                         #                        h_option="GA", local_sizes=local_sizes)
                         
                         bias_n <- sum(res_DCSMSE_n$beta_hat) - sum(beta_true)
                         coverage_n  <- abs(sum(res_DCSMSE_n$beta_hat) - res_DCSMSE_n$bias - sum(beta_true)) < z * res_DCSMSE_n$se
                         
                         #bias_m <- sum(res_DCSMSE_m$beta_hat) - sum(beta_true)
                         #coverage_m  <- abs(sum(res_DCSMSE_m$beta_hat) - res_DCSMSE_m$bias - sum(beta_true)) < z * res_DCSMSE_m$se
                         
                         #bias_ga <- sum(res_DCSMSE_ga$beta_hat) - sum(beta_true)
                         #coverage_ga  <- abs(sum(res_DCSMSE_ga$beta_hat) - res_DCSMSE_ga$bias - sum(beta_true)) < z * res_DCSMSE_ga$se
                         res_bias <- c(res_bias, bias_n)
                         res_cr <- c(res_cr, coverage_n)
                         
                       }
                       
                       # return(c(logmn, bias_n, bias_m, bias_ga, coverage_n, coverage_m, coverage_ga))
                       return(c(logmn, res_bias, res_cr))
                     }
  # methods_names <- c("n", "m", "ga")
  
  results_all_3 <- rbind(results_all_3, c(apply(results, 2, mean), apply(results[, 2:7], 2, var)))
  methods_names <- as.character(c_h_seq)
  colnames(results_all_3) <- c("logmn",
                               paste(methods_names, "_bias" , sep=""), 
                               paste(methods_names, "_coverage", sep=""),
                               paste(methods_names, "_var" , sep=""))
  print(results_all_3[-1, ])
}

save(results_all_3, file="results_all_3.RData")

results_all_cv <- rep(NA, 5)
p <- 1
logmn_seq <-  seq(1.40, 2.0, length=7)
for(logmn_id in c(2, 3, length(logmn_seq))){
  
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L)
  beta_true <- rep(1/sqrt(p), p)

  results_cv <- foreach(rep=1:100, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                   .errorhandling = 'remove') %dorng% {
                     
                     data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
                     x <- data_MSE$x[1:sum(local_sizes), ]
                     y <- data_MSE$y[1:sum(local_sizes)]
                     beta_0 <- "none"
                     res_bias <- c()
                     res_cr <- c()
                     
                     # method 2: AvgSMSE, h=h_n, h_m, \sqrt{h_nh_m}

                      res_DCSMSE_n <- DC_SMSE(x, y, 
                                              beta_init=beta_0,
                                              cv=T,
                                              c_h=1,
                                              c_h_seq=c(0.01, 0.05, 0.1, 0.5, 1, 5, 10), 
                                              lambda_h=1, 
                                              lambda_opt=F, 
                                              h_option="n", 
                                              local_sizes=local_sizes,
                                              max_iter_GD=100,
                                              stepsizeDecrease=T)
        
                       
                       bias_n <- sum(res_DCSMSE_n$beta_hat) - sum(beta_true)
                       coverage_n  <- abs(sum(res_DCSMSE_n$beta_hat) - res_DCSMSE_n$bias - sum(beta_true)) < z * res_DCSMSE_n$se
                       c_h_n <- res_DCSMSE_n$c_h 
                       
                       res_bias <- c(res_bias, bias_n)
                       res_cr <- c(res_cr, coverage_n)
                       
                       return(c(logmn, res_bias, res_cr, c_h_n))
                       
                     }
                     
                     # return(c(logmn, bias_n, bias_m, bias_ga, coverage_n, coverage_m, coverage_ga))
                     
    results_all_cv <- rbind(results_all_cv, c(apply(results_cv, 2, mean), var(results_cv[ ,2])))
    # colnames(results_all_cv) <- 
    print(results_all_cv[-1, ])

}
results_cv[1:100, ]



results_all_m <- rep(NA, 7)

set.seed(1)

#N <- floor(m ^ {1.6})
#data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
for(logmn_id in 1:length(logmn_seq)){
  
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L) # local size on each machine
  c_h_seq <- c(1)
  #lambda_h_seq <- c(1, 10, 20, 30) 
  beta_true <- rep(1/sqrt(p), p)
  #beta_true[seq(2, p, 2)] <- -1/sqrt(p)
  noiseType <- "norm"
  confidenceLevel <- 0.05
  z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
  
  cl <- makeCluster(6)
  registerDoParallel(cl, cores=6)
  results <- foreach(rep=1:100, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                     .errorhandling = 'remove') %dorng% {
                       
                       data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
                       x <- data_MSE$x[1:sum(local_sizes), ]
                       y <- data_MSE$y[1:sum(local_sizes)]
                       beta_0 <- "none"
                       res_bias <- c()
                       res_cr <- c()
                       
                       # method 2: AvgSMSE, h=h_n, h_m, \sqrt{h_nh_m}
                       #for(c_h in c_h_seq){
                         
                         #lambda_opt = ifelse(lambda_h==-1, T, F)
                         #lambda_h = ifelse(lambda_h==-1, 1, lambda_h)
                         
                         res_DCSMSE_n <- DC_SMSE(x, y, beta_init=beta_0, 
                                                 c_h=3.2,
                                                 lambda_h=1, 
                                                 lambda_opt=F, 
                                                 h_option="n", 
                                                 local_sizes=local_sizes)
                         res_DCSMSE_m <- DC_SMSE(x, y, beta_init=beta_0, c_h=1,lambda_h=1,  lambda_opt=F, 
                                                 h_option="m", local_sizes=local_sizes)
                         #res_DCSMSE_ga <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                         #                        h_option="GA", local_sizes=local_sizes)
                         
                         bias_n <- sum(res_DCSMSE_n$beta_hat) - sum(beta_true)
                         coverage_n  <- abs(sum(res_DCSMSE_n$beta_hat) - res_DCSMSE_n$bias - sum(beta_true)) < z * res_DCSMSE_n$se
                         
                         bias_m <- sum(res_DCSMSE_m$beta_hat) - sum(beta_true)
                         coverage_m  <- abs(sum(res_DCSMSE_m$beta_hat) - res_DCSMSE_m$bias - sum(beta_true)) < z * res_DCSMSE_m$se
                         
                         #bias_ga <- sum(res_DCSMSE_ga$beta_hat) - sum(beta_true)
                         #coverage_ga  <- abs(sum(res_DCSMSE_ga$beta_hat) - res_DCSMSE_ga$bias - sum(beta_true)) < z * res_DCSMSE_ga$se
                         res_bias <- c(res_bias, bias_n, bias_m)
                         res_cr <- c(res_cr, coverage_n, coverage_m)
                         
                      # }
                       
                       # return(c(logmn, bias_n, bias_m, bias_ga, coverage_n, coverage_m, coverage_ga))
                       return(c(logmn, res_bias, res_cr))
                     }
  # methods_names <- c("n", "m", "ga")
  
  results_all_m <- rbind(results_all_m, c(apply(results, 2, mean), apply(results[, c(2, 3)], 2, var)))
  methods_names <- c("n_3.5", "m")
  colnames(results_all_m) <- c("logmn",
                             paste(methods_names, "_bias" , sep=""), 
                             paste(methods_names, "_coverage", sep=""),
                             paste(methods_names, "_var" , sep=""))
  print(results_all_m[-1, ])
}



results_all_0 <- rep(NA, 16)

set.seed(1)

#N <- floor(m ^ {1.6})
#data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
for(logmn_id in 1:length(logmn_seq)){
  
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L) # local size on each machine
  c_h_seq <- c(0.5, 1, 2, 3, 4)
  #lambda_h_seq <- c(1, 10, 20, 30) 
  beta_true <- rep(1/sqrt(p), p)
  beta_true[seq(2, p, 2)] <- -1/sqrt(p)
  noiseType <- "norm"
  confidenceLevel <- 0.05
  z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
  
  cl <- makeCluster(6)
  registerDoParallel(cl, cores=6)
  results <- foreach(rep=1:100, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                     .errorhandling = 'remove') %dorng% {
                       
                       data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
                       x <- data_MSE$x[1:sum(local_sizes), ]
                       y <- data_MSE$y[1:sum(local_sizes)]
                       beta_0 <- "none"
                       res_bias <- c()
                       res_cr <- c()
                       
                       # method 2: AvgSMSE, h=h_n, h_m, \sqrt{h_nh_m}
                       for(c_h in c_h_seq){
                         
                         #lambda_opt = ifelse(lambda_h==-1, T, F)
                         #lambda_h = ifelse(lambda_h==-1, 1, lambda_h)
                         
                         res_DCSMSE_n <- DC_SMSE(x, y, beta_init=beta_0, 
                                                 c_h=c_h,
                                                 lambda_h=1, 
                                                 lambda_opt=F, 
                                                 h_option="n", 
                                                 local_sizes=local_sizes)
                         #res_DCSMSE_m <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                         #                        h_option="m", local_sizes=local_sizes)
                         #res_DCSMSE_ga <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                         #                        h_option="GA", local_sizes=local_sizes)
                         
                         bias_n <- sum(res_DCSMSE_n$beta_hat) - sum(beta_true)
                         coverage_n  <- abs(sum(res_DCSMSE_n$beta_hat) - res_DCSMSE_n$bias - sum(beta_true)) < z * res_DCSMSE_n$se
                         
                         #bias_m <- sum(res_DCSMSE_m$beta_hat) - sum(beta_true)
                         #coverage_m  <- abs(sum(res_DCSMSE_m$beta_hat) - res_DCSMSE_m$bias - sum(beta_true)) < z * res_DCSMSE_m$se
                         
                         #bias_ga <- sum(res_DCSMSE_ga$beta_hat) - sum(beta_true)
                         #coverage_ga  <- abs(sum(res_DCSMSE_ga$beta_hat) - res_DCSMSE_ga$bias - sum(beta_true)) < z * res_DCSMSE_ga$se
                         res_bias <- c(res_bias, bias_n)
                         res_cr <- c(res_cr, coverage_n)
                         
                       }
                       
                       # return(c(logmn, bias_n, bias_m, bias_ga, coverage_n, coverage_m, coverage_ga))
                       return(c(logmn, res_bias, res_cr))
                     }
  # methods_names <- c("n", "m", "ga")
  
  results_all_0 <- rbind(results_all_0, c(apply(results, 2, mean), apply(results[, 2:6], 2, var)))
  methods_names <- as.character(c_h_seq)
  colnames(results_all_0) <- c("logmn",
                               paste(methods_names, "_bias" , sep=""), 
                               paste(methods_names, "_coverage", sep=""),
                               paste(methods_names, "_var" , sep=""))
  print(results_all_0[-1, ])
}





