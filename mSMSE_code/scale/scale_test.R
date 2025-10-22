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
noiseType <- noise_seq[noiseID]
repID <- (taskId - noiseID) / length(noise_seq) + 1


logmn_seq <- seq(1.50, 1.9, length=5)
m <- 1000

beta_true <- rep(1/sqrt(p), p)

lambda <- 1 
alpha <- 2 # kernel order
errorTolerence <- 0.001
maxIterations <- 7
confidenceLevel <- 0.05
z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
c_h <- 1.0

results_all <- rep(NA, 7)

set.seed(0)

#N <- floor(m ^ {1.6})
#data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
for(logmn_id in 1:length(logmn_seq)){
  
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L) # local size on each machine
  c_h <- 1
  beta_true <- rep(1/sqrt(p), p)
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
      
      # method 2: AvgSMSE, h=h_n, h_m, \sqrt{h_nh_m}
      
      res_DCSMSE_n <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                            h_option="n", local_sizes=local_sizes)
      res_DCSMSE_m <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                              h_option="m", local_sizes=local_sizes)
      res_DCSMSE_ga <- DC_SMSE(x, y, beta_init=beta_0, c_h=c_h,
                              h_option="GA", local_sizes=local_sizes)
      
      bias_n <- sum(res_DCSMSE_n$beta_hat) - sum(beta_true)
      coverage_n  <- abs(sum(res_DCSMSE_n$beta_hat) - res_DCSMSE_n$bias - sum(beta_true)) < z * res_DCSMSE_n$se
      
      bias_m <- sum(res_DCSMSE_m$beta_hat) - sum(beta_true)
      coverage_m  <- abs(sum(res_DCSMSE_m$beta_hat) - res_DCSMSE_m$bias - sum(beta_true)) < z * res_DCSMSE_m$se
      
      bias_ga <- sum(res_DCSMSE_ga$beta_hat) - sum(beta_true)
      coverage_ga  <- abs(sum(res_DCSMSE_ga$beta_hat) - res_DCSMSE_ga$bias - sum(beta_true)) < z * res_DCSMSE_ga$se
   
      return(c(logmn, bias_n, bias_m, bias_ga, coverage_n, coverage_m, coverage_ga))
   }
  methods_names <- c("n", "m", "ga")
  names(results) <- c("logmn",
                      paste(methods_names, "_bias" , sep=""), 
                      paste(methods_names, "_coverage", sep=""))
  results_all <- rbind(results_all, apply(results, 2, mean))
  print(results_all[-1, ])

}







