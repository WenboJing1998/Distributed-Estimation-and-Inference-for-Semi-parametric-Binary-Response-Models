 ### Sensitivity to lambda, p = 1
library(psych, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(MASS, warn.conflicts = F, quietly = T)


#os <- "C://Users//WENBO JING//Documents//SMSE//p=10//"
#source(paste(os, "MSE_functions.R", sep=""))

source("MSE_functions.R")

p <- 1
lambda_seq <- c(-1, 1, 10, 30)

logmn_seq <- seq(1.5, 1.9, length=3)
m <- 1000

beta_true <- rep(1/sqrt(p), p)

noiseType <- "norm" 
alpha <- 2 # kernel order

results_all <- rep(NA, 5)

set.seed(0)

for(logmn_id in 1:length(logmn_seq))
{
  # logmn <- 1.5
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L) # local size on each machine
  c_h_init_seq <- c(1, 2, 3)
  
  cl <- makeCluster(6)
  registerDoParallel(cl, cores=6)
  results <- foreach(rep=1:200, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
                     .errorhandling = 'remove') %dorng% {
                       
    data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
    x <- data_MSE$x[1:sum(local_sizes), ]
    y <- data_MSE$y[1:sum(local_sizes)]
    
    lambda_seq <- c(1, 10, 30, -1)
    res_each_rep <- c()
    
    for(c_h in seq(1, 1, 0)){
      for(lambda in lambda_seq){

        lambda_opt <- ifelse(lambda==-1, T, F)
        lambda_h <- ifelse(lambda==-1, 1, lambda)

        results_mSMSE <- mSMSE(x, y, local_sizes, v0=rep(1, p), 
                               minIterations=3, maxIterations=10, lambda_h=lambda_h,
                                lambda_opt = lambda_opt, c_h=c_h, c_h_init=c_h_init_seq,
                                verbose = T)
        coverage_4 <- apply(as.matrix(results_mSMSE$interval_each_iteration[ ,1:5]), 2, function(x) (sum(beta_true) < x[2]) & (sum(beta_true) > x[1]))
        lambda_h <-  results_mSMSE$lambda_h
        res_each_rep <- c(res_each_rep, coverage_4[length(coverage_4)])
        if(lambda==-1) res_each_rep <- c(res_each_rep, lambda_h)
      }
    }
    names(res_each_rep) <- c(paste("cr_lambda=", lambda_seq, sep=""), "lambda_star")
    return(res_each_rep)
  }
  
  stopCluster(cl)
  
  results_all <- rbind(results_all, apply(results, 2, mean))
  print(results_all[-1, ])
}


save(results_all, file=paste("sensi_results_", "p=", p, ".Rdata", sep=""))

# Add column names to the results and save.

methods_names <- c(paste("round", 1:5, sep=""))
colnames(results) <- c(paste(methods_names, "_beta" , sep=""), 
                       paste(methods_names, "_coverage", sep=""),
                       "lambda_h")

save(results, file=paste("Results//Sensi_results_", "p=", p, "_lambda=", lambda, "_rep=", repID, ".Rdata", sep=""))


################ results #################

lambda <- 20

os <- "C://Users//WENBO JING//Documents//SMSE//Sensitivity//Results_p=1//"
p <- 1

logmn_seq <- seq(1.35, 1.75, length=5)

results_summary <-  matrix(NA, ncol = 16,
                           nrow = length(logmn_seq),
                           dimnames=list(as.character(round(logmn_seq, 2)),
                                         c(
                                           "bias-1","bias-2","bias-3","bias-4","bias-5",
                                           "var-1","var-2","var-3","var-4","var-5",
                                           "coverage.rate-1","coverage.rate-2","coverage.rate-3","coverage.rate-4","coverage.rate-5",
                                           "lambda_h"

                                         )))


beta_true <- rep(1/sqrt(p), p)
for (logmn_id in 1:length(logmn_seq)) {
  results_time <- rep(NA, 11)
  for (repID in 1:150) {
    load(paste(os, "Sensi_results_", "p=", p, "_lambda=", lambda, "_rep=", repID, ".Rdata", sep=""))
    results_time <- rbind(results_time, results[logmn_id, ])

  }
  results_time <- results_time[-1, ]
  results_summary[logmn_id, c(1:5)] <- apply(results_time[,1:5], 2, mean) - sum(beta_true)
  results_summary[logmn_id, c(6:10)] <- apply(results_time[,1:5], 2, function(x) mean((x-mean(x))^2))
  results_summary[logmn_id, c(11:15)] <- apply(results_time[ ,6:10], 2, mean)

}

print_table <- function(results_total){


  ## Bias Var Cr

  out1 <- c()
  out2 <- c()
  for(name in rownames(results_total)){
    record <- results_total[name, ]
    row1 <- round(record[11:15], 2)

    out1 <- paste(out1, paste(name, paste(row1, collapse=" & "), sep=" & "), "\\\\", "\n" , sep=" ")

  }

  cat(out1)

  return(1)
}

print_table(results_summary)





# 
# 
# 
