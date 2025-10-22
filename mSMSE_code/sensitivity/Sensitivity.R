 ### Sensitivity to lambda
library(psych, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(MASS, warn.conflicts = F, quietly = T)


#os <- "C://Users//WENBO JING//Documents//SMSE//p=10//"
#source(paste(os, "MSE_functions.R", sep=""))

source("MSE_functions.R")

p <- 10
lambda_seq <- c(1, 10, 30, -1)

# lambda <- -1

taskId <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#taskId <- 798
taskId <- as.numeric(taskId)
lambdaID <- taskId %% length(lambda_seq)
if(lambdaID == 0) lambdaID <- length(lambda_seq) 
lambda <- lambda_seq[lambdaID]
repID <- (taskId - lambdaID) / length(lambda_seq) + 1


logmn_seq <- seq(1.5, 1.9, length=3)
m <- 1000

beta_true <- rep(1/sqrt(p), p)

noiseType <- "norm" 
alpha <- 2 # kernel order
errorTolerence <- 0.001
maxIterations <- 10
confidenceLevel <- 0.05
z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)

lambda_opt <- ifelse(lambda==-1, T, F)
lambda_h <- ifelse(lambda==-1, 1, lambda)

c_h <- 1.25
c_h_init_seq <- c(1, 2, 3)

results <- matrix(NA, nrow=length(logmn_seq), ncol=3)

set.seed(taskId)

# N <- floor(m ^ {1.76})

for(logmn_id in 1:length(logmn_seq)){
  
  logmn <- logmn_seq[logmn_id]
  print(logmn)
  
  n <- floor(m ^ {logmn})
  L <- floor(n/m)
  local_sizes <- rep(m, L) # local size on each machine
  data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)

  x <- data_MSE$x[1:sum(local_sizes), ]
  y <- data_MSE$y[1:sum(local_sizes)]
  
  beta_0 <- "SMSE"

  # method 4: proposed multi-round SMSE
  
  results_mSMSE <- mSMSE(x, y, local_sizes, v0=rep(1, p), 
                         minIterations=3, maxIterations=10, lambda_h=lambda_h,
                         lambda_opt = lambda_opt, c_h=c_h, c_h_init=c_h_init_seq,
                         verbose = T)
  coverage_4 <- apply(as.matrix(results_mSMSE$interval_each_iteration), 2, function(x) (sum(beta_true) < x[2]) & (sum(beta_true) > x[1]))
  lambda_h <-  results_mSMSE$lambda_h
  results[logmn_id, ] <- c(logmn, coverage_4[length(coverage_4)], lambda_h)
}

# Add column names to the results and save.

# methods_names <- c(paste("round", 1:5, sep=""))
#colnames(results) <- c(paste(methods_names, "_beta" , sep=""), 
#                       paste(methods_names, "_coverage", sep=""),
#                       "lambda_h")

colnames(results) <- c("logmn", "coverage", "lambda_h")

save(results, file=paste("results//sensi_results_", "p=", p, "_lambda=", lambda, "_rep=", repID, ".Rdata", sep=""))






################ results #################

lambda <- -1

os <- "C://NYU//mSMSE//Sensitivity//results//"
p <- 10

logmn_seq <- seq(1.5, 1.9, length=3)

# results_summary <-  matrix(NA, ncol = 16,
#                            nrow = length(logmn_seq),
#                            dimnames=list(as.character(round(logmn_seq, 2)),
#                                          c(
#                                            "bias-1","bias-2","bias-3","bias-4","bias-5",
#                                            "var-1","var-2","var-3","var-4","var-5",
#                                            "coverage.rate-1","coverage.rate-2","coverage.rate-3","coverage.rate-4","coverage.rate-5",
#                                            "lambda_h"
# 
#                                          )))
results_summary <-  matrix(NA, ncol = 3,
                           nrow = length(logmn_seq))

beta_true <- rep(1/sqrt(p), p)
for (logmn_id in 1:length(logmn_seq)) {
  results_all <- rep(NA, 3)
  for (repID in 1:200) {
    load(paste(os, "sensi_results_", "p=", p, "_lambda=", lambda, "_rep=", repID, ".Rdata", sep=""))
    results_all <- rbind(results_all, results[logmn_id, ])

  }
  results_all <- results_all[-1, ]
  #results_summary[logmn_id, c(1:5)] <- apply(results_time[,1:5], 2, mean) - sum(beta_true)
  #results_summary[logmn_id, c(6:10)] <- apply(results_time[,1:5], 2, function(x) mean((x-mean(x))^2))
  #results_summary[logmn_id, c(11:15)] <- apply(results_time[ ,6:10], 2, mean)
  results_summary[logmn_id, ] <- apply(results_all, 2, mean)
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
