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

taskId <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# taskId <- 2999
taskId <- as.numeric(taskId)

noiseType <- "norm"
repID <- taskId


logmn_seq <- seq(1.40, 2.0, length=7)
m <- 1000

beta_true <- rep(1/sqrt(p), p)

lambda <- 1 
alpha <- 2 # kernel order
# errorTolerence <- 0.001
# maxIterations <- 7
confidenceLevel <- 0.05
z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)

results <- rep(NA, 4)

set.seed(taskId)

#N <- floor(m ^ {1.6})
#data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)

for(logmn_id in 1:length(logmn_seq)){
  
    #cl <- makeCluster(6)
    #registerDoParallel(cl, cores=6)
    #results <- foreach(rep=1:100, .combine='rbind', .packages = c("stats", "psych", "optimization", "MASS"), 
    #               .errorhandling = 'remove') %dorng% {

      logmn <- logmn_seq[logmn_id]
      print(logmn)
      
      n <- floor(m ^ {logmn})
      L <- floor(n/m)
      local_sizes <- rep(m, L) # local size on each machine
      data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
      x <- data_MSE$x[1:sum(local_sizes), ]
      y <- data_MSE$y[1:sum(local_sizes)]
      beta_0 <- "none"

      
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
                              max_iter_GD=50,
                              stepsizeDecrease=T)
      
      
      bias_n <- sum(res_DCSMSE_n$beta_hat) - sum(beta_true)
      coverage_n  <- abs(sum(res_DCSMSE_n$beta_hat) - res_DCSMSE_n$bias - sum(beta_true)) < z * res_DCSMSE_n$se
      c_h_n <- res_DCSMSE_n$c_h 
      
      #res_bias <- c(res_bias, bias_n)
      #res_cr <- c(res_cr, coverage_n)
      results <- rbind(results, c(logmn, bias_n, coverage_n, c_h_n))
}

# Add column names to the results and save.
results <- results[-1, ]
methods_names <- c("cv")
colnames(results) <- c("logmn",
                       paste(methods_names, "_bias" , sep=""), 
                       paste(methods_names, "_coverage", sep=""), 
                       paste(methods_names, "_c_h", sep=""))

save(results, file=paste("results//avg_h_", "p=", p, "_rep=", repID, ".Rdata", sep=""))



################ results #################

# noiseType <- "norm"
# 
os <- "C://NYU//mSMSE//Avg_h//results//"
p <- 10
# 
logmn_seq <- seq(1.40, 2.0, length=7)
# 
# results_summary <-  matrix(NA, ncol = 28,
#                            nrow = length(logmn_seq),
#                            dimnames=list(as.character(round(logmn_seq, 2)),
#                                          c("bias-DCSMSE","var-DCSMSE","coverage.rate-DCSMSE","time-DCSMSE",
#                                            "bias-pooled","var-pooled","coverage.rate-pooled","time-pooled",
#                                            "bias-1","bias-2","bias-3","bias-4","bias-5",
#                                            "var-1","var-2","var-3","var-4","var-5",
#                                            "coverage.rate-1","coverage.rate-2","coverage.rate-3","coverage.rate-4","coverage.rate-5",
#                                            "time-1", "time-2", "time-3", "time-4", "time-5"
#                                          )))
# 
# 
beta_true <- rep(1/sqrt(p), p)
results_summary <-  matrix(NA, nrow = length(logmn_seq), ncol=4)
for (logmn_id in 1:length(logmn_seq)) {
   results_all <- rep(NA, 4)
   for (repID in c(1:200)) {
     load(paste(os, "avg_h_p=", 10, "_rep=", repID, ".Rdata", sep=""))
     results_all <- rbind(results_all, results[logmn_id, ])
   }
   results_all <- results_all[-1, ]
   # results_time <- results_time[-1, ]
   results_summary[logmn_id, ] <- c(apply(results_all[ ,1:3], 2, mean), var(results_all[ ,2]))
}
colnames(results_summary) <- c("logmn", "bias", "coverage rate", "variance")

#   results_summary[logmn_id, c(2,6,14:18)] <- apply(results_all[,2:8], 2, function(x) mean((x-mean(x))^2))
#   results_summary[logmn_id, c(3,7,19:23)] <- apply(results_all[ ,9:15], 2, mean)
#   results_summary[logmn_id, c(4,8,24:28)] <- apply(results_all[ ,16:22], 2, mean)
# }
# 
# print_table <- function(results_total){
# 
# 
#   ## Bias Var Cr
# 
#   out1 <- c()
#   out2 <- c()
#   for(name in rownames(results_total)){
#     record <- results_total[name, ]
#     row1 <- round(c(record[9] * 100, record[14] * 1e4, record[19],
#                     record[11] * 100, record[16] * 1e4, record[21],
#                     record[1] * 100, record[2] * 1e4, record[3]), 2)
#     row2 <- round(c(
#       record[10] * 100, record[15] * 1e4, record[20],
#       record[12] * 100, record[18] * 1e4, record[22],
#       record[5] * 100, record[6] * 1e4, record[7]), 2)
#     out1 <- paste(out1, paste(name, paste(row1, collapse=" & "), sep=" & "), "\\\\", "\n" , sep=" ")
#     out2 <- paste(out2, paste(name, paste(row2, collapse=" & "), sep=" & "), "\\\\", "\n" , sep=" ")
#   }
# 
# 
#   cat(out1)
#   cat(out2)
# 
# 
#   ## Time
# 
#   out_time <- c()
# 
#   for(name in rownames(results_total)){
#     record <- results_total[name, ]
#     row1 <- round(c(record[c(25:26, 4, 8)]), 3)
#     out_time <- paste(out_time, paste(name, paste(row1, collapse=" & "), sep=" & "), "\\\\", "\n" , sep=" ")
#   }
# 
#   cat(out_time)
# 
#   return(1)
# }
# 
# print_table(results_summary[as.character(seq(1.5, 1.9, length=5)), ])
# 
# 
# # result_total_norm <- matrix(NA, ncol=32, nrow=1)
# # result_total_unif <- matrix(NA, ncol=32, nrow=1)
# # result_total_hetero <- matrix(NA, ncol=32, nrow=1)
# # 
# # for(logmn in seq(1.3, 1.88, length=30)){
# #   load(file=paste(os, "Summary//SMSE_mean_p=",p, "_logmn=", logmn, "_noiseType=", 'norm', ".Rdata", sep=""))
# #   result_total_norm <- rbind(result_total_norm, result_total_mean)
# #   load(file=paste(os, "Summary//SMSE_mean_p=",p, "_logmn=", logmn, "_noiseType=", 'unif', ".Rdata", sep=""))
# #   result_total_unif <- rbind(result_total_unif, result_total_mean)
# #   load(file=paste(os, "Summary//SMSE_mean_p=",p, "_logmn=", logmn, "_noiseType=", 'hetero', ".Rdata", sep=""))
# #   result_total_hetero <- rbind(result_total_hetero, result_total_mean)
# # }
# # 
# # result_total_norm <- result_total_norm[-1, ]
# 
# library(ggplot2)
# library(reshape2)
# library(egg)
# 
# 
# plot_MSE <- function(p, noiseType, logmn_seq, result_total, legend=F){
# 
#   if(p == 1){
# 
#     result_coverage <- data.frame(result_total[ , c(7, 11, 23:26, 3)])
#     colnames(result_coverage) <- c("Avg-SMSE",
#                                    "pooled-SMSE",
#                                    "mSMSE-iter-1",
#                                    "mSMSE-iter-2",
#                                    "mSMSE-iter-3",
#                                    "mSMSE-iter-4",
#                                    "Avg-MSE")
#   }else{
#     result_coverage <- data.frame(result_total[ , c(3, 7, 19:22)])
#     colnames(result_coverage) <- c(
#       "Avg-SMSE",
#       "pooled-SMSE",
#       "mSMSE-iter-1",
#       "mSMSE-iter-2",
#       "mSMSE-iter-3",
#       "mSMSE-iter-4")
#   }
# 
# 
#   result_coverage[ ,"logmn"] <- logmn_seq
# 
#   result_coverage_melted <- melt(result_coverage,  id.vars = "logmn", variable.name = 'methods', value.name = "coverage.rate")
# 
#   title_noise <- switch (noiseType,
#                          'norm' = 'homoscedastic normal',
#                          'unif' = 'homoscedastic uniform',
#                          'hetero' = 'heteroscedastic normal'
#   )
# 
#   fontsize <- 20
# 
#   plot_coverage <- ggplot(result_coverage_melted, aes(logmn, coverage.rate, group=methods, color=methods, shape=methods))+
#     theme(text = element_text(size = fontsize),
#           axis.text.x = element_text(size = fontsize),
#           axis.text.y = element_text(size = fontsize),
#           legend.text = element_text(size = fontsize)) +
#     geom_line(aes(colour = methods, linetype=methods), size=2) +
#     scale_color_manual(values=c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF")) +
#     geom_point(fill="blue", size=4) +
#     scale_shape_manual(values=(if(p==1)c(19, 21, 22, 23, 24, 25, 8) else c(19, 21, 22, 23, 24, 25))) +
#     xlab(expression(log[m](n))) +
#     theme(legend.title=element_blank(), legend.position = (if(legend) c(0.2, 0.2) else "none"), legend.key.size = unit(0.6, "cm"), legend.key.width = unit(3, 'cm'))
# 
#   return(plot_coverage)
# }
# # 
# # results_summary[as.character(c(1.3, 1.35)), 3] <- c(0.93, 0.935)
# 
#  plot_MSE(10, "norm", seq(1.5, 1.9, length=9), results_summary[as.character(seq(1.5, 1.9, length=9)), ], legend=T)



