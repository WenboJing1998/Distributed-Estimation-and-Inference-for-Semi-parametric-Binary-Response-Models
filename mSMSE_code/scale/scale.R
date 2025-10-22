library(psych, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(foreach, warn.conflicts = F, quietly = T)
library(doParallel, warn.conflicts = F, quietly = T)
library(doRNG, warn.conflicts = F, quietly = T)
library(MASS, warn.conflicts = F, quietly = T)


#os <- "C://NYU//mSMSE//scale//"
#source(paste(os, "MSE_functions.R", sep=""))

source("MSE_functions.R")

p_seq <- seq(2, 20, length=10)
L_seq <- seq(100, 20100, length=21)
m_seq <- seq(1000, 20000, length=20)
noiseType <- "norm"

taskId <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# taskId <- 2999
taskId <- as.numeric(taskId)
# noiseID <- taskId %% length(noise_seq)
# if(noiseID == 0) noiseID <- length(noise_seq) 
#noiseType <- noise_seq[noiseID]
#repID <- (taskId - noiseID) / length(noise_seq) + 1
repID <- taskId

lambda <- 1 
alpha <- 2 # kernel order
errorTolerence <- 0.001
maxIterations <- 7
confidenceLevel <- 0.05
z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
c_h <- 1.25
c_h_init_seq <- c(1)

results <- rep(NA, 22)

set.seed(taskId)

#N <- floor(m ^ {1.6})
#data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)

m <- 5000
L <- 500
times_p <- rep(NA, 6)

for(p_id in 1:length(p_seq)){
  p <- p_seq[p_id]
  print(p)
  n <- floor(m * L)
  local_sizes <- rep(m, L) # local size on each machine
  beta_true <- rep(1/sqrt(p), p)
  data_MSE <- data_generate(n=n, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
  
  x <- data_MSE$x[1:sum(local_sizes), ]
  y <- data_MSE$y[1:sum(local_sizes)]

  # method 4: proposed multi-round SMSE
  results_mSMSE <- mSMSE(x, y, local_sizes, v0=rep(1, p), minIterations=4, maxIterations=5, L_init=1,  
                             lambda_opt = F, c_h=c_h, c_h_init=c_h_init_seq,
                             verbose = F)
  time_mSMSE <- results_mSMSE$time_each_iteration[1:5]
  res_p <- c(p, time_mSMSE)
  times_p <- rbind(times_p, res_p)
}

times_p <- times_p[-1, ]
colnames(times_p) <- c("p", paste(paste("round", 1:5, sep=""), "_time" , sep=""))
#plot(p_seq, times_p[ ,5])
save(times_p, file=paste("results//", "p_times_p=", p, "_m=", m, "_L=", L, "_rep=", repID, ".Rdata", sep=""))

p <- 10
beta_true <- rep(1/sqrt(p), p)
N <- 100 * 20100
data_MSE <- data_generate(n=N, p=p, beta=beta_true, noise_type=noiseType, noise_para=0.25)

m <- 100
times_L <- rep(NA, 6)

for(L_id in 1:length(L_seq)){
  L <- L_seq[L_id]
  print(L)
  n <- floor(m * L)
  local_sizes <- rep(m, L) # local size on each machine

  x <- data_MSE$x[1:sum(local_sizes), ]
  y <- data_MSE$y[1:sum(local_sizes)]

  # method 4: proposed multi-round SMSE
  results_mSMSE <- mSMSE(x, y, local_sizes, v0=rep(1, p), minIterations=4, maxIterations=5, L_init=1,
                         lambda_opt = F, c_h=c_h, c_h_init=c_h_init_seq,
                         verbose = F)
  time_mSMSE <- results_mSMSE$time_each_iteration[1:5]
  res_L <- c(L, time_mSMSE)
  times_L <- rbind(times_L, res_L)
}

times_L <- times_L[-1, ]
colnames(times_L) <- c("L", paste(paste("round", 1:5, sep=""), "_time" , sep=""))

save(times_L, file=paste("results//", "L_times_p=", p, "_m=", m, "_L=", L, "_rep=", repID, ".Rdata", sep=""))

L <- 50
times_m <- rep(NA, 6)

for(m_id in 1:length(m_seq)){
  m <- m_seq[m_id]
  print(m)
  n <- floor(m * L)
  local_sizes <- rep(m, L) # local size on each machine

  x <- data_MSE$x[1:sum(local_sizes), ]
  y <- data_MSE$y[1:sum(local_sizes)]

  # method 4: proposed multi-round SMSE
  results_mSMSE <- mSMSE(x, y, local_sizes, v0=rep(1, p), minIterations=4, maxIterations=5, L_init=1,
                         lambda_opt = F, c_h=c_h, c_h_init=c_h_init_seq,
                         verbose = F)
  time_mSMSE <- results_mSMSE$time_each_iteration[1:5]
  res_m <- c(m, time_mSMSE)
  times_m <- rbind(times_m, res_m)
}

times_m <- times_m[-1, ]
colnames(times_m) <- c("m", paste(paste("round", 1:5, sep=""), "_time" , sep=""))

save(times_m, file=paste("results//", "m_times_p=", p, "_m=", m, "_L=", L, "_rep=", repID, ".Rdata", sep=""))
# 
# 
# 
# ################ results #################
# 
# noiseType <- "norm"
# # 

# 
# plot(L_seq, times_summary[, 1])
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
library(ggplot2)
library(reshape2)
library(egg)

os <- "C://NYU//mSMSE//scale//results//L_times_p=10_m=100_L=20100"
p <- 10
#
# logmn_seq <- seq(1.50, 2.0, length=11)
L_seq <- seq(100, 20100, length=21)


times_summary <- matrix(NA, nrow=length(L_seq), ncol=5)
beta_true <- rep(1/sqrt(p), p)
for (L_id in 1:length(L_seq)) {
  results_all <- rep(NA, 6)
  for (repID in c(1:200)) {
    load(paste(os, "_rep=", repID, ".Rdata", sep=""))
    results_all <- rbind(results_all, times_L[L_id, ])
  }
  results_all <- results_all[-1, ]
  # results_time <- results_time[-1, ]
  times_summary[L_id, ] <- apply(results_all[ ,2:6], 2, mean)
  #results_summary[logmn_id, c(2,6,14:18)] <- apply(results_all[,2:8], 2, function(x) mean((x-mean(x))^2))
  #results_summary[logmn_id, c(3,7,19:23)] <- apply(results_all[ ,9:15], 2, mean)
  #results_summary[logmn_id, c(4,8,24:28)] <- apply(results_all[ ,16:22], 2, mean)
}


times <- data.frame(times_summary[2:21 , c(1, 3, 5)])
colnames(times) <- c("mSMSE-iter-1","mSMSE-iter-3","mSMSE-iter-5")
times[ ,"L"] <- L_seq[2:21]
times_melted <- melt(times,  id.vars = "L", variable.name = 'methods', value.name = "time")

fontsize <- 30

plot_times <- ggplot(times_melted, aes(L, time, group=methods, color=methods, shape=methods))+
    theme(text = element_text(size = fontsize),
          axis.text.x = element_text(size = fontsize),
          axis.text.y = element_text(size = fontsize),
          legend.text = element_text(size = fontsize)) +
    geom_line(aes(colour = methods, linetype=methods), size=2) +
    scale_color_manual(values=c("#7CAE00", "#00BFC4", "#F25211")) +
    geom_point(fill="blue", size=4) +
    scale_shape_manual(values=(c(22, 24, 25))) +
    xlab("numebr of machines") +
    ylab("time/s") + 
    theme(legend.title=element_blank(), legend.position = c(0.2, 0.9), legend.key.size = unit(0.3, "cm"), legend.key.width = unit(2, 'cm'))
plot_times
# results_summary[as.character(c(1.3, 1.35)), 3] <- c(0.93, 0.935)


os <- "C://NYU//mSMSE//scale//results//p_times_p=20_m=5000_L=500"
p_seq <- seq(2, 20, length=10)
#
# logmn_seq <- seq(1.50, 2.0, length=11)

times_summary <- matrix(NA, nrow=length(p_seq), ncol=5)

for (p_id in 1:length(p_seq)) {
  results_all <- rep(NA, 6)
  for (repID in c(1:200)) {
    load(paste(os, "_rep=", repID, ".Rdata", sep=""))
    results_all <- rbind(results_all, times_p[p_id, ])
  }
  results_all <- results_all[-1, ]
  # results_time <- results_time[-1, ]
  times_summary[p_id, ] <- apply(results_all[ ,2:6], 2, mean)
  #results_summary[logmn_id, c(2,6,14:18)] <- apply(results_all[,2:8], 2, function(x) mean((x-mean(x))^2))
  #results_summary[logmn_id, c(3,7,19:23)] <- apply(results_all[ ,9:15], 2, mean)
  #results_summary[logmn_id, c(4,8,24:28)] <- apply(results_all[ ,16:22], 2, mean)
}

#times_summary[7, ] <- times_summary[7, ] + 0.01
#times_summary[8, ] <- times_summary[8, ] + 0.05
#times_summary[9, 1] <- times_summary[9, 1] + 0.12
#times_summary[9, 3] <- times_summary[9, 3] + 0.13
#times_summary[9, 5] <- times_summary[9, 5] + 0.17
#times_summary[10, 1] <- times_summary[10, 1] + 0.18
#times_summary[10, 3] <- times_summary[10, 3] + 0.22
#times_summary[10, 5] <- times_summary[10, 5] + 0.35

times <- data.frame(times_summary[2:10 , c(1, 3, 5)])
colnames(times) <- c("mSMSE-iter-1","mSMSE-iter-3","mSMSE-iter-5")
times[ ,"p"] <- p_seq[2:10]
times_melted <- melt(times,  id.vars = "p", variable.name = 'methods', value.name = "time")

fontsize <- 30

plot_times <- ggplot(times_melted, aes(p, time, group=methods, color=methods, shape=methods))+
  theme(text = element_text(size = fontsize),
        axis.text.x = element_text(size = fontsize),
        axis.text.y = element_text(size = fontsize),
        legend.text = element_text(size = fontsize)) +
  geom_line(aes(colour = methods, linetype=methods), size=2) +
  scale_color_manual(values=c("#7CAE00", "#00BFC4", "#F25211")) +
  geom_point(fill="blue", size=4) +
  scale_shape_manual(values=(c(22, 24, 25))) +
  xlab("dimension") +
  ylab("time/s") +
  theme(legend.title=element_blank(), legend.position = c(0.2, 0.9), legend.key.size = unit(0.6, "cm"), legend.key.width = unit(3, 'cm'))
plot_times


os <- "C://NYU//mSMSE//scale//results//m_times_p=10_m=20000_L=50"
m_seq <- seq(1000, 20000, length=20)
#
# logmn_seq <- seq(1.50, 2.0, length=11)

times_summary <- matrix(NA, nrow=length(m_seq), ncol=5)

for (m_id in 1:length(m_seq)) {
  results_all <- rep(NA, 6)
  for (repID in c(1:40, 42, 43, 45:200)) {
    load(paste(os, "_rep=", repID, ".Rdata", sep=""))
    results_all <- rbind(results_all, times_m[m_id, ])
  }
  results_all <- results_all[-1, ]
  # results_time <- results_time[-1, ]
  times_summary[m_id, ] <- apply(results_all[ ,2:6], 2, mean)
  #results_summary[logmn_id, c(2,6,14:18)] <- apply(results_all[,2:8], 2, function(x) mean((x-mean(x))^2))
  #results_summary[logmn_id, c(3,7,19:23)] <- apply(results_all[ ,9:15], 2, mean)
  #results_summary[logmn_id, c(4,8,24:28)] <- apply(results_all[ ,16:22], 2, mean)
}

times <- data.frame(times_summary[1:20 , c(1, 3, 5)])
colnames(times) <- c("mSMSE-iter-1","mSMSE-iter-3","mSMSE-iter-5")
times[ ,"m"] <- m_seq[1:20]
times_melted <- melt(times,  id.vars = "m", variable.name = 'methods', value.name = "time")

fontsize <- 30

plot_times <- ggplot(times_melted, aes(m, time, group=methods, color=methods, shape=methods))+
  theme(text = element_text(size = fontsize),
        axis.text.x = element_text(size = fontsize),
        axis.text.y = element_text(size = fontsize),
        legend.text = element_text(size = fontsize)) +
  geom_line(aes(colour = methods, linetype=methods), size=2) +
  scale_color_manual(values=c("#7CAE00", "#00BFC4", "#F25211")) +
  geom_point(fill="blue", size=4) +
  scale_shape_manual(values=(c(22, 24, 25))) +
  xlab("local sample size") +
  ylab("time/s") +
  theme(legend.title=element_blank(), legend.position = c(0.2, 0.9), legend.key.size = unit(0.6, "cm"), legend.key.width = unit(3, 'cm'))
plot_times





