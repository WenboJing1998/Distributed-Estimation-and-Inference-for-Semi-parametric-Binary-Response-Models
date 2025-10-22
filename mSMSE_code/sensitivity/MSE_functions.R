## data generation
# noise_para: the standard-deviance of the distribution of the noise

data_generate <- function(n, p, beta, noise_type, noise_para){
  X <- rnorm(n, 0, 1)
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
  cov <- (0.5)^exponent
  Z <- mvrnorm(n = n, 
               mu = rep(0, p),  
               Sigma = cov) 
  colnames(Z) <- paste("Z", 1:p)
  XZbeta <- X + Z %*% beta 
  if(p==1){
    epsilon <- switch (noise_type,
                 "unif" = runif(n, min=-sqrt(3)*noise_para, max=sqrt(3)*noise_para),
                 "norm" = rnorm(n, mean=0, sd=noise_para),
                 "hetero" = noise_para * ((1 + 0.5 * ( Z[ ,1]) ^ 2 )) * rnorm(n, 0, 1) / 1.66) # let sd(u) = noise_para
  }else{
    epsilon <- switch (noise_type,
                 "unif" = runif(n, min=-sqrt(3)*noise_para, max=sqrt(3)*noise_para),
                 "norm" = rnorm(n, mean=0, sd=noise_para),
                 "hetero" = noise_para * ((1 + 0.5 * (Z[ ,1] + Z[ ,2]) ^ 2 )) * rnorm(n, 0, 1) / 3.28)  # let sd(u) = noise_para
  }
  y <- 2 * ((XZbeta + epsilon) >= 0) - 1
  list(y=y, x=cbind(X, Z), epsilon=epsilon)
}

#library("MASS")
# p = 10
# data = data_generate(n=5000000, p=10, beta=rep(1/sqrt(p), p), noise_type="norm", noise_para=0.25)
# noise_type = "norm"
# n=5000000
# noise_para = 0.25

#sd(data$epsilon)


# second-order kernel
kerf <- function(v) (abs(v) < 1) * (0.5 + 15/16 * (v - 2/3 * v^3 + 1/5 * v^5)) + 1*(v>=1)
dkerf <- function(v) (abs(v) < 1) * (15/16 * (1 - v^2)^2)
d2kerf <- function(v) (abs(v) < 1) * (15/8 * (1 - v^2)*(-2*v))


# Total search algorithm for p=1

TotalSearch <- function(f, lb, ub, len){
   value_min <- 1e8
   optimizer <- lb
   point <- lb
   while(point < ub){
     # print(point)
     point <- point + len
     fpoint <- f(point)
     if(fpoint < value_min){
       optimizer <- point
       value_min <- fpoint
     }
   }
   return(optimizer)
}


## MSE by Manski, p=1

Ave_MSE <- function(x, y, local_sizes){
  
  p <- ncol(x) - 1
  n <- nrow(x)
  
  if(sum(local_sizes)!=n) stop("Sizes do not match!")
  
  L <- length(local_sizes)
  sizes_cumsum <- c(0, cumsum(local_sizes))
  
  beta_list <- matrix(NA, nrow = p, ncol = L)
  
  time_each_machine <- c()

  for (l in 1:L) {
      
      time_start <- Sys.time()
      
      S_n <- function(beta) mean( - y[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1])] * 
                                    (x[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1]), ] %*% c(1, beta) > 0))
      
      beta_hat <- TotalSearch(S_n, lb=0.5, ub=1.5, len=2e-3)
      
      beta_list[ ,l] <- beta_hat
      
      time_end <- Sys.time()
      
      time_each_machine <- c(time_each_machine,  time_end - time_start)
  }
    


  time_average_start <- Sys.time()
  
  beta_hat <- apply(local_sizes * t(beta_list) / n, 2, sum)
  
  time_average_end <- Sys.time()
  
  time_average <- time_average_end - time_average_start 
  
  time_total <- max(time_each_machine) + time_average
  
  ## Inference
  
  se <- (1/sqrt(L)) * sd(beta_list[1, ])
  
  return(list(beta_hat=beta_hat,
              se=se,
              beta_list=beta_list,
              time=time_total))
}



## Smoothed MSE proposed by Horowitz, solved by Gradient Descent or Path Following method by Feng.(2019)

Smoothed_MSE <- function(x, 
                         y,
                         beta_init="none",
                         c_h=1,
                         kernelOrder=2, 
                         delta=0.1,
                         inference=T){
  
  
  p <- ncol(x) - 1
  n <- nrow(x)
  
  
  alpha <- kernelOrder
  kerf <- function(v) (abs(v) < 1) * (0.5 + 15/16 * (v - 2/3 * v^3 + 1/5 * v^5)) + 1*(v>=1)
  dkerf <- function(v) (abs(v) < 1) * (15/16 * (1 - v^2)^2)
  d2kerf <- function(v) (abs(v) < 1) * (15/8 * (1 - v^2)*(-2*v))
  
  #lambda_h <- lambda_init
  #lambda_star <- lambda_h
  
  h_n <- c_h * (1 / n) ^ {1 / (2 * alpha + 1)}
  h_delta <- c_h * (1 / n) ^ {delta / (2 * alpha + 1)}
  
  if(all(beta_init=="none")){
    beta_init <- rep(0, p)
  }
  
  if(p == 1){
    time_start <- Sys.time()
    res <- GradientDescent(beta_0 = beta_init, x, y, h_n=h_n)
    beta_hat <- res$beta
    #S_n <- function(beta) mean( - y * kerf(x %*% c(1, beta) / h_n))
    #beta_hat <- TotalSearch(S_n, lb=-2, ub=2, len=1e-3)
    time_end <- Sys.time()
    time <- time_end - time_start
  }else{
    time_start <- Sys.time()
    res <- GradientDescent(beta_0 = beta_init, x, y, h_n=h_n)
    beta_hat <- res$beta
    time_end <- Sys.time()
    time <- time_end - time_start
  }
  
  if(!inference)  return(list(beta_hat=beta_hat, 
                              time=time))
  
  
  ## Make inference for sum(beta)
  
  xbeta <- x %*% c(1, beta_hat)
  
  if(p == 1){
    W_hat <- h_delta ^ {-alpha-1} * mean( c( y * dkerf(xbeta / h_delta)) * x[ , -1])
    V_hat <- (1 / h_n^2) * mean( c(- y * d2kerf(xbeta / h_n)) * (x[ ,-1])^2)
    Va_hat <- (1 / h_n) * mean((dkerf(xbeta / h_n))^2 * (x[ ,-1])^2)
  }else{
    W_hat <- h_delta ^ {-alpha-1} * apply(c(y * dkerf(xbeta / h_delta)) * x[ , -1], 2, mean)
    
    V_hat <- matrix(0, nrow=p, ncol=p)
    yd2H <- c(y * d2kerf(xbeta / h_n))
    
    Va_hat <- matrix(0, nrow=p, ncol=p)
    ydH2 <- (dkerf(xbeta / h_n))^2
    
    for (i in 1:n) {
      V_hat <- V_hat + (1 / (n * h_n^2)) * -yd2H[i] * (x[i, -1]) %*% t(x[i, -1])
      Va_hat <- Va_hat + (1 / (n * h_n)) * ydH2[i] * (x[i, -1]) %*% t(x[i, -1])
    }
  }
  
  
  V_inv <- solve(V_hat)
  ViW <- V_inv %*% W_hat
  ViVaVi <- V_inv %*% Va_hat %*% V_inv
  
  se <- 1 / (sqrt(n * h_n)) * sqrt(sum(ViVaVi))
  bias <- h_n ^ (alpha) * sum(V_inv %*% W_hat)

  return(list(beta_hat=beta_hat,
              bias=bias, 
              se=se,
              time=time)) #,
              #lambda_star=lambda_star))
}


## Divide and Conquer SMSE
DC_SMSE <- function(x, 
                    y,
                    local_sizes,
                    beta_init="none",
                    h_option="n",
                    c_h=1,
                    kernelOrder=2, 
                    delta=0.1){
  p <- ncol(x) - 1
  n <- nrow(x)
  
  if(sum(local_sizes)!=n) stop("Sizes do not match!")
  
  L <- length(local_sizes)
  sizes_cumsum <- c(0, cumsum(local_sizes))
  
  alpha <- kernelOrder 

  kerf <- function(v) (abs(v) < 1) * (0.5 + 15/16 * (v - 2/3 * v^3 + 1/5 * v^5)) + 1*(v>=1)
  dkerf <- function(v) (abs(v) < 1) * (15/16 * (1 - v^2)^2)
  d2kerf <- function(v) (abs(v) < 1) * (15/8 * (1 - v^2)*(-2*v))
  
  #lambda_h <- lambda_init
  #lambda_star <- lambda_h
  
  beta_list <- matrix(NA, nrow = p, ncol = L)
  
  time_each_machine <- c()
  
  if(all(beta_init=="none")){
    beta_init <- rep(0, p)
  }
  
  if(p == 1){
    
    if(h_option=="n"){
      h_n <- c_h * (1 / n) ^ {1 / (2 * alpha + 1)}
    }else{
      h_n <- c_h * (1 / local_sizes[l]) ^ {1 / (2 * alpha + 1)}
    }
    
    for (l in 1:L) {
      time_start <- Sys.time()
      res <- GradientDescent(beta_0 = beta_init,
                             x=x[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1]), ],
                             y=y[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1])],
                             h_n=h_n)
      time_end <- Sys.time()
      beta_hat <- res$beta
      #S_n <- function(beta) mean( - y[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1])] * 
                                      #kerf(x[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1]), ] %*% c(1, beta) / h_n))
      
      #beta_hat <- TotalSearch(S_n, lb=-2, ub=2, len=1e-3)
      beta_list[ ,l] <- beta_hat
      time_each_machine <- c(time_each_machine,  time_end - time_start)
    }
      
      
    
  }else{
    
      for (l in 1:L) {
        if(h_option=="n"){
          h_n <- c_h * (1 / n) ^ {1 / (2 * alpha + 1)}
        }else{
          h_n <- c_h * (1/ local_sizes[l]) ^ {1 / (2 * alpha + 1)}
        }
        time_start <- Sys.time()
        res <- GradientDescent(beta_0 = beta_init,
                                    x=x[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1]), ],
                                    y=y[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1])],
                                    h_n=h_n)
        time_end <- Sys.time()
        beta_hat <- res$beta
        beta_list[ ,l] <- beta_hat
        time_each_machine <- c(time_each_machine, time_end - time_start)
      }
  }  
    
  time_average_start <- Sys.time()
  beta_hat <- apply(local_sizes * t(beta_list) / n, 2, sum)
  time_average_end <- Sys.time()
  time_average <- time_average_end - time_average_start 
  time_total <- max(time_each_machine) + time_average
  
  ## Inference
  
  h_n_infer <- c_h * (1 / n) ^ {1 / (2 * alpha + 1)}
  h_delta <- c_h * (1 / n) ^ { delta / (2 * alpha + 1)}
  
  xbeta <- x %*% c(1, beta_hat)
  
  if(p == 1)
  {
    W_hat <- h_delta ^ {-alpha-1} * mean( c( y * dkerf(xbeta / h_delta)) * x[ , -1])
    V_hat <- (1/h_n_infer^2) * mean( c(- y * d2kerf(xbeta / h_n_infer)) * (x[ ,-1])^2)
    Va_hat <- (1 / h_n_infer) * mean((dkerf(xbeta / h_n_infer))^2 * (x[ ,-1])^2)
    
  }else{
    W_hat <- h_delta ^ {-alpha-1} * apply(c(y * dkerf(xbeta / h_delta)) * x[ , -1], 2, mean)
    V_hat <- matrix(0, nrow=p, ncol=p)
    yd2H <- c(y * d2kerf(xbeta / h_n_infer))
    Va_hat <- matrix(0, nrow = p, ncol = p)
    ydH2 <- (dkerf(xbeta / h_n_infer))^2
    for (i in 1:n) 
    {
      V_hat <- V_hat + (1/(n * h_n_infer^2)) * -yd2H[i] * (x[i, -1]) %*% t(x[i, -1])
      Va_hat <- Va_hat + (1 / (n * h_n_infer)) * ydH2[i] * (x[i, -1]) %*% t(x[i, -1])
    }
  }
  
  V_inv <- solve(V_hat)
  ViW <- V_inv %*% W_hat
  ViVaVi <- V_inv %*% Va_hat %*% V_inv
  # lambda_star <- lambda_init
  
  se <- 1 / (sqrt(n * h_n)) * sqrt(sum(ViVaVi))
  bias <- (h_n ^ (alpha)) * sum(V_inv %*% W_hat)
  
  return(list(beta_hat=beta_hat,
              beta_list=beta_list,
              bias=bias, 
              se=se,
              time=time_total,
              time_each_machine=time_each_machine))#,
              #lambda_star=lambda_star))
}

# p <- 10
# m <- 500
# logmn <- 1.8
# n <- floor(m ^ {logmn})
# L <- floor(n / m)
# errorTolerence <- 0.001
# maxIterations <- 10
# confidenceLevel <- 0.05
# z <- qnorm( 1 - confidenceLevel / 2, mean=0, sd=1)
# lambda <- 1 
# alpha <- 2 # kernel order
# beta_true <- rep(1/sqrt(p), p)
# local_sizes <- rep(m, L) # local size on each machine
# noiseType = "norm"
# 
# data = data_generate(n, p, beta=beta_true, noise_type=noiseType, noise_para=0.25)
# x <- data$x[1:sum(local_sizes), ]
# y <- data$y[1:sum(local_sizes)]
# # minIterations = 30
# v0 <- rep(1, p)
# kernelOrder = 2
# lambda_init=1
# delta=0.1
# beta_init="SMSE"
# lambda_opt = F
# minIterations = 5



## Our proposed multi-round method
mSMSE <- function(x,
                  y, 
                  local_sizes,
                  v0,
                  beta_init="SMSE",
                  kernelOrder=2,
                  c_h=1,
                  maxIterations=10,
                  errorTolerence=0.001,
                  minIterations=5,
                  lambda_opt=F,
                  lambda_h=1,
                  c_h_init_seq=c(1, 2, 3),
                  delta=0.1,
                  confidenceLevel=0.05,
                  verbose=F)
{
  n <- dim(x)[1]
  p <- dim(x)[2] - 1

  if(sum(local_sizes)!=n) stop("Sizes do not match!")
  
  L <- length(local_sizes)
  sizes_cumsum <- c(0, cumsum(local_sizes))
  
  
  ## build kernel function by order
  alpha <- kernelOrder

  kerf <- function(v) (abs(v) < 1) * (0.5 + 15/16 * (v - 2/3 * v^3 + 1/5 * v^5)) + 1*(v>=1)
  dkerf <- function(v) (abs(v) < 1) * (15/16 * (1 - v^2)^2)
  d2kerf <- function(v) (abs(v) < 1) * (15/8 * (1 - v^2)*(-2*v))
  
  z <- qnorm( (1 - confidenceLevel/2) , mean=0, sd=1)
  
  m <- local_sizes[1]
  Convergence <- FALSE
  initial_machine <- 1
  
  #Compute the gradient U and the Hessian V on a certain set
  
  UV_machine <- function(x, y, beta, h, h_delta)
  {
    time_start <- Sys.time()
    n <- nrow(x)
    p <- ncol(x) - 1
    xbeta <- x %*% c(1, beta)
    
    U <- 0
    V <- 0
    
    if(p == 1){
      
      U <- (1/h) * mean(c(- y * dkerf(xbeta / h)) * x[ , -1])
      V <- (1/h^2) * mean( c(- y * d2kerf(xbeta / h)) * (x[, -1])^2) 
      time_end <- Sys.time()
      
      U_hat <- h_delta ^ {-alpha-1} * mean( c( y * dkerf(xbeta / h_delta)) * x[ , -1])
      Vs_hat <- (1 / h_t) * mean((dkerf(xbeta / h))^2 * (x[ ,-1])^2)
      
    }else{
      
      U <- (1/h) * apply(c(- y * dkerf(xbeta / h)) * x[ , -1], 2, mean)
      V <- matrix(0, nrow=p, ncol=p)
      yd2H <- c(y * d2kerf(xbeta / h))
      
      for (i in 1:n) {
        V <- V + (1/(n * h^2)) * -yd2H[i] * (x[i, -1]) %*% t(x[i, -1])
      }
      
      # Compute the estimated matrices for inference. Not added in the total time.
      
      U_hat <- h_delta ^ {-alpha-1} * apply(c(y * dkerf(xbeta / h_delta)) * x[ , -1], 2, mean)
      
      Vs_hat <- matrix(0, nrow=p, ncol=p)
      dH2 <- (dkerf(xbeta / h))^2
      for (i in 1:n) {
        Vs_hat <- Vs_hat + (1 / (n * h)) * dH2[i] * (x[i, -1]) %*% t(x[i, -1])
      }
      time_end <- Sys.time()
    }
    return(list(U=U, V=V, U_hat=U_hat, Vs_hat=Vs_hat, time=time_end-time_start))
  }
  
  while(!Convergence){
    
  ## Compute the initial estimator beta_0
    if(all(beta_init=="SMSE")){
      beta_0 <- rep(0, p)
      initial_seq <- list()
      score_init_seq <- c()
      for(c_h_init in c_h_init_seq){
          initial <- Smoothed_MSE(x[(sizes_cumsum[initial_machine] + 1):(sizes_cumsum[initial_machine+1]), ], 
                              y[(sizes_cumsum[initial_machine] + 1):(sizes_cumsum[initial_machine+1])], 
                              beta_init=beta_0,
                              c_h=c_h_init,
                              #lambda_init=1,
                              inference=F)
          initial_seq <- append(initial_seq, list(initial$beta_hat))
          score_init <- if(p != 1) mean(y == sign(x[ ,1] + x[ ,-1] %*% initial$beta_hat)) else
            mean(y == sign(x[ ,1] + x[ ,-1] * initial$beta_hat))
          score_init_seq <- c(score_init_seq, score_init)
      }
    
      beta_0 <- initial_seq[[which.max(score_init_seq)]]
    
      if(verbose) print(paste("The initial value: ", toString(beta_0)))
    
      time_initial <- initial$time
    
    }#else{
      #    beta_0 <- beta_init
      #    time_initial <- 0
      #}
  
  #m = 500
  #mean(y[1:m] != sign(x[1:m, ] %*% c(1, beta_0)))
  #mean(y[1:m] != sign(x[1:m, ] %*% c(1, beta_true)))
  
  continueCondition <- TRUE
  iter <- 0
  beta_old <- beta_0
  beta_each_iteration <- matrix(NA, nrow=p, ncol=1)
  time_each_iteration <- c()
  interval_each_iteration <- matrix(NA, nrow=2, ncol=1)
  FindLambdaStar <- FALSE
  
  h_t_n <- c_h * (1 / n) ^ {1 / (2 * alpha + 1)}

  while (continueCondition) {
    
    iter <- iter + 1
    
    if(verbose) print(paste("Begin Iteration ", iter))
    
    # inference on beta_old
    
      
      h_t_m <- c_h * (1 / m) ^ {2 ^ iter / (3 * alpha)}
      h_t <- max(h_t_n, h_t_m) 
      h_delta <- c_h * (1 / n) ^ { delta / (2 * alpha + 1)}
      
      time_each_machine <- c()
      time_average <- 0
      U_t <- 0
      V_t <- 0
      U_hat <- 0
      Vs_hat <- 0
      
      for(l in 1:L)
      {
        UV_l <- UV_machine(x=x[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1]), ],
                           y=y[(sizes_cumsum[l] + 1): (sizes_cumsum[l+1])],
                           beta=beta_old,
                           h=h_t,
                           h_delta=h_delta)
        
        time_each_machine <- c(time_each_machine, UV_l$time) 
        time_average_start <- Sys.time()
        
        U_t <- U_t + (local_sizes[l] / n) * UV_l$U
        V_t <- V_t + (local_sizes[l] / n) * UV_l$V
        
        U_hat <- U_hat + (local_sizes[l] / n) * UV_l$U_hat
        Vs_hat <- Vs_hat + (local_sizes[l] / n) * UV_l$Vs_hat
        
        time_average_end <- Sys.time()
        time_average <- time_average  +  time_average_end - time_average_start
      }
       # Check the invertibility of V_t
      time_update_start <- Sys.time()
      
      if(p != 1){
        if(determinant(V_t)$modulus > -100){
          V_inverse <- solve(V_t)
        }else{
          V_inverse <- solve(V_t + 1e-5 * diag(p))
        }
      }else{
        V_inverse <- if(V_t!=0) 1 / V_t else 100
      }
      
      
      beta_new <- beta_old - V_inverse %*% U_t
      time_update_end <- Sys.time()
      time_update <- time_update_end - time_update_start
  
      if(verbose) print(paste("The estimated value in iteration", iter, "is ", toString(beta_new)))
      
      ViU <- V_inverse %*% U_hat
      ViVsVi <- V_inverse %*% Vs_hat %*% V_inverse
      time_iter <-  max(time_each_machine) + time_average + time_update
      beta_each_iteration <- cbind(beta_each_iteration, beta_new)
      time_each_iteration <- c(time_each_iteration, time_iter)
    

    # Construct prediction interval for sum(beta)
    interval_center <-  v0 %*% (beta_new - h_t^{alpha} * ViU)
    # interval_center <-  v0 %*% (beta_new - sqrt(lambda_h / (n * h_t)) * ViU)
    # interval_center <- sum(beta_new) - sqrt(lambda_h / (n * h_t)) * sum(ViU)
    # if(sum(ViVsVi) < 0) stop("The variance becomes negative!")
    interval_length <- sqrt((1 / (n * h_t)) * (v0 %*% ViVsVi %*% matrix(v0, ncol=1))) * z
    interval_left <- interval_center - interval_length
    interval_right <- interval_center + interval_length
    
    #if(iter == 1){
    #interval_center_0 <-  v0 %*% (beta_0 - h_t^{alpha} * ViU)
    #interval_length_0 <- interval_length
    #interval_left_0 <- interval_center_0 - interval_length
    #interval_right_0 <- interval_center_0 + interval_length
    #}
    
    interval_each_iteration <- cbind(interval_each_iteration, c(interval_left, interval_right))
  
    continueCondition <- ((sqrt(sum((beta_new - beta_old)^2)) > (errorTolerence * sqrt(p))) & (iter <= maxIterations)) | (iter <= minIterations)

    beta_old <- beta_new
    
    
    if(!continueCondition & !FindLambdaStar){
        lambda_h <- if(lambda_opt) (v0 %*% ViVsVi %*% matrix(v0, ncol=1)) / (2 * alpha * (v0 %*% ViU) ^ 2) 
                    else lambda_h
        h_t_n <- (lambda_h / n) ^ {1 / (2 * alpha + 1)}
        FindLambdaStar <- TRUE
        continueCondition <- TRUE
      }
    }
  
  
  time_each_iteration <- time_initial + cumsum(time_each_iteration)
  
  Convergence <- (iter <= maxIterations) | (initial_machine >= L)
  
  if(!Convergence) initial_machine <- initial_machine + 1
  # if(iter > maxIterations){
  #   beta_each_iteration[ ,1] = beta_0
  #   interval_each_iteration[ ,1] =  c(interval_left_0, interval_right_0)
  #   score <- c()
  #   for(t in 1:iter){
  #     score <- c(score, mean(y == sign(x[ ,1] + x[ ,-1] %*% beta_each_iteration[ ,t])))
  #   }
  #   best_id <- which.max(score)
  #   for(t in 1:(iter+1)){
  #     beta_each_iteration[ ,t] <- beta_each_iteration[ ,best_id]
  #     interval_each_iteration[ ,t] <- interval_each_iteration[ ,best_id]
  #   }
  }
  
  
  result <- list(beta_each_iteration = beta_each_iteration[ ,-1],
                 interval_each_iteration = interval_each_iteration[ , -1],
                 time_each_iteration = time_each_iteration,
                 IterationsUsed = iter,
                 FindLambdaStar = FindLambdaStar,
                 lambda_h = lambda_h,
                 beta_0 = beta_0,
                 n=n,
                 p=p,
                 m=m,
                 alpha=alpha)
  
  return(result)
  
}


Grad <- function(x,
                 y,
                 h,
                 beta){
  
  
  n <- dim(x)[1]
  xbeta <- x %*% beta
  if(p==2){
    grad <- (1/h) * mean(c(-y * dkerf(xbeta / h)) * x[ , -1]) 
  }else{
    grad <- (1/h) * apply(c(-y * dkerf(xbeta / h)) * x[ , -1], 2, mean)
  }
  return(grad)
  
}



GradientDescent <- function(x,
                            y,
                            h_n,
                            beta_0,
                            maxIterations=100,
                            errorTolerence=1e-4){
  p <- dim(x)[2] - 1
  
  kerf <- function(v) (abs(v) < 1) * (0.5 + 15/16 * (v - 2/3 * v^3 + 1/5 * v^5)) + 1*(v>=1)
  dkerf <- function(v) (abs(v) < 1) * (15/16 * (1 - v^2)^2)
  d2kerf <- function(v) (abs(v) < 1) * (15/8 * (1 - v^2)*(-2*v))
  
  beta_old <- beta_0
  xbeta_old <- x %*% c(1, beta_0)
  value_old <- mean( - y * kerf(xbeta_old / h_n))
  
  iter <- 0 
  
  ContinueCondition <- TRUE
  
  
  
  while(ContinueCondition) {
    
    iter <- iter + 1
    
    # print(iter)
    stepsize <- 1
    
    
    if(p == 1) 
    {
      grad <- (1/h_n) * mean( c( - y * dkerf(xbeta_old / h_n)) * x[ , -1]) 
      
    }else{
      
      grad <- (1/h_n) * apply(c(- y * dkerf(xbeta_old / h_n)) * x[ , -1], 2, mean)
      
    }
    
    direction <- -grad
    beta_new <- beta_old + stepsize * direction
    xbeta_new <- x %*% c(1, beta_new)
    value_new <- mean( - y * kerf(xbeta_new / h_n))
    while(value_new > value_old + 1e-6) 
    {
      stepsize <- stepsize * 0.5
      # print(stepsize)
      beta_new <- beta_old + stepsize * direction
      xbeta_new <- x %*% c(1, beta_new)
      value_new <- mean( - y * kerf(xbeta_new / h_n))
    }
    
    ContinueCondition <- ((sqrt(sum((beta_new - beta_old)^2))) > sqrt(p) * errorTolerence) & (iter < maxIterations) 
    
    beta_old <- beta_new
    xbeta_old <- xbeta_new
    value_old <- value_new
  }
  
  res <- list(beta=beta_new,
              value=value_new)
  
  return(res)
}






SubOpti <- function(lambda, beta, grad){
  subopti <- rep(0, length(beta))
  idx_0 <- (beta == 0) & (abs(grad) >= lambda)
  if(any(idx_0)){ 
    temp <- min(abs(grad[idx_0] + lambda), abs(grad[idx_0] - lambda))
    subopti[idx_0] <- temp
  }
  idx_1 <- (beta != 0)
  subopti[idx_1] <- abs(grad[idx_1] + lambda * sign(beta[idx_1]))
  return(max(subopti))
}



SoftTreshold <- function(lambda, eta, beta, grad){
  beta_new <- rep(0, length(beta))
  beta_bar <- beta - eta*grad
  idx <- abs(beta_bar) > lambda*eta
  beta_new[idx] = sign(beta_bar[idx])*(abs(beta_bar)[idx] - lambda*eta)
  return(beta_new)
}


ProxGrad <- function(x,
                     y,
                     h,
                     beta,
                     epsilon,
                     lambda,
                     eta,
                     maxiter = 10000,
                     stage = 'NA'){
  #initialization
  grad <- Grad(x, y, h, beta)
  subopti <- SubOpti(lambda, beta[-1], grad)
  idx <- 0
  while(subopti > epsilon & idx < maxiter){
    idx <- idx + 1
    beta[-1] <- SoftTreshold(lambda, eta, beta[-1], grad)
    grad <- Grad(x, y, h, beta)
    subopti <- SubOpti(lambda, beta[-1], grad)
  }
  
  
  #if(idx == maxiter)
  # print(paste('Maximum iteration reached with sub-optimality ', subopti, " at stage ", stage, sep=""))
  
  return(beta)
}


PathFollowing <- function(x, y, h, lambda_tgt,
                          epsilon_final = 0.0001,
                          beta_0,
                          stages = 10,
                          nu = 0.25,
                          eta = 0.5,
                          maxiter_it = 1000,
                          maxiter_final = 2000){
  
  
  beta <- beta_0 
  lambda <- max(abs(Grad(x, y, h, beta)))
  
  #print(lambda)
  
  if(stages == 0){
    
    beta <- ProxGrad(x,
                     y,
                     h,
                     beta,
                     epsilon_final,
                     lambda_tgt,
                     eta,
                     maxiter = maxiter_final, 
                     stage=toString(stage))
    
    value <- mean( - y * kerf(x %*% beta / h))
    
    return(list(beta=beta, value=value))
    
  }
  
  phi <- (lambda_tgt/lambda)^(1.0/stages)
  
  #print(phi)
  
  
  
  for(stage in 1:stages){
    
    #print(stage)
    
    lambda <- lambda * phi
    
    #print(lambda)
    epsilon <- nu * lambda
    
    beta <- ProxGrad(x,
                     y,
                     h,
                     beta,
                     epsilon,
                     lambda,
                     eta,
                     maxiter = maxiter_it, 
                     stage=toString(stage))
  }
  
  beta <- ProxGrad(x, y, h, beta, epsilon_final, lambda_tgt, eta, maxiter=maxiter_final, stage = 'Final')
  
  value <- mean( - y * kerf(x %*% beta / h))
  
  return(list(beta=beta, value=value))
  
}




