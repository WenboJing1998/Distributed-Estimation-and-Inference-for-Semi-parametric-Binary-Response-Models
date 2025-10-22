import numpy as np
import cvxpy as cp
import os
import itertools
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
from scipy.stats import norm
from pathfollowing import *
from scipy.linalg import toeplitz


def data_generate(n, p, beta, noise_type, noise_para):
    rho = 0.5
    cov_mat = toeplitz(rho ** (np.arange(0, p + 1, 1)), rho ** (np.arange(0, p + 1, 1)))
    Covariates = np.array(np.random.multivariate_normal([0]*(p + 1), cov_mat, n)).reshape(n, (p + 1))
    X = Covariates[:, 0].reshape(n)
    Z = Covariates[:, 1:]
    XZbeta = X + Z @ beta 
    epsilon = np.random.normal(0, noise_para, n) if noise_type == "homo" else \
               noise_para * ((1 + 2 * (XZbeta) ** 2 )) * np.random.normal(0, noise_para, n) 
    y = (2 * ((XZbeta + epsilon) >= 0) - 1).reshape(n)
    return {'X': X,
            'Y': y,
            'Z': Z,
            'true_coeff' : beta}


def mSMSE_hd_single_iteration_cv(X, 
                                 y, 
                                 Z, 
                                 local_sizes,
                                 beta_old,
                                 c_h_seq,
                                 lambda_seq,
                                 K=5,
                                 kernelOrder=6
                                 ):
       
    n, p = Z.shape[0], Z.shape[1]

    if(sum(local_sizes)!=n):
        Warning("Sizes do not match!")
  
    L = len(local_sizes)
    # sizes_cumsum = np.insert(np.cumsum(local_sizes), 0, 0)

    ## build kernel function by order
    alpha = kernelOrder

    if(kernelOrder==6):
        dkerf = lambda v: (abs(v) < 1) * ((15/16) * (1-v**2)**2) * ((315/128)*(1-(22/3)*(v**2)+(143/15)*(v**4)))
        d2kerf = lambda v: (abs(v) < 1) * ((315/256) * v * (-35 + 189 * v ** 2 - 297 * v ** 4 + 143 * v ** 6))
    elif(kernelOrder==2):
        dkerf = lambda v: (abs(v) < 1) * (105/64 * (1 - 5 * (v ** 2) + 7 * (v ** 4) - 3 * (v ** 6)))
        d2kerf = lambda v: (abs(v) < 1) * (105/64 * (- 10 * v + 28 * (v ** 3) - 18 * (v ** 5)))
  
    num_h = len(c_h_seq)
    num_lambda = len(lambda_seq)
    cv_errors = np.zeros((num_h, num_lambda))
  # time_each_iteration = []
    for i, j in itertools.product(range(num_h), range(num_lambda)):

        c_h = c_h_seq[i]
        lambda_t = lambda_seq[j]

        h_t = c_h * (np.log(p) / n) ** (1 / (2 * alpha + 1))
        
        # print(f"Begin CV for bandwidth {h_t} and penalty parameter {lambda_t}:")

        cv_error = []

        kfolds = list(KFold(5).split(X, y, Z))

        for k in range(K):

            train_id = kfolds[k][0]
            test_id = kfolds[k][1]

            len_1 = int((1 - (1 / K)) * local_sizes[0])
            Z_1 = Z[train_id[:len_1], :]
            XZbeta_1 = X[train_id[:len_1]] + Z_1 @ beta_old
            yd2H_1 = -y[train_id[:len_1]] * d2kerf(XZbeta_1 / h_t)
            V_1 = (1 / (len_1 * (h_t ** 2))) * (yd2H_1 * np.transpose(Z_1)) @ Z_1

            
            XZbeta = X[train_id] + Z[train_id, :] @ beta_old
            U_t = np.mean((1 / h_t) * np.transpose((- y[train_id] * \
                                dkerf(XZbeta / h_t)) * np.transpose(Z[train_id])), axis=0)
            
            G_t = (V_1 @ beta_old - U_t).reshape(p, 1)


            beta_variable = cp.Variable((p, 1))
            u_variable = cp.Variable((p, 1))
            objective = cp.Minimize(cp.sum(u_variable))

            # lambda_t = max(abs(V_1 @ beta.reshape(p, 1) - G_t))

            constraints = [beta_variable >= -u_variable,
                            beta_variable <= u_variable,
                            V_1 @ beta_variable - G_t <= lambda_t,
                            V_1 @ beta_variable - G_t >= -lambda_t]
            
            prob = cp.Problem(objective, constraints)

            prob.solve(verbose=False)

            if(prob.status == "optimal"):
                beta_k = np.copy(beta_variable.value.reshape(p))
                cv_error.append(np.mean(y[test_id] != np.sign(X[test_id] + Z[test_id] @ beta_k)))
            else:
                cv_error.append(100)

        cv_errors[i, j] = np.mean(np.array(cv_error))

    opt_id = np.unravel_index(np.argmin(cv_errors), shape=(num_h, num_lambda))
    c_h_opt = c_h_seq[opt_id[0]]
    lambda_opt = lambda_seq[opt_id[1]]
    cv_error_min = cv_errors[opt_id[0], opt_id[1]]
    return c_h_opt, lambda_opt, cv_error_min, cv_errors


def mSMSE_hd(X, y, Z, 
             local_sizes,
             beta,
             cv=True,
             c_h_seq=[5, 10],
             lambda_seq=[0.01, 0.05],
             beta_0="pf",
             kernelOrder=6,
             maxIterations=10,
             errorTolerence=5e-4,
             #minIterations=5,
             c_init_h=1,
             c_init_lambda=0.5):
       
  n, p = Z.shape[0], Z.shape[1]

  if(sum(local_sizes)!=n):
    Warning("Sizes do not match!")
  
  L = len(local_sizes)
  # sizes_cumsum = np.insert(np.cumsum(local_sizes), 0, 0)

  ## build kernel function by order
  alpha = kernelOrder

  if(kernelOrder==6):
    dkerf = lambda v: (abs(v) < 1) * ((15/16) * (1-v**2)**2) * ((315/128)*(1-(22/3)*(v**2)+(143/15)*(v**4)))
    d2kerf = lambda v: (abs(v) < 1) * ((315/256) * v * (-35 + 189 * v ** 2 - 297 * v ** 4 + 143 * v ** 6))
  elif(kernelOrder==2):
    dkerf = lambda v: (abs(v) < 1) * (105/64 * (1 - 5 * (v ** 2) + 7 * (v ** 4) - 3 * (v ** 6)))
    d2kerf = lambda v: (abs(v) < 1) * (105/64 * (- 10 * v + 28 * (v ** 3) - 18 * (v ** 5)))
  
  
  ## Compute the initial estimator beta_0
  
  if(beta_0=="pf"): 

    h_1 = c_init_h * (s * np.log(p) / local_sizes[0]) ** (1 / 5)
    lambda_1 = c_init_lambda * np.sqrt(np.log(p)/(local_sizes[0] * h_1))
    
    data_first = {"X": X[:local_sizes[0]], "Y": y[:local_sizes[0]], "Z": Z[:local_sizes[0], :], 'true_coeff' : beta}
    beta_0 = -PathFollowing(data_first, h_1, norm.pdf,  lambda_1) 
    
    print(f"The initial error: {np.sum((beta_0 - beta)**2)}")
    # time_initial = initial$time
  #else:
    # print(f"The initial error: {np.sum((beta_0 - beta)**2)}")
    # time_initial <- 0

  continueCondition = True
  iter = 0
  beta_old = np.copy(beta_0)
  # cv_error_min_old = np.inf

  # beta_each_iteration = np.zeros(p, 1)
  # time_each_iteration = []
  error = []
  # beta_diff = []
  # cv_error = []
  beta_hats = []
  training_error = []

  while (continueCondition):

    print(f"Begin Iteration {iter}:")

    if(cv):
      c_h, lambda_t, cv_error_min_new, _ = mSMSE_hd_single_iteration_cv(X, y, Z, 
                                                                        local_sizes, 
                                                                        beta_old, 
                                                                        c_h_seq, 
                                                                        lambda_seq, 
                                                                        K=5, 
                                                                        kernelOrder=kernelOrder)
      # cv_error.append(cv_error_min_new)
      # if(cv_error_min_new > cv_error_min_old + 0.01):
      #   break
      # else:
      #   cv_error_min_old = np.copy(cv_error_min_new)
    else:
       c_h = c_h_seq[0]
       lambda_t = lambda_seq[0]

    h_t = c_h * (np.log(p) / n) ** (1 / (2 * alpha + 1))

    print(f"Iteration:{iter}, badwidth:{c_h}, lambda:{lambda_t}")

    XZbeta_1 = X[:local_sizes[0]] + Z[:local_sizes[0], :] @ beta_old
    yd2H_1 = -y[:local_sizes[0]] * d2kerf(XZbeta_1 / h_t)
    V_1 = (1 / (local_sizes[0] * (h_t ** 2))) * (yd2H_1 * np.transpose(Z[:local_sizes[0], :])) @ Z[:local_sizes[0], :]

    XZbeta = X + Z @ beta_old
    U_t = np.mean((1 / h_t) * np.transpose((- y * dkerf(XZbeta / h_t)) * np.transpose(Z)), axis=0)
    # for l in range(L):
    #    XZbeta = X[(sizes_cumsum[l]):(sizes_cumsum[l+1])] + Z[(sizes_cumsum[l]):(sizes_cumsum[l+1]), :] @ beta_old
    #    U_l = np.mean((1 / h_t) * np.transpose((- y[(sizes_cumsum[l]):(sizes_cumsum[l+1])] * \
    #                            dkerf(XZbeta / h_t)) * np.transpose(Z[(sizes_cumsum[l]):(sizes_cumsum[l+1])])), axis=0)
    #   U_t += (local_sizes[l] / n) * U_l

    G_t = (V_1 @ beta_old - U_t).reshape(p, 1)

    beta_variable = cp.Variable((p, 1))
    u_variable = cp.Variable((p, 1))
    objective = cp.Minimize(cp.sum(u_variable))
    lambda_t = max(abs(V_1 @ beta.reshape(p, 1) - G_t))
    constraints = [beta_variable >= -u_variable,
                    beta_variable <= u_variable,
                    V_1 @ beta_variable - G_t <= lambda_t,
                    V_1 @ beta_variable - G_t >= -lambda_t]
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose=False)
    beta_new = np.copy(beta_variable.value.reshape(p))

    training_error.append(np.mean(y != np.sign(X + Z @ beta_new)))
    #beta_diff.append(np.sum((beta_new - beta_old) ** 2))

    #if(iter > 4):
    #  if(beta_diff[-1] > beta_diff[-2]):
    #    break

    continueCondition = np.sum((beta_new - beta_old) ** 2) > errorTolerence and \
                        iter < maxIterations
    
    error.append(np.sum((beta_new - beta) ** 2)) 
    beta_hats.append(beta_new)
    #print(f"Error: {error[-1]}")
    #break
    #beta_old = np.copy(beta_new)
    #error.append(np.sum((beta_new - beta) ** 2))
    #print(f"Error: {error[-1]} CV Error:{cv_error[-1]} train_error={training_error[-1]}")
    print(f"Error: {error[-1]} training_error={training_error[-1]}")
    iter += 1
    beta_old = np.copy(beta_new)


  if(iter > maxIterations):
     best_id = int(np.argmin(np.array(training_error)))
     beta_hat = beta_hats[best_id]
     final_error = np.sum((beta_hat - beta) ** 2)
     print(f"Final Error: {final_error}")
     #result = {"error":error, 
     #          "final_error":final_error,
     #         "beta_hat":beta_hat,
     #         "beta_hats":beta_hats, 
     #         "beta_0": beta_0}
              #"beta_diff": beta_diff#,
              #"cv_error":cv_error}
  else:
     beta_hat = beta_new
     final_error = np.sum((beta_hat - beta) ** 2)
     print(f"Final Error: {np.sum((beta_hat - beta) ** 2)}")

  result = {"error":error, 
            "final_error":final_error,
            "beta_hat":beta_hat,
            "beta_hats":beta_hats, 
            "beta_0": beta_0}
            #"beta_diff": beta_d
  return(result)


def pf_cv(X, y, Z, h_seq, c_lambda_seq, K=5):
       
    n, p = Z.shape[0], Z.shape[1]
    num_h = len(h_seq)
    num_lambda = len(c_lambda_seq)
    cv_errors = np.zeros((num_h, num_lambda))
  # time_each_iteration = []
    for i, j in itertools.product(range(num_h), range(num_lambda)):

        h_pf = h_seq[i]
        c_lambda = c_lambda_seq[j]

        lambda_pf = c_lambda * np.sqrt(np.log(p) / (n * h_pf))
        
        # print(f"Begin CV for bandwidth {h_t} and penalty parameter {lambda_t}:")

        cv_error = []

        kfolds = list(KFold(5).split(X, y, Z))

        for k in range(K):

            train_id = kfolds[k][0]
            test_id = kfolds[k][1]

            data_train = {"X": X[train_id], "Y": y[train_id], "Z": Z[train_id, :]}
            beta_pf_k = -PathFollowing(data_train, h_pf, norm.pdf, lambda_pf)

            cv_error.append(np.mean(y[test_id] != np.sign(X[test_id] + Z[test_id] @ beta_pf_k)))

        cv_errors[i, j] = np.mean(np.array(cv_error))

    opt_id = np.unravel_index(np.argmin(cv_errors), shape=(num_h, num_lambda))
    h_opt = h_seq[opt_id[0]]
    c_lambda_opt = c_lambda_seq[opt_id[1]]
    cv_error_min = cv_errors[opt_id[0], opt_id[1]]

    return h_opt, c_lambda_opt, cv_error_min, cv_errors

L_seq = np.arange(6, 26, 2)

taskId = os.environ.get("SLURM_ARRAY_TASK_ID")
# taskId <- 2999
taskId = int(taskId) - 1
# np.savetxt(f"taskID.txt", np.array([taskId]))

m = 400
L = int(L_seq[taskId])
n = m * L
s = 10
p = 500
beta = np.array([-1 / np.sqrt(s)] * (s) + [0] * (p-s))
noise_type = "homo"
noise_para = 0.5
local_sizes = np.array([m] * L)
dkerf = lambda v: (abs(v) < 1) * ((15 / 16) * (1 - v ** 2) ** 2) * ((315 / 128) * (1 - (22 / 3) * (v ** 2) + (143 / 15) * (v ** 4)))
d2kerf = lambda v: (abs(v) < 1) * ((315 / 256) * v * (-35 + 189 * v ** 2 - 297 * v ** 4 + 143 * v ** 6))

num_rep = 500
errors = []

np.random.seed(taskId)

for rep in range(num_rep):

    print(f"Begin reptition:{rep}")
    # np.savetxt(f"rep.txt", np.array([rep]))

    data = data_generate(n, p, beta, noise_type, noise_para)
    X, y, Z = data["X"], data["Y"], data["Z"]
    n, p = Z.shape[0], Z.shape[1]

    # 1. Ave-DC

    sizes_cumsum = np.insert(np.cumsum(local_sizes), 0, 0)
    beta_each = np.zeros((L, p)) 
    for l in range(L):
        data_ell = {"X": X[(sizes_cumsum[l]):(sizes_cumsum[l+1])], "Y": y[(sizes_cumsum[l]):(sizes_cumsum[l+1])], "Z": Z[(sizes_cumsum[l]):(sizes_cumsum[l+1]), :], 'true_coeff' : beta}
        h_l, c_lambda_l, _, _ = pf_cv(X[(sizes_cumsum[l]):(sizes_cumsum[l+1])], \
                                      y[(sizes_cumsum[l]):(sizes_cumsum[l+1])], \
                                      Z[(sizes_cumsum[l]):(sizes_cumsum[l+1]), :], 
                                      h_seq=[0.5, 1.0], 
                                      c_lambda_seq=[0.5, 1], 
                                      K=5)
        # print(f"{h_l}, {c_lambda_l}")
        beta_l = -PathFollowing(data_ell, h_l, norm.pdf, c_lambda_l * np.sqrt(np.log(p) / (local_sizes[l] * h_l)))
        beta_each[l, :] = np.copy(beta_l)

    beta_DC = np.mean(beta_each, axis=0)
    error_DC = np.sum((beta_DC - beta)**2)

    # 2. mSMSE-hd

    res_mSMSE = mSMSE_hd(X, y, Z, local_sizes, beta, cv=False,
                         c_h_seq=[7], lambda_seq=[0.01], 
                         maxIterations=10, errorTolerence=1e-3,
                         c_init_h=1, c_init_lambda=0.1)
    error_mSMSE = res_mSMSE["final_error"]

    # 3. pooled-pf
    data_pf = {"X": X, "Y": y, "Z": Z}
    h_opt, c_lambda_opt, _, _ = pf_cv(X, y, Z, h_seq=[0.5, 0.8], c_lambda_seq=[0.5, 0.8], K=5)
    beta_pf = -PathFollowing(data_pf, h_opt, norm.pdf, c_lambda_opt * np.sqrt(np.log(p) / (n * h_opt)))
    error_pooled = np.sum((beta_pf - beta)**2)

    print(f"Errors: DC {error_DC}, mSMSE {error_mSMSE}, pooled-pf {error_pooled}")

    errors.append([error_DC, error_mSMSE, error_pooled])

np.savetxt(f"results/L2err_n={n}_m={m}_p={p}.txt", np.array(errors))

# errors_load = np.loadtxt(f"L2err_n={n}_s={s}_p={p}.txt")
# np.mean(np.array(errors), axis=0)