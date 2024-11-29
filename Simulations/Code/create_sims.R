#====================================
# Run simulation scenarios
# Programmer: SM Wu
# Date updated: 2024/11/25
#====================================

rm(list=ls())

### Workflow
library(tidyverse)
library(ggpubr)
library(gtools)

wd <- "~/Documents/Github/WRPC/"
wd <- "/n/netscratch/stephenson_lab/Lab/stephwu18/WRPC/"
code_dir <- "Simulations/Code/"
data_dir <- "Simulations/Data/"
res_dir <- "Simulations/Results/"

# Data simulation code
source(paste0(wd, "Simulations/Code/simulate_data_wrpc.R"))

# # Read in data
# load(paste0(wd, data_dir, "w1_sim_pop.RData"))
# load(paste0(wd, data_dir, "w1_sim_samp.RData"))
# data_vars <- sim_samp

#================= Simulate data ===============================================
#================= Scen 1 & 2: 2 strata, X ~ C =================================

# Create scenario results folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_1/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}
dir_path <- paste0(wd, data_dir, "scen_2/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### Generate population
# Population size and strata dimensions
N <- 40000; H <- 2; S <- 2; J <- 30; R <- 5 
N_h <- c(15000, 25000)
N_s <- c(30000, 10000)

# Generate C ~ S
K <- 3  
formula_c <- "~ s_all"
V_unique <- data.frame(s_all = as.factor(1:S))
# Corresponds to an overall true_pi ~= (0.199, 0.552, 0.249)
# Matrix of global class assignment probabilities for each level of s_i
pi_global_mat <- matrix(c(0.1, 0.2, 0.7,   # global class probs for s_i=1
                          0.8, 0.1, 0.1),  # global class probs for s_i=2
                        byrow = TRUE, nrow = S, ncol = K)
beta_mat_c <- get_betas_c(pi_mat = pi_global_mat, formula_c = formula_c, 
                          V_unique = V_unique)

# Generate H 
h_all <- c(rep(1, times = N_h[1]), rep(2, times = N_h[2]))

# Global-local allocation probability of global assignment
nu <- matrix(c(0.7, 0.3), byrow = FALSE, nrow = H, ncol = J)

# Generate X ~ C
formula_x_global <- "~ c_all"
formula_x_local <- "~ h_all"
modal_theta_prob <- 0.85
V_unique <- expand.grid(c_all = as.factor(1:K), h_all = as.factor(1:H))
# Matrix of global profiles
profiles_global <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
                                               rep(3, times = 0.5 * J)),
                                        C2 = c(rep(4, times = 0.2 * J), 
                                               rep(2, times = 0.8 * J)),
                                        C3 = c(rep(5, times = 0.3 * J), 
                                               rep(4, times = 0.4 * J),
                                               rep(1, times = 0.3 * J))))
# List of beta matrices for each item j
beta_list_x_global <- get_betas_x(profiles = profiles_global, R = R, 
                                  modal_theta_prob = modal_theta_prob, 
                                  formula_x = formula_x_global, 
                                  V_unique = V_unique)
# List of matrices of local profiles for each subpopulation
profiles_local <- as.matrix(data.frame(H1 = c(rep(1, times = J/2),
                                              rep(2, times = J/2)),
                                       H2 = c(rep(3, times = J/3),
                                              rep(4, times = J/3),
                                              rep(5, times = J/3))))
# List over subpopulation of list of matrices of beta matrices for each item j
beta_list_x_local <- get_betas_x(profiles = profiles_local, R = R, 
                                 modal_theta_prob = modal_theta_prob, 
                                 formula_x = formula_x_local, 
                                 V_unique = V_unique)

# Simulate population
pop_seed <- 1  # Set seed
sim_pop <- simulate_pop_wrpc_v2(N = N, H = H, S = S, J = J, R = R, K = K, 
                                N_h = N_h, N_s = N_s,
                                modal_theta_prob = modal_theta_prob,
                                formula_c = formula_c, 
                                formula_x_global = formula_x_global, 
                                formula_x_local = formula_x_local, 
                                beta_mat_c = beta_mat_c, 
                                beta_list_x_global = beta_list_x_global,
                                beta_list_x_local = beta_list_x_local,
                                pop_seed = pop_seed, save_res = TRUE, 
                                save_path = paste0(wd, data_dir, "scen_1/"))
# Duplicate for scenario 2
save(sim_pop, file = paste0(wd, data_dir, "scen_2/sim_pop.RData"))

### Stratified sampling only affecting pi
num_samps <- 100
for (i in 1:num_samps) {
  # Scenario 1: Sample size 2000
  sim_samp <- simulate_samp_wrpc(sim_pop = sim_pop, samp_size = 2000, strat = TRUE, 
                                 strat_dist = c(0.3, 0.7), clust = FALSE, 
                                 samp_seed = i, save_res = TRUE, 
                                 save_path = paste0(wd, data_dir, "scen_1/", i))
  
  # Scenario 2: Sample size 400
  sim_samp <- simulate_samp_wrpc(sim_pop = sim_pop, samp_size = 400, strat = TRUE, 
                                 strat_dist = c(0.3, 0.7), clust = FALSE, 
                                 samp_seed = i, save_res = TRUE, 
                                 save_path = paste0(wd, data_dir, "scen_2/", i))
}


#================= Scen 11 & 12: 2 strata, X ~ C + S ===========================

### Generate population
# Population size and strata dimensions
N <- 40000; H <- 2; S <- 2; J <- 30; R <- 5 
N_h <- c(15000, 25000)
N_s <- c(20000, 20000)

# Generate C ~ S
K <- 3  
formula_c <- "~ s_all"
V_unique <- data.frame(s_all = as.factor(1:S))
# Corresponds to an overall true_pi ~= (0.199, 0.552, 0.249)
# Matrix of global class assignment probabilities for each level of s_i
pi_global_mat <- matrix(c(0.3, 0.5, 0.2,   # global class probs for s_i=1
                          0.1, 0.3, 0.6),  # global class probs for s_i=2
                        byrow = TRUE, nrow = S, ncol = K)
beta_mat_c <- get_betas_c(pi_mat = pi_global_mat, formula_c = formula_c, 
                          V_unique = V_unique)

# Generate H 
h_all <- c(rep(1, times = N_h[1]), rep(2, times = N_h[2]))

# Generate S
s_all <- c(rep(1, times = N_s[1]), rep(2, times = N_s[2]))

# Global-local allocation probability of global assignment
nu <- matrix(c(0.7, 0.3), byrow = FALSE, nrow = H, ncol = J)

# Generate X ~ C + S for global patterns
formula_x_global <- "~ c_all"
modal_theta_prob <- 0.85
V_unique <- expand.grid(c_all = as.factor(1:K), h_all = as.factor(1:H))
# Matrix of global profiles
profiles_global <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
                                               rep(3, times = 0.5 * J)),
                                        C2 = c(rep(4, times = 0.2 * J), 
                                               rep(2, times = 0.8 * J)),
                                        C3 = c(rep(5, times = 0.3 * J), 
                                               rep(4, times = 0.4 * J),
                                               rep(1, times = 0.3 * J))))
# List of beta matrices for each item j
beta_list_x_temp <- get_betas_x(profiles = profiles_global, R = R, 
                                  modal_theta_prob = modal_theta_prob, 
                                  formula_x = formula_x_global, 
                                  V_unique = V_unique)
# Add in coefficients for strata s_all, updating formula_x and beta_list_x
formula_x_global <- "~ c_all + s_all + c_all:s_all"
beta_list_x_1_28 <- lapply(1:(J-2), function(j)
  cbind(beta_list_x_temp[[j]], s_all2 = rep(0, 5),
        `c_all2:s_all2` = rep(0, 5), `c_all3:s_all2` = rep(0, 5)))
# Items 29-30 are affected as follows: s_all=2 associated with r=4,5 for k=1,2,3
beta_list_x_29_30 <- lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], s_all2 = c(0, -0.5, -0.5, 0.5, 0.5), 
        `c_all2:s_all2` = c(0, -0.5, -0.5, 0.5, 0.5), 
        `c_all3:s_all2` = c(0, -0.5, -0.5, 0.5, 0.5)))
beta_list_x_global <- c(beta_list_x_1_28, beta_list_x_29_30)

# Generate X ~ H + S for local patterns
formula_x_local <- "~ h_all"
# List of matrices of local profiles for each subpopulation
profiles_local <- as.matrix(data.frame(H1 = c(rep(1, times = J/2),
                                              rep(2, times = J/2)),
                                       H2 = c(rep(3, times = J/3),
                                              rep(4, times = J/3),
                                              rep(5, times = J/3))))
# List over subpopulation of list of matrices of beta matrices for each item j
beta_list_x_temp <- get_betas_x(profiles = profiles_local, R = R, 
                                 modal_theta_prob = modal_theta_prob, 
                                 formula_x = formula_x_local, 
                                 V_unique = V_unique)
# Add in coefficients for strata s_all, updating formula_x and beta_list_x
formula_x_local <- "~ h_all + s_all + h_all:s_all"
beta_list_x_1_28 <- lapply(1:(J-2), function(j)
  cbind(beta_list_x_temp[[j]], s_all2 = rep(0, 5), `h_all2:s_all2` = rep(0, 5)))
# Items 29-30 are affected as follows: s_all=2 associated with r=3,4 for h=1,2
beta_list_x_29_30 <- lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], s_all2 = c(0, -0.5, 0.5, 0.5, -0.5), 
        `h_all2:s_all2` = c(0, -0.5, 0.5, 0.5, -0.5)))
beta_list_x_local <- c(beta_list_x_1_28, beta_list_x_29_30)

# Simulate population
pop_seed <- 1  # Set seed
sim_pop <- simulate_pop_wrpc_v2(N = N, H = H, S = S, J = J, R = R, K = K, 
                                N_h = N_h, N_s = N_s,
                                modal_theta_prob = modal_theta_prob,
                                formula_c = formula_c, 
                                formula_x_global = formula_x_global, 
                                formula_x_local = formula_x_local, 
                                beta_mat_c = beta_mat_c, 
                                beta_list_x_global = beta_list_x_global,
                                beta_list_x_local = beta_list_x_local,
                                pop_seed = pop_seed, save_res = TRUE, 
                                save_path = paste0(wd, data_dir, "scen_11/"))
# Duplicate for scenario 12
save(sim_pop, file = paste0(wd, data_dir, "scen_12/sim_pop.RData"))

### Stratified sampling only affecting pi
num_samps <- 100
for (i in 1:num_samps) {
  # Scenario 11: Sample size 2000
  sim_samp <- simulate_samp_wrpc(sim_pop = sim_pop, samp_size = 2000, strat = TRUE, 
                                 strat_dist = c(0.3, 0.7), clust = FALSE, 
                                 samp_seed = i, save_res = TRUE, 
                                 save_path = paste0(wd, data_dir, "scen_11/", i))
  
  # Scenario 12: Sample size 400
  sim_samp <- simulate_samp_wrpc(sim_pop = sim_pop, samp_size = 400, strat = TRUE, 
                                 strat_dist = c(0.3, 0.7), clust = FALSE, 
                                 samp_seed = i, save_res = TRUE, 
                                 save_path = paste0(wd, data_dir, "scen_12/", i))
}


#================= Scen 21 & 22: 2 strata, X ~ C + S (V2) ======================
# Create scenario results folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_21/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}
dir_path <- paste0(wd, data_dir, "scen_22/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### Generate population
# Population size and strata dimensions
N <- 40000; H <- 2; S <- 2; J <- 30; R <- 5 
N_h <- c(15000, 25000)
N_s <- c(30000, 10000)

# Generate C ~ S
K <- 3  
formula_c <- "~ s_all"
V_unique <- data.frame(s_all = as.factor(1:S))
# Corresponds to an overall true_pi ~= (0.199, 0.552, 0.249)
# Matrix of global class assignment probabilities for each level of s_i
pi_global_mat <- matrix(c(0.1, 0.2, 0.7,   # global class probs for s_i=1
                          0.8, 0.1, 0.1),  # global class probs for s_i=2
                        byrow = TRUE, nrow = S, ncol = K)
beta_mat_c <- get_betas_c(pi_mat = pi_global_mat, formula_c = formula_c, 
                          V_unique = V_unique)

# Generate H 
h_all <- c(rep(1, times = N_h[1]), rep(2, times = N_h[2]))

# Generate S
s_all <- rep(c(rep(1, times = (N_s[1]/4)), rep(2, times = (N_s[2]/4))), 4)

# Global-local allocation probability of global assignment
nu <- matrix(c(0.7, 0.3), byrow = FALSE, nrow = H, ncol = J)

# Generate X ~ C + S for global patterns
formula_x_global <- "~ c_all"
modal_theta_prob <- 0.85
V_unique <- expand.grid(c_all = as.factor(1:K), h_all = as.factor(1:H))
# Matrix of global profiles
profiles_global <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
                                               rep(3, times = 0.5 * J)),
                                        C2 = c(rep(4, times = 0.2 * J), 
                                               rep(2, times = 0.8 * J)),
                                        C3 = c(rep(5, times = 0.3 * J), 
                                               rep(4, times = 0.4 * J),
                                               rep(1, times = 0.3 * J))))
# List of beta matrices for each item j
beta_list_x_temp <- get_betas_x(profiles = profiles_global, R = R, 
                                modal_theta_prob = modal_theta_prob, 
                                formula_x = formula_x_global, 
                                V_unique = V_unique)
# Add in coefficients for strata s_all, updating formula_x and beta_list_x
formula_x_global <- "~ c_all + s_all + c_all:s_all"
beta_list_x_1_27 <- lapply(1:(J-3), function(j)
  cbind(beta_list_x_temp[[j]], s_all2 = rep(0, 5),
        `c_all2:s_all2` = rep(0, 5), `c_all3:s_all2` = rep(0, 5)))
# Items 28-30 are affected as follows: s_all=2 associated with r=4,5 for k=1,
# with r=3,4,5 for k=2, and with r=3,4,5 for k=3
beta_list_x_28_30 <- lapply((J-2):J, function(j) 
  cbind(beta_list_x_temp[[j]], s_all2 = c(0, 0, -0.5, 1, 1), 
        `c_all2:s_all2` = c(0, -2, 1, 0.5, 0.5), 
        `c_all3:s_all2` = c(0, 0, 2, 0.5, 1)))
beta_list_x_global <- c(beta_list_x_1_27, beta_list_x_28_30)

# Generate X ~ H + S for local patterns
formula_x_local <- "~ h_all"
# List of matrices of local profiles for each subpopulation
profiles_local <- as.matrix(data.frame(H1 = c(rep(1, times = J/2),
                                              rep(2, times = J/2)),
                                       H2 = c(rep(3, times = J/3),
                                              rep(4, times = J/3),
                                              rep(5, times = J/3))))
# List over subpopulation of list of matrices of beta matrices for each item j
beta_list_x_temp <- get_betas_x(profiles = profiles_local, R = R, 
                                modal_theta_prob = modal_theta_prob, 
                                formula_x = formula_x_local, 
                                V_unique = V_unique)
# Add in coefficients for strata s_all, updating formula_x and beta_list_x
formula_x_local <- "~ h_all + s_all + h_all:s_all"
beta_list_x_1_27 <- lapply(1:(J-3), function(j)
  cbind(beta_list_x_temp[[j]], s_all2 = rep(0, 5), `h_all2:s_all2` = rep(0, 5)))
# Items 28-30 are affected as follows: s_all=2 associated with r=3 for h=1 and 
# with r=4 for h=2
beta_list_x_28_30 <- lapply((J-2):J, function(j) 
  cbind(beta_list_x_temp[[j]], s_all2 = c(0, -0.5, 1, 0, 0), 
        `h_all2:s_all2` = c(0, 0, 0, 1, -0.5)))
beta_list_x_local <- c(beta_list_x_1_27, beta_list_x_28_30)

# Simulate population
pop_seed <- 1  # Set seed
sim_pop <- simulate_pop_wrpc_v2(N = N, H = H, S = S, J = J, R = R, K = K, 
                                N_h = N_h, N_s = N_s,
                                modal_theta_prob = modal_theta_prob,
                                formula_c = formula_c, 
                                formula_x_global = formula_x_global, 
                                formula_x_local = formula_x_local, 
                                beta_mat_c = beta_mat_c, 
                                beta_list_x_global = beta_list_x_global,
                                beta_list_x_local = beta_list_x_local,
                                pop_seed = pop_seed, save_res = TRUE, 
                                save_path = paste0(wd, data_dir, "scen_21/"))
# Duplicate for scenario 22
save(sim_pop, file = paste0(wd, data_dir, "scen_22/sim_pop.RData"))

### Stratified sampling only affecting pi
num_samps <- 100
for (i in 1:num_samps) {
  # Scenario 21: Sample size 2000
  sim_samp <- simulate_samp_wrpc(sim_pop = sim_pop, samp_size = 2000, strat = TRUE, 
                                 strat_dist = c(0.3, 0.7), clust = FALSE, 
                                 samp_seed = i, save_res = TRUE, 
                                 save_path = paste0(wd, data_dir, "scen_21/", i))
  
  # Scenario 22: Sample size 400
  sim_samp <- simulate_samp_wrpc(sim_pop = sim_pop, samp_size = 400, strat = TRUE, 
                                 strat_dist = c(0.3, 0.7), clust = FALSE, 
                                 samp_seed = i, save_res = TRUE, 
                                 save_path = paste0(wd, data_dir, "scen_22/", i))
}
# data_vars <- sim_samp

#================== Check simulated values =====================================
# Check distribution of c_i
# Population
prop.table(table(sim_pop$true_Ci))
# Sample
prop.table(table(sim_samp$true_Ci))

# Check cross-tab between s_i and h_i
# Population
table(sim_pop$true_Si, sim_pop$true_Hi)
# Sample 
table(sim_samp$true_Si, sim_samp$true_Hi)

# Observed global thetas
# item 1: 1, 4, 5 (similar in both)
# Population
t(sapply(1:K, function(k) 
  prop.table(table(sim_pop$X_data[sim_pop$true_Ci == k & 
                                    sim_pop$true_Gij[, 1] == 1, 1]))))
# Sample
sapply(1:K, function(k) 
  prop.table(table(sim_samp$X_data[sim_samp$true_Ci == k & 
                                     sim_samp$true_Gij[, 1] == 1, 1])))
# item 30 has association with s_all: 3, 2, 1 (sample worse than population)
# Population
t(sapply(1:K, function(k) 
  prop.table(table(sim_pop$X_data[sim_pop$true_Ci == k & 
                                    sim_pop$true_Gij[, 30] == 1, 30]))))
# Sample: 
sapply(1:K, function(k) 
  prop.table(table(sim_samp$X_data[sim_samp$true_Ci == k & 
                                     sim_samp$true_Gij[, 30] == 1, 30])))


# Observed local thetas
# item 1: 1, 3 (should be similar to population)
# Population
t(sapply(1:H, function(h) 
  prop.table(table(sim_pop$X_data[sim_pop$true_Hi == h & 
                                    sim_pop$true_Gij[, 1] == 0, 1]))))
# Sample: 
t(sapply(1:H, function(h) 
  prop.table(table(sim_samp$X_data[sim_samp$true_Hi == h & 
                                     sim_samp$true_Gij[, 1] == 0, 1]))))
# item 30 has association with s_all: 2, 5 (should be worse in sample)
# Population
t(sapply(1:H, function(h) 
  prop.table(table(sim_pop$X_data[sim_pop$true_Hi == h  & 
                                    sim_pop$true_Gij[, 30] == 0, 30]))))
# Sample
t(sapply(1:H, function(h) 
  prop.table(table(sim_samp$X_data[sim_samp$true_Hi == h & 
                                     sim_samp$true_Gij[, 30] == 0, 30]))))

