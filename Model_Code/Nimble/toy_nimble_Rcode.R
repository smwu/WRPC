#====================================
# Convert the nimble RPC model into R
# 2024/10/27
#====================================

rm(list=ls())

### Workflow
library(tidyverse)
library(ggpubr)
library(gtools)

wd <- "~/Documents/Github/WRPC/"
code_dir <- "Code/Nimble/"
data_dir <- "Code/Nimble/"
res_dir <- "Code/Nimble/"

# Data simulation code
source(paste0(wd, "Code/simulate_data_wrpc.R"))

# # Read in data
# load(paste0(wd, data_dir, "w1_sim_pop.RData"))
# load(paste0(wd, data_dir, "w1_sim_samp.RData"))
# data_vars <- sim_samp

#================= Simulate data ===============================================
### Generate population
# Population size and strata dimensions
N <- 10000; H <- 2; S <- 2; J <- 10; R <- 4 
N_h <- c(5000, 5000)
N_s <- c(6500, 3500)

# Generate C ~ S
K <- 3  
formula_c <- "~ s_all"
V_unique <- data.frame(s_all = as.factor(1:S))
# Corresponds to an overall true_pi ~= (0.253, 0.522, 0.225)
# Matrix of global class assignment probabilities for each level of s_i
pi_global_mat <- matrix(c(0.3, 0.5, 0.2,   # global class probs for s_i=1
                          0.1, 0.6, 0.3),  # global class probs for s_i=2
                        byrow = TRUE, nrow = S, ncol = K)
beta_mat_c <- get_betas_c(pi_mat = pi_global_mat, formula_c = formula_c, 
                          V_unique = V_unique)

# Generate H 
h_all <- c(rep(1, times = N_h[1]), rep(2, times = N_h[2]))

# Global-local allocation probability of global assignment
nu <- matrix(c(0.3, 0.7), byrow = FALSE, nrow = H, ncol = J)

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
                                        C3 = c(rep(3, times = 0.3 * J), 
                                               rep(4, times = 0.4 * J),
                                               rep(1, times = 0.3 * J))))
# List of beta matrices for each item j
beta_list_x_global <- get_betas_x(profiles = profiles_global, R = R, 
                                  modal_theta_prob = modal_theta_prob, 
                                  formula_x = formula_x_global, 
                                  V_unique = V_unique)
# List of matrices of local profiles for each subpopulation
profiles_local <- as.matrix(data.frame(H1 = rep(1, times = J), 
                                       H2 = rep(4, times = J)))
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
                                save_path = paste0(wd, res_dir, "w1_"))

### Stratified sampling only affecting pi
sim_samp <- simulate_samp_wrpc(sim_pop = sim_pop, samp_size = 1000, strat = TRUE, 
                               strat_dist = c(0.3, 0.7), clust = FALSE, 
                               samp_seed = 1, save_res = TRUE, 
                               save_path = paste0(wd, res_dir, "w1_"))
data_vars <- sim_samp

#=======================================

# Define true parameters
parameter_true <- list(
  c_all = sim_pop$true_Ci,
  g_mat = sim_pop$true_Gij,
  pi = as.vector(sim_pop$true_pi_global),
  theta_global = sim_pop$true_global_thetas,
  theta_local = sim_pop$true_local_thetas,
  nu = sim_pop$true_nu
)

# Sample vs pop parameters
samp_inds <- sim_samp$samp_ind
pop_inds <- 1:sim_pop$N
# nu[h=2, j=1]
# population
mean(parameter_true$g_mat[sim_pop$true_Hi == 2, 1])
# sample
mean(parameter_true$g_mat[sim_pop$true_Hi == 2 & pop_inds %in% samp_inds, 1])


#================ (1) Run model ================================================

# Fixed K
K <- 3

# Run WRPC
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
stratum_id <- NULL      # Stratum indicators, nx1
cluster_id <- NULL                   # No clustering
sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
h_all <- data_vars$true_Hi           # Subpopulation indicators, nx1

print("Read in data")
# Obtain dimensions
n <- dim(x_mat)[1]        # Number of individuals
J <- dim(x_mat)[2]        # Number of exposure items
R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
             function(x) length(unique(x)))  
R <- max(R_j)             # Maximum number of exposure categories across items
H <- length(unique(h_all)) # Number of subpopulations
# Obtain normalized weights
kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1

### Set hyperparameters
a <- b <- 1
alpha <- rep(1, K) / K
eta_global <- eta_local <- matrix(0.01, nrow = J, ncol = R)
for (j in 1:J) {
  eta_global[j, 1:R_j[j]] <- rep(1, R_j[j])
  eta_local[j, 1:R_j[j]] <- rep(1, R_j[j])
}

### Initialize model
# Initialize at true values
nu <- data_vars$true_nu
pi <- as.vector(data_vars$true_pi_global)
c_all <- data_vars$true_Ci
theta_global <- data_vars$true_global_thetas
theta_local <- data_vars$true_local_thetas
g_mat <- data_vars$true_Gij

# # Initialize at different values
# nu <- matrix(0.5, nrow = H, ncol = J)
# pi <- rep(1, times = K) / K
# c_all <- rep(1:K, length.out = nrow(x_mat))
# theta_global <- array(1 / R, dim = c(J, K, R))
# theta_local <- array(1 / R, dim = c(J, H, R))
# g_mat <- matrix(rbinom(n = nrow(x_mat), size = 1, prob = 0.5), nrow = nrow(x_mat), ncol = J)


# Initialize at random values
# # Prior for global-local assignment probability nu
# nu <- t(sapply(1:H, function(h) 
#   stats::rbeta(n = J, shape1 = alpha_nu[h], shape2 = beta[h])))
# # Prior for pi
# pi <- c(gtools::rdirichlet(n = 1, alpha = alpha))
# # Initialize global class assignment, c, for individuals
# c_all <- apply(stats::rmultinom(n = n, size = 1, prob = pi), 2, 
#                function(x) which(x == 1))
# # Prior for global thetas
# theta_global <- array(0, dim = c(J, K, R))
# for (j in 1:J) {
#   for (k in 1:K) {
#     theta_global[j, k, ] <- 
#       c(gtools::rdirichlet(n = 1, alpha = eta_global[j, ]))
#   }
# }
# # Prior for local thetas
# theta_local <- array(0, dim = c(J, H, R))
# for (j in 1:J) {
#   for (h in 1:H) {
#     theta_local[j, h, ] <- 
#       c(gtools::rdirichlet(n = 1, alpha = eta_local[j, ]))
#   }
# }
# # Initialize global-local assignment, g, for individuals and items
# g_mat <- matrix(NA, nrow = n, ncol = J)
# for (i in 1:n) {
#   g_mat[i, ] <- rbinom(n = J, size = 1, prob = nu[h_all[i], ])
# }


### Run MCMC
n_runs <- 20000
burn <- 0
thin <- 1
update <- n_runs / 10
switch <- 50
# Number of MCMC iterations to store
n_storage <- ceiling(n_runs / thin) 
# Initialize variables
nu_MCMC <- array(NA, dim = c(n_storage, H, J))
pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
theta_global_MCMC <- array(NA, dim = c(n_storage, J, K, R))
theta_local_MCMC <- array(NA, dim = c(n_storage, J, H, R))
c_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
g_mat_MCMC <- array(NA, dim = c(n_storage, n, J))


# Begin runtime tracker
start_time <- Sys.time()
set.seed(20241027)

### Update parameters and variables
for (m in 1:n_runs) {
  
  ### Update g_mat 
  for (i in 1:n) {
    h_i <- h_all[i]
    p_i <- numeric(J)
    # Posterior global-local assignment probability for each item
    for (j in 1:J) {
      global_prob <- nu[h_i, j] * theta_global[j, c_all[i], x_mat[i, j]]
      local_prob <- (1 - nu[h_i, j]) * theta_local[j, h_i, x_mat[i, j]]
      p_i[j] <- global_prob / (global_prob + local_prob)
    }
    # For each individual, assign items to global or local
    g_mat[i, ] <- rbinom(n = J, size = 1, prob = p_i)
  }
  
  ### Update pi
  alpha_post <- numeric(K)
  for (k in 1:K) {
    # Add sum of normalized weights for those assigned to class k, equiv. to
    # weighted number of individuals assigned to each class
    alpha_post[k] <- alpha[k] + sum(w_all[c_all == k])
  }
  pi <- c(gtools::rdirichlet(n = 1, alpha = alpha_post))
  
  ### Update c_all
  # Posterior class membership probabilities
  pred_global_class_probs <- matrix(NA, nrow = n, ncol = K) 
  # Individual log-likelihood for each class
  log_cond_c <- matrix(NA, nrow = n, ncol = K)       
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K) {
      # Calculate theta component of individual log-likelihood assuming class k
      log_global_theta_comp_k <- 0
      for (j in 1:J) {
        # Add in global thetas for items with global assignment
        if (g_mat[i, j] == 1) {
          log_global_theta_comp_k <- log_global_theta_comp_k + 
            log(theta_global[j, k, x_mat[i, j]])
        }
      }
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi[k]) + log_global_theta_comp_k
    }
    # Calculate p(c_i=k|-) = p(x,c_i=k) / p(x)
    pred_global_class_probs[i, ] <- 
      exp(log_cond_c[i, ] - matrixStats::logSumExp(log_cond_c[i, ], na.rm = TRUE))
    # Update class assignment using the posterior probabilities
    c_all[i] <- which(stats::rmultinom(n = 1, size = 1, 
                                       prob = pred_global_class_probs[i, ]) == 1)
  }
  
  ### Update theta_global
  eta_global_post <- numeric(R)
  for (j in 1:J) {
    for (k in 1:K) {
      for (r in 1:R) {
        # Add sum of normalized weights for those assigned to class k with x_ij = r
        eta_global_post[r] <- eta_global[r] + 
          sum(w_all[(g_mat[, j] == 1) & (c_all == k) & (x_mat[, j] == r)])
      }
      theta_global[j, k, ] <- c(gtools::rdirichlet(n = 1, alpha = eta_global_post))
    }
  }
  
  ### Update theta_global
  eta_local_post <- numeric(R)
  for (j in 1:J) {
    for (h in 1:H) {
      for (r in 1:R) {
        # Add sum of normalized weights for those assigned to class k with x_ij = r
        eta_local_post[r] <- eta_local[r] + 
          sum(w_all[(g_mat[, j] == 0) & (h_all == h) & (x_mat[, j] == r)])
      }
      theta_local[j, h, ] <- c(gtools::rdirichlet(n = 1, alpha = eta_local_post))
    }
  }
  
  ### Update nu
  for (h in 1:H) {
    for (j in 1:J) {
      a_post <- a + sum(w_all[(g_mat[, j] == 1) & (h_all == h)])
      b_post <- b + sum(w_all[(g_mat[, j] == 0) & (h_all == h)])
      nu[h, j] <- stats::rbeta(n = 1, shape1 = a_post, shape2 = b_post)
    }
  }
  
  ### Store posterior values based on thinning 
  if (m %% thin == 0) {
    m_thin <- m / thin
    nu_MCMC[m_thin, , ] <- nu
    g_mat_MCMC[m_thin, , ] <- g_mat
    pi_MCMC[m_thin, ] <- pi
    theta_global_MCMC[m_thin, , , ] <- theta_global
    theta_local_MCMC[m_thin, , , ] <- theta_local
    c_all_MCMC[m_thin, ] <- c_all
  }
  
  ### Relabel classes every 10 iterations to encourage mixing 
  if (m %% switch == 0) {
    new_order_c <- gtools::permute(1:K)      # New permuted labels
    new_c_all <- numeric(K)        # Class assignments with new labels
    for (k in 1:K) {
      new_c_all[c_all == k] <- new_order_c[k]
    }
    c_all <- new_c_all             # Relabel class assignments
    pi <- pi[new_order_c]            # Relabel class probabilities
    theta_global <- theta_global[, new_order_c, ]  # Relabel item category probabilities
  }
  
  # Print out progress 
  if (m %% update == 0) {
    print(paste0("Iteration ", m, " completed!"))
  }
}

# Discard burn-in
warmup <- ceiling(burn / thin)
if (warmup > 0) {
  nu_MCMC <- nu_MCMC[-(1:warmup), , ]
  pi_MCMC <- pi_MCMC[-(1:warmup), ]
  theta_global_MCMC <- theta_global_MCMC[-(1:warmup), , , ]
  theta_local_MCMC <- theta_local_MCMC[-(1:warmup), , , ]
  c_all_MCMC <- c_all_MCMC[-(1:warmup), ]
  g_mat_MCMC <- g_mat_MCMC[-(1:warmup), , ]
}

MCMC_out <- list(nu_MCMC = nu_MCMC, pi_MCMC = pi_MCMC, 
                 theta_global_MCMC = theta_global_MCMC, 
                 theta_local_MCMC = theta_local_MCMC,
                 c_all_MCMC = c_all_MCMC, g_mat_MCMC = g_mat_MCMC)

# End runtime tracker
end_time <- Sys.time()


### Post-processing to address label switching
class_cutoff_global <- 0.05
### Global classes
# Get median number of classes with >= cutoff% of individuals, over all iterations
M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
K_med <- round(stats::median(rowSums(MCMC_out$pi_MCMC >= class_cutoff_global)))

# Cluster individuals into reduced number of classes using agglomerative clustering
# Calculate pairwise distance matrix using Hamming distance: proportion of
# iterations where two individuals have differing class assignments
distMat <- e1071::hamming.distance(t(MCMC_out$c_all_MCMC))
# Hierarchical clustering dendrogram
dendrogram_global <- stats::hclust(stats::as.dist(distMat), method = "complete") 
# Group individuals into K_med classes
red_c_all <- stats::cutree(dendrogram_global, k = K_med)        
# Modify classes if any classes are less than the cutoff percentage
class_prop <- prop.table(table(red_c_all))
if (any(class_prop < class_cutoff_global)) {
  # Get classes that are too small
  small <- which(class_prop < class_cutoff_global)
  # Group individuals into a larger number of classes 
  red_c_all_temp <- stats::cutree(dendrogram, k = K_med + length(small))
  red_c_all <- red_c_all_temp
  class_prop_temp <- prop.table(table(red_c_all_temp))
  # Get updated classes that are too small
  small_temp <- sort(which(class_prop_temp < class_cutoff_global))
  for (small_c in 1:length(small_temp)) {
    c_ind <- small_temp[small_c]
    class_small <- which(red_c_all_temp == c_ind)
    # Get nearest class
    inds <- 1:length(class_prop_temp)
    class_dist <- sapply(inds, function(x) 
      mean(distMat[class_small, which(red_c_all_temp == x)]))
    # Set small class distance to Inf
    class_dist[small_temp] <- Inf
    nearest <- which.min(class_dist[-c_ind])
    red_c_all[red_c_all_temp == c_ind] <- nearest
  }
  class_prop <- prop.table(table(red_c_all))
  K_med <- length(class_prop)
  # # Check class sizes
  # prop.table(table(red_c_all))
}
# Get unique reduced classes to aid relabeling
unique_red_classes <- unique(red_c_all)

# For each iteration, relabel new classes using the most common old class assignment
relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
for (k in 1:K_med) {
  red_class <- unique_red_classes[k]
  relabel_red_classes[, k] <- 
    apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == red_class]), 1, get_mode)
}

# Reduce and reorder parameter estimates using new classes
pi <- matrix(NA, nrow = M, ncol = K_med)
theta_global <- array(NA, dim = c(M, J, K_med, R))
for (m in 1:M) {
  iter_order <- relabel_red_classes[m, ]
  pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
  pi[m, ] <- pi_temp / sum(pi_temp)
  theta_global[m, , , ] <- MCMC_out$theta_global_MCMC[m, , iter_order, ]
}

post_MCMC_out <- list(K_med = K_med, pi = pi, 
                      theta_global = theta_global, 
                      dendrogram_global = dendrogram_global)


### Get estimates
## Identify unique classes using modal exposure categories 
# Posterior median estimate for theta across iterations
theta_global_med_temp <- apply(post_MCMC_out$theta_global, c(2, 3, 4), 
                               function(x) stats::median(x, na.rm = TRUE))
# Posterior modal exposure categories for each exposure item and reduced class
theta_global_modes <- apply(theta_global_med_temp, c(1, 2), which.max)
# Identify unique classes
unique_global_classes <- which(!duplicated(theta_global_modes, MARGIN = 2))
# Number of unique classes
K_red <- length(unique_global_classes)

## Use new classes to adjust and re-normalize posterior samples 
# Combine duplicated classes and re-normalize pi to sum to 1
M <- dim(post_MCMC_out$pi)[1]                # Number of iterations
## Global pi
# Obtain pi for unique classes
pi_red <- post_MCMC_out$pi[, unique_global_classes, drop = FALSE] 
# Reduce classes if there are duplicated classes
if (K_red < dim(post_MCMC_out$pi)[2]) {  
  for (k in 1:K_red) {
    # Find duplicated modal theta patterns
    dupes_k <- apply(theta_global_modes, 2, function(x)
      identical(x,theta_global_modes[, unique_global_classes[k]]))
    # Combine class proportions for all duplicated patterns together
    pi_red[, k] <- 
      apply(as.matrix(post_MCMC_out$pi[, dupes_k]), 1, sum, na.rm = TRUE)
  }
}
# Re-normalize to ensure pi sums to 1 for each iteration
pi_red = pi_red / rowSums(pi_red, na.rm = TRUE)

# Get posterior parameter samples for unique classes for theta 
# Global theta: MxJxKxR
theta_global_red <- post_MCMC_out$theta_global[, , unique_global_classes, , 
                                               drop = FALSE]
theta_global_red <- plyr::aaply(theta_global_red, c(1, 2, 3), function(x) 
  x / sum(x, na.rm = TRUE), .drop = FALSE) # Re-normalize
# Local theta: MxJxHxR
theta_local_red <- plyr::aaply(MCMC_out$theta_local_MCMC, c(1, 2, 3), function(x) 
  x / sum(x, na.rm = TRUE), .drop = FALSE) # Re-normalize

## Posterior median estimates 
# Pi
pi_med <- apply(pi_red, 2, function(x) 
  stats::median(x, na.rm = TRUE))
pi_med <- pi_med / sum(pi_med, na.rm = TRUE)  # Re-normalize

# Theta
theta_global_med <- apply(theta_global_red, c(2, 3, 4), function(x) 
  stats::median(x, na.rm = TRUE))
theta_global_med <- plyr::aaply(theta_global_med, c(1, 2), function(x) 
  x / sum(x, na.rm = TRUE), .drop = FALSE)  # Re-normalize
theta_local_med <- apply(theta_local_red, c(2, 3, 4), function(x) 
  stats::median(x, na.rm = TRUE))
theta_local_med <- plyr::aaply(theta_local_med, c(1, 2), function(x) 
  x / sum(x, na.rm = TRUE), .drop = FALSE)  # Re-normalize
# Theta modes
theta_global_modes <- apply(theta_global_med, c(1, 2), which.max)
theta_local_modes <- apply(theta_local_med, c(1, 2), which.max)

# Nu
nu_red <- MCMC_out$nu_MCMC
nu_med <- apply(nu_red, c(2, 3), function(x) stats::median(x, na.rm = TRUE))

### Update c using unique classes and posterior estimates 
# Global-local assignments to use for obtaining c_all and l_all
g_mat <- MCMC_out$g_mat_MCMC[M, , ]

# Posterior class membership probabilities
pred_global_class_probs <- matrix(NA, nrow = n, ncol = K_red) 
# Individual log-likelihood for each class
log_cond_c <- matrix(NA, nrow = n, ncol = K_red)       
# Calculate posterior class membership, p(c_i=k|-), for each class k
for (i in 1:n) {
  for (k in 1:K_red) {
    # Calculate theta component of individual log-likelihood assuming class k
    log_global_theta_comp_k <- 0
    for (j in 1:J) {
      # Add in global thetas for items with global assignment
      if (g_mat[i, j] == 1) {
        log_global_theta_comp_k <- log_global_theta_comp_k + 
          log(theta_global_med[j, k, x_mat[i, j]])
      }
    }
    # Individual log-likelihood for class k
    log_cond_c[i, k] <- log(pi_med[k]) + log_global_theta_comp_k
  }
  # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
  pred_global_class_probs[i, ] <- 
    exp(log_cond_c[i, ] - matrixStats::logSumExp(log_cond_c[i, ], na.rm = TRUE))
  # Update class assignment using the posterior probabilities
  c_all[i] <- which(stats::rmultinom(n = 1, size = 1, 
                                     prob = pred_global_class_probs[i, ]) == 1)
}

## Update g using posterior estimates
# Final global-local assignments for each item
for (i in 1:n) {
  h_i <- h_all[i]
  p_i <- numeric(J)
  # Posterior global-local assignment probability for each item
  for (j in 1:J) {
    global_prob <- nu_med[h_i, j] * theta_global_med[j, c_all[i], x_mat[i, j]]
    local_prob <- (1 - nu_med[h_i, j]) * theta_local_med[j, h_i, x_mat[i, j]]
    p_i[j] <- global_prob / (global_prob + local_prob)
  }
  # For each individual, assign items to global or local
  g_mat[i, ] <- rbinom(n = J, size = 1, prob = p_i)
}

estimates <- list(K_red = K_red, pi_red = pi_red, 
                  theta_global_red = theta_global_red,
                  theta_local_red = theta_local_red, 
                  pi_med = pi_med, 
                  theta_global_med = theta_global_med, 
                  theta_local_med = theta_local_med, 
                  theta_global_modes = theta_global_modes,
                  theta_local_modes = theta_local_modes, c_all = c_all, 
                  pred_global_class_probs = pred_global_class_probs, 
                  g_mat = g_mat, nu_med = nu_med, nu_red = nu_red)
# save(estimates, file = paste0(wd, res_dir, "res_10000_diff_init.RData"))
save(estimates, file = paste0(wd, res_dir, "res_w1_10000.RData"))


#================ (2) Examine results ==========================================
ordering <- c(2, 1, 3)
theta_global_plot <- theta_global_red[, , ordering, ]
theta_local_plot <- theta_local_red
pi_plot <- pi_red[, ordering]
nu_plot <- nu_red
c_all_plot <- c_all_MCMC
g_mat_plot <- g_mat_MCMC
c_plot <- ifelse(c_all == ordering[1], 1, ifelse(c_all == ordering[2], 2, 3))
# par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))


## plot trace of theta_global
j <- 1
r <- 1
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Item probability: item ", j, ", global theta"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  plot(theta_global_plot[, j, k, r], type = "l", ylim = c(0,1), xlab = "", ylab = "")
  abline(h = parameter_true$theta_global[j, k, r], col = "red")
}
theta_global_modes

## plot trace of theta_local
j <- 1
r <- 1
mult.fig(mfrow = c(H, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Item probability: item ", j, ", local theta"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (h in 1:H) {
  plot(theta_local_plot[, j, h, r], type = "l", ylim = c(0,1), xlab = "", ylab = "")
  abline(h = parameter_true$theta_local[j, h, r], col = "red")
}
theta_local_modes
summary(theta_local_plot[, 1, 2, 1])

## plot trace of pi
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for pi"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  plot(pi_plot[, k], type = "l", ylim = c(0,1), xlab = "", ylab = "")
  abline(h = parameter_true$pi[k], col = "red")
}

## plot trace of nu
j <- 1
mult.fig(mfrow = c(H, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for nu"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (h in 1:H) {
  plot(nu_plot[, h, j], type = "l", ylim = c(0,1), xlab = "", ylab = "")
  abline(h = parameter_true$nu[h, j], col = "red")
}
summary(nu_red[, 2, 1])

## look at c_all --------
index <- c(1, 500, 1000)
mult.fig(mfrow = c(3, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main = paste("c_all assignment"),
         cex = 0.8, marP = - c(0, 1, 2, 0))
for (i in index) {
  plot(c_all_plot[, i], type = "l", xlab = "", ylab = "")
  abline(h = parameter_true$c_all[i], col = "red")
}
prop.table(table(c_all_plot[, 1]))
prop.table(table(c_all_plot[, 500]))
prop.table(table(c_all_plot[, 1000]))
# Check misclassification
table(sim_samp$true_Ci, c_plot)
prop.table(table(sim_samp$true_Ci, c_plot))

## look at g_mat --------
i <- 1
var_index <- c(1, 5, 10)
mult.fig(mfrow = c(3, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main = paste("g_mat assignment"),
         cex = 0.8, marP = - c(0, 1, 2, 0))
for (j in var_index) {
  plot(g_mat_plot[, i, j], type = "l", xlab = "", ylab = "")
  abline(h = parameter_true$g_mat[i, j], col = "red")
}
prop.table(table(g_mat_plot[, 1, 1]))
prop.table(table(g_mat_plot[, 1, 5]))
prop.table(table(g_mat_plot[, 1, 10]))
# Check misclassification
table(c(sim_samp$true_Gij), c(g_mat))
prop.table(table(c(sim_samp$true_Gij), c(g_mat)))


#===================== Run adaptive sampler ====================================

# Overfitted K
K <- 30

# Run WRPC
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
stratum_id <- NULL      # Stratum indicators, nx1
cluster_id <- NULL                   # No clustering
sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
h_all <- data_vars$true_Hi           # Subpopulation indicators, nx1

print("Read in data")
# Obtain dimensions
n <- dim(x_mat)[1]        # Number of individuals
J <- dim(x_mat)[2]        # Number of exposure items
R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
             function(x) length(unique(x)))  
R <- max(R_j)             # Maximum number of exposure categories across items
H <- length(unique(h_all)) # Number of subpopulations
# Obtain normalized weights
kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1

### Set hyperparameters
a <- b <- 1
alpha <- rep(1, K) / K
eta_global <- eta_local <- matrix(0.01, nrow = J, ncol = R)
for (j in 1:J) {
  eta_global[j, 1:R_j[j]] <- rep(1, R_j[j])
  eta_local[j, 1:R_j[j]] <- rep(1, R_j[j])
}

### Initialize model
# # Initialize at different values
# nu <- matrix(0.5, nrow = H, ncol = J)
# pi <- rep(1, times = K) / K
# c_all <- rep(1:K, length.out = nrow(x_mat))
# theta_global <- array(1 / R, dim = c(J, K, R))
# theta_local <- array(1 / R, dim = c(J, H, R))
# g_mat <- matrix(rbinom(n = nrow(x_mat), size = 1, prob = 0.5), nrow = nrow(x_mat), ncol = J)


## Initialize at random values
# Prior for global-local assignment probability nu
nu <- t(sapply(1:H, function(h)
  stats::rbeta(n = J, shape1 = a, shape2 = b)))
# Prior for pi
pi <- c(gtools::rdirichlet(n = 1, alpha = alpha))
# Initialize global class assignment, c, for individuals
c_all <- apply(stats::rmultinom(n = n, size = 1, prob = pi), 2,
               function(x) which(x == 1))
# Prior for global thetas
theta_global <- array(0, dim = c(J, K, R))
for (j in 1:J) {
  for (k in 1:K) {
    theta_global[j, k, ] <-
      c(gtools::rdirichlet(n = 1, alpha = eta_global[j, ]))
  }
}
# Prior for local thetas
theta_local <- array(0, dim = c(J, H, R))
for (j in 1:J) {
  for (h in 1:H) {
    theta_local[j, h, ] <-
      c(gtools::rdirichlet(n = 1, alpha = eta_local[j, ]))
  }
}
# Initialize global-local assignment, g, for individuals and items
g_mat <- matrix(NA, nrow = n, ncol = J)
for (i in 1:n) {
  g_mat[i, ] <- rbinom(n = J, size = 1, prob = nu[h_all[i], ])
}


### Run MCMC
n_runs <- 2000
burn <- 0
thin <- 1
update <- n_runs / 10
switch <- 50
# Number of MCMC iterations to store
n_storage <- ceiling(n_runs / thin) 
# Initialize variables
nu_MCMC <- array(NA, dim = c(n_storage, H, J))
pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
theta_global_MCMC <- array(NA, dim = c(n_storage, J, K, R))
theta_local_MCMC <- array(NA, dim = c(n_storage, J, H, R))
c_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
g_mat_MCMC <- array(NA, dim = c(n_storage, n, J))


# Begin runtime tracker
start_time <- Sys.time()
set.seed(20241027)

### Update parameters and variables
for (m in 1:n_runs) {
  
  ### Update g_mat 
  for (i in 1:n) {
    h_i <- h_all[i]
    p_i <- numeric(J)
    # Posterior global-local assignment probability for each item
    for (j in 1:J) {
      global_prob <- nu[h_i, j] * theta_global[j, c_all[i], x_mat[i, j]]
      local_prob <- (1 - nu[h_i, j]) * theta_local[j, h_i, x_mat[i, j]]
      p_i[j] <- global_prob / (global_prob + local_prob)
    }
    # For each individual, assign items to global or local
    g_mat[i, ] <- rbinom(n = J, size = 1, prob = p_i)
  }
  
  ### Update pi
  alpha_post <- numeric(K)
  for (k in 1:K) {
    # Add sum of normalized weights for those assigned to class k, equiv. to
    # weighted number of individuals assigned to each class
    alpha_post[k] <- alpha[k] + sum(w_all[c_all == k])
  }
  pi <- c(gtools::rdirichlet(n = 1, alpha = alpha_post))
  
  ### Update c_all
  # Posterior class membership probabilities
  pred_global_class_probs <- matrix(NA, nrow = n, ncol = K) 
  # Individual log-likelihood for each class
  log_cond_c <- matrix(NA, nrow = n, ncol = K)       
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K) {
      # Calculate theta component of individual log-likelihood assuming class k
      log_global_theta_comp_k <- 0
      for (j in 1:J) {
        # Add in global thetas for items with global assignment
        if (g_mat[i, j] == 1) {
          log_global_theta_comp_k <- log_global_theta_comp_k + 
            log(theta_global[j, k, x_mat[i, j]])
        }
      }
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi[k]) + log_global_theta_comp_k
    }
    # Calculate p(c_i=k|-) = p(x,c_i=k) / p(x)
    pred_global_class_probs[i, ] <- 
      exp(log_cond_c[i, ] - matrixStats::logSumExp(log_cond_c[i, ], na.rm = TRUE))
    # Update class assignment using the posterior probabilities
    c_all[i] <- which(stats::rmultinom(n = 1, size = 1, 
                                       prob = pred_global_class_probs[i, ]) == 1)
  }
  
  ### Update theta_global
  eta_global_post <- numeric(R)
  for (j in 1:J) {
    for (k in 1:K) {
      for (r in 1:R) {
        # Add sum of normalized weights for those assigned to class k with x_ij = r
        eta_global_post[r] <- eta_global[j, r] + 
          sum(w_all[(g_mat[, j] == 1) & (c_all == k) & (x_mat[, j] == r)])
      }
      theta_global[j, k, ] <- c(gtools::rdirichlet(n = 1, alpha = eta_global_post))
    }
  }
  
  ### Update theta_global
  eta_local_post <- numeric(R)
  for (j in 1:J) {
    for (h in 1:H) {
      for (r in 1:R) {
        # Add sum of normalized weights for those assigned to class k with x_ij = r
        eta_local_post[r] <- eta_local[j, r] + 
          sum(w_all[(g_mat[, j] == 0) & (h_all == h) & (x_mat[, j] == r)])
      }
      theta_local[j, h, ] <- c(gtools::rdirichlet(n = 1, alpha = eta_local_post))
    }
  }
  
  ### Update nu
  for (h in 1:H) {
    for (j in 1:J) {
      a_post <- a + sum(w_all[(g_mat[, j] == 1) & (h_all == h)])
      b_post <- b + sum(w_all[(g_mat[, j] == 0) & (h_all == h)])
      nu[h, j] <- stats::rbeta(n = 1, shape1 = a_post, shape2 = b_post)
    }
  }
  
  ### Store posterior values based on thinning 
  if (m %% thin == 0) {
    m_thin <- m / thin
    nu_MCMC[m_thin, , ] <- nu
    g_mat_MCMC[m_thin, , ] <- g_mat
    pi_MCMC[m_thin, ] <- pi
    theta_global_MCMC[m_thin, , , ] <- theta_global
    theta_local_MCMC[m_thin, , , ] <- theta_local
    c_all_MCMC[m_thin, ] <- c_all
  }
  
  ### Relabel classes every 10 iterations to encourage mixing 
  if (m %% switch == 0) {
    new_order_c <- gtools::permute(1:K)      # New permuted labels
    new_c_all <- numeric(K)        # Class assignments with new labels
    for (k in 1:K) {
      new_c_all[c_all == k] <- new_order_c[k]
    }
    c_all <- new_c_all             # Relabel class assignments
    pi <- pi[new_order_c]            # Relabel class probabilities
    theta_global <- theta_global[, new_order_c, ]  # Relabel item category probabilities
  }
  
  # Print out progress 
  if (m %% update == 0) {
    print(paste0("Iteration ", m, " completed!"))
  }
}

# Discard burn-in
warmup <- ceiling(burn / thin)
if (warmup > 0) {
  nu_MCMC <- nu_MCMC[-(1:warmup), , ]
  pi_MCMC <- pi_MCMC[-(1:warmup), ]
  theta_global_MCMC <- theta_global_MCMC[-(1:warmup), , , ]
  theta_local_MCMC <- theta_local_MCMC[-(1:warmup), , , ]
  c_all_MCMC <- c_all_MCMC[-(1:warmup), ]
  g_mat_MCMC <- g_mat_MCMC[-(1:warmup), , ]
}

MCMC_out <- list(nu_MCMC = nu_MCMC, pi_MCMC = pi_MCMC, 
                 theta_global_MCMC = theta_global_MCMC, 
                 theta_local_MCMC = theta_local_MCMC,
                 c_all_MCMC = c_all_MCMC, g_mat_MCMC = g_mat_MCMC)

# End runtime tracker
end_time <- Sys.time()


### Global classes
class_cutoff_global <- 0.05
# Get median number of classes with >= cutoff% of individuals, over all iterations
M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
K_med <- round(stats::median(rowSums(MCMC_out$pi_MCMC >= class_cutoff_global)))

# Check estimated number of classes
K_med
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
hist(rowSums(MCMC_out$pi_MCMC >= class_cutoff_global))

#============== Plotting sanity checks
### Post-processing to address label switching
# Cluster individuals into reduced number of classes using agglomerative clustering
# Calculate pairwise distance matrix using Hamming distance: proportion of
# iterations where two individuals have differing class assignments
distMat <- e1071::hamming.distance(t(MCMC_out$c_all_MCMC))
# Hierarchical clustering dendrogram
dendrogram_global <- stats::hclust(stats::as.dist(distMat), method = "complete") 
# Group individuals into K_med classes
red_c_all <- stats::cutree(dendrogram_global, k = K_med)        
# Modify classes if any classes are less than the cutoff percentage
class_prop <- prop.table(table(red_c_all))
if (any(class_prop < class_cutoff_global)) {
  # Get classes that are too small
  small <- which(class_prop < class_cutoff_global)
  # Group individuals into a larger number of classes 
  red_c_all_temp <- stats::cutree(dendrogram, k = K_med + length(small))
  red_c_all <- red_c_all_temp
  class_prop_temp <- prop.table(table(red_c_all_temp))
  # Get updated classes that are too small
  small_temp <- sort(which(class_prop_temp < class_cutoff_global))
  for (small_c in 1:length(small_temp)) {
    c_ind <- small_temp[small_c]
    class_small <- which(red_c_all_temp == c_ind)
    # Get nearest class
    inds <- 1:length(class_prop_temp)
    class_dist <- sapply(inds, function(x) 
      mean(distMat[class_small, which(red_c_all_temp == x)]))
    # Set small class distance to Inf
    class_dist[small_temp] <- Inf
    nearest <- which.min(class_dist[-c_ind])
    red_c_all[red_c_all_temp == c_ind] <- nearest
  }
  class_prop <- prop.table(table(red_c_all))
  K_med <- length(class_prop)
  # # Check class sizes
  # prop.table(table(red_c_all))
}
# Get unique reduced classes to aid relabeling
unique_red_classes <- unique(red_c_all)

# For each iteration, relabel new classes using the most common old class assignment
relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
for (k in 1:K_med) {
  red_class <- unique_red_classes[k]
  relabel_red_classes[, k] <- 
    apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == red_class]), 1, get_mode)
}

# Reduce and reorder parameter estimates using new classes
pi <- matrix(NA, nrow = M, ncol = K_med)
theta_global <- array(NA, dim = c(M, J, K_med, R))
for (m in 1:M) {
  iter_order <- relabel_red_classes[m, ]
  pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
  pi[m, ] <- pi_temp / sum(pi_temp)
  theta_global[m, , , ] <- MCMC_out$theta_global_MCMC[m, , iter_order, ]
}

post_MCMC_out <- list(K_med = K_med, pi = pi, 
                      theta_global = theta_global, 
                      dendrogram_global = dendrogram_global)

## plot trace of pi for subsetted classes that are large
ordering <- c(3, 1, 2)
pi_plot <- post_MCMC_out$pi[, ordering]
mult.fig(mfrow = c(K_med, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for pi"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K_med) {
  plot(pi_plot[, k], type = "l", ylim = c(0,1), xlab = "", ylab = "")
  abline(h = parameter_true$pi[k], col = "red")
}
