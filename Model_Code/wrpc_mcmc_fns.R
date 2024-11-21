#==============================================
# WRPC MCMC Helper Functions
# Helper function for running WRPC MCMC sampler
# Author: Stephanie Wu
# Date Updated: 2024/11/02
#==============================================


init_WRPC <- function(K, n, J, R, H, h_all, alpha, eta_global, eta_local, a, b) {
  # Prior for global-local assignment probability nu
  nu <- matrix(NA, nrow = H, ncol = J)
  for (h in 1:H) {
    for (j in 1:J) {
      nu[h, j] <- stats::rbeta(n = 1, shape1 = a[h, j], shape2 = b[h, j])
    }
  }
  # nu <- t(sapply(1:H, function(h)
  #   stats::rbeta(n = J, shape1 = a, shape2 = b)))
  # # Debugging
  # print("nu: ")
  # print(nu)
  
  # Initialize global-local assignment, g, for individuals and items
  g_mat <- matrix(NA, nrow = n, ncol = J)
  for (i in 1:n) {
    g_mat[i, ] <- rbinom(n = J, size = 1, prob = nu[h_all[i], ])
    # # Debugging
    # if (sum(is.na(g_mat[i])) > 0) {
    #   print(nu[h_all[i], ])
    # }
  }
  
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
        c(gtools::rdirichlet(n = 1, alpha = eta_local[j, h, ]))
    }
  }
  
  # Return parameters
  RPC_params <- list(pi = pi, theta_global = theta_global, 
                     theta_local = theta_local, nu = nu, c_all = c_all, 
                     g_mat = g_mat)
  return(RPC_params)
}



run_MCMC_WRPC <- function(RPC_params, n, K, H, J, R, w_all, h_all, x_mat, g_mat,
                         alpha, eta_global, eta_local, a, b, n_runs, burn, thin, 
                         update, switch) {
  
  # Number of MCMC iterations to store
  n_storage <- floor(n_runs / thin) 
  
  # Initialize variables
  nu_MCMC <- array(NA, dim = c(n_storage, H, J))
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_global_MCMC <- array(NA, dim = c(n_storage, J, K, R))
  theta_local_MCMC <- array(NA, dim = c(n_storage, J, H, R))
  c_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
  g_mat_MCMC <- array(NA, dim = c(n_storage, n, J))
  
  # Initialized values
  nu <- RPC_params$nu
  pi <- RPC_params$pi
  theta_global <- RPC_params$theta_global
  theta_local <- RPC_params$theta_local
  c_all <- RPC_params$c_all
  g_mat <- RPC_params$g_mat
  
  ### Update parameters and variables
  for (m in 1:n_runs) {
    # print(paste0("iter: ", m))
    
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
      # # Debugging
      # if (sum(is.na(g_mat[i])) > 0) {
      #   print(pi)
      # }
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
    
    ### Update theta_local
    eta_local_post <- numeric(R)
    for (j in 1:J) {
      for (h in 1:H) {
        for (r in 1:R) {
          # Add sum of normalized weights for those assigned to class k with x_ij = r
          eta_local_post[r] <- eta_local[j, h, r] + 
            sum(w_all[(g_mat[, j] == 0) & (h_all == h) & (x_mat[, j] == r)])
        }
        theta_local[j, h, ] <- c(gtools::rdirichlet(n = 1, alpha = eta_local_post))
      }
    }
    
    ### Update nu
    for (h in 1:H) {
      for (j in 1:J) {
        a_post <- a[h, j] + sum(w_all[(g_mat[, j] == 1) & (h_all == h)])
        b_post <- b[h, j] + sum(w_all[(g_mat[, j] == 0) & (h_all == h)])
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
    
    ### Relabel classes every couple of iterations to encourage mixing 
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
  warmup <- floor(burn / thin)
  if (warmup > 0) {
    nu_MCMC <- nu_MCMC[-(1:warmup), , , drop = FALSE]
    pi_MCMC <- pi_MCMC[-(1:warmup), , drop = FALSE]
    theta_global_MCMC <- theta_global_MCMC[-(1:warmup), , , , drop = FALSE]
    theta_local_MCMC <- theta_local_MCMC[-(1:warmup), , , , drop = FALSE]
    c_all_MCMC <- c_all_MCMC[-(1:warmup), ]
    g_mat_MCMC <- g_mat_MCMC[-(1:warmup), , , drop = FALSE]
  }
  
  # Return samples
  MCMC_out <- list(nu_MCMC = nu_MCMC, pi_MCMC = pi_MCMC, 
                   theta_global_MCMC = theta_global_MCMC, 
                   theta_local_MCMC = theta_local_MCMC,
                   c_all_MCMC = c_all_MCMC, g_mat_MCMC = g_mat_MCMC)
  return(MCMC_out)
}


run_MCMC_Rcpp_WRPC <- function(RPC_params, n, K, H, J, R, w_all, h_all, x_mat, 
                               g_mat, alpha, eta_global, eta_local, a, b, n_runs, 
                               burn, thin, update, switch) {
  
  # Number of MCMC iterations to store
  n_storage <- floor(n_runs / thin) 
  
  # Initialize variables
  nu_MCMC <- array(NA, dim = c(n_storage, H, J))
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_global_MCMC <- array(NA, dim = c(n_storage, J, K, R))
  theta_local_MCMC <- array(NA, dim = c(n_storage, J, H, R))
  c_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
  g_mat_MCMC <- array(NA, dim = c(n_storage, n, J))
  
  # Initialized values
  nu <- RPC_params$nu
  pi <- RPC_params$pi
  theta_global <- RPC_params$theta_global
  theta_local <- RPC_params$theta_local
  c_all <- RPC_params$c_all
  g_mat <- RPC_params$g_mat
  
  ### Update parameters and variables
  for (m in 1:n_runs) {
    # print(paste0("iter: ", m))
    
    ### Update g_mat 
    g_mat <- update_g_WRPC(g_mat = g_mat, J = J, n = n, h_all = h_all, c_all = c_all,
                            nu = nu, x_mat = x_mat, theta_global = theta_global,
                            theta_local = theta_local)
    
    ### Update pi
    pi <- update_pi_WRPC(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
    
    ### Update c_all
    c_all <- update_c_WRPC(c_all = c_all, n = n, K = K, J = J,
                           theta_global = theta_global, x_mat = x_mat,
                           pi = pi, g_mat = g_mat)
    
    ### Update theta_global
    theta_global <- update_theta_global_WRPC(theta_global = theta_global, J = J,
                                             K = K, R = R, eta_global = eta_global,
                                             w_all = w_all, c_all = c_all,
                                             x_mat = x_mat, g_mat = g_mat)
    
    ### Update theta_local
    theta_local <- update_theta_local_WRPC(theta_local = theta_local, H = H,
                                            J = J, R = R, eta_local = eta_local,
                                            w_all = w_all, h_all = h_all,
                                            x_mat = x_mat, g_mat = g_mat)
    
    ### Update nu
    nu <- update_nu_WRPC(nu = nu, H = H, J = J, n = n, a = a, b = b, 
                         w_all = w_all, h_all = h_all, g_mat = g_mat)
    
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
    
    ### Relabel classes every couple of iterations to encourage mixing 
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
  warmup <- floor(burn / thin)
  if (warmup > 0) {
    nu_MCMC <- nu_MCMC[-(1:warmup), , , drop = FALSE]
    pi_MCMC <- pi_MCMC[-(1:warmup), , drop = FALSE]
    theta_global_MCMC <- theta_global_MCMC[-(1:warmup), , , , drop = FALSE]
    theta_local_MCMC <- theta_local_MCMC[-(1:warmup), , , , drop = FALSE]
    c_all_MCMC <- c_all_MCMC[-(1:warmup), ]
    g_mat_MCMC <- g_mat_MCMC[-(1:warmup), , , drop = FALSE]
  }
  
  # Return samples
  MCMC_out <- list(nu_MCMC = nu_MCMC, pi_MCMC = pi_MCMC, 
                   theta_global_MCMC = theta_global_MCMC, 
                   theta_local_MCMC = theta_local_MCMC,
                   c_all_MCMC = c_all_MCMC, g_mat_MCMC = g_mat_MCMC)
  return(MCMC_out)
}


post_process_WRPC <- function(MCMC_out, J, R, H, class_cutoff_global) {
  
  ### Global classes
  # Get median number of classes with >= cutoff% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(stats::median(rowSums(MCMC_out$pi_MCMC >= class_cutoff_global)))
  
  ### Cluster individuals into reduced number of classes using agglomerative clustering
  # Calculate pairwise distance matrix using Hamming distance: proportion of
  # iterations where two individuals have differing class assignments
  distMat <- e1071::hamming.distance(t(MCMC_out$c_all_MCMC))
  # Hierarchical clustering dendrogram
  dendrogram_global <- stats::hclust(stats::as.dist(distMat), method = "complete") 
  # Group individuals into K_med classes
  red_c_all <- stats::cutree(dendrogram_global, k = K_med)       
  
  ### Modify classes if any classes are less than the cutoff percentage
  class_prop <- prop.table(table(red_c_all))
  if (any(class_prop < class_cutoff_global)) {
    # Get classes that are too small
    small <- which(class_prop < class_cutoff_global)
    # Group individuals into a larger number of classes 
    red_c_all_temp <- stats::cutree(dendrogram_global, k = K_med + length(small))
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
  
  ### For each iteration, relabel new classes using the most common old class assignment
  relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
  for (k in 1:K_med) {
    red_class <- unique_red_classes[k]
    relabel_red_classes[, k] <- 
      apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == red_class]), 1, get_mode)
  }
  
  ### Reduce and reorder parameter estimates using new classes
  pi <- matrix(NA, nrow = M, ncol = K_med)
  theta_global <- array(NA, dim = c(M, J, K_med, R))
  for (m in 1:M) {
    iter_order <- relabel_red_classes[m, ]
    pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
    pi[m, ] <- pi_temp / sum(pi_temp)
    theta_global[m, , , ] <- MCMC_out$theta_global_MCMC[m, , iter_order, ]
  }
  
  ### Return variables
  post_MCMC_out <- list(K_med = K_med, pi = pi, theta_global = theta_global, 
                        dendrogram_global = dendrogram_global)
  return(post_MCMC_out)
  # plot(dendrogram_global, labels = FALSE)
  # rect.hclust(dendrogram_global, k = K_med)
}



get_estimates_WRPC <- function(MCMC_out, post_MCMC_out, n, J, H, R, x_mat) {
  
  ### Identify unique classes using modal exposure categories 
  # Posterior median estimate for theta across iterations
  theta_global_med_temp <- apply(post_MCMC_out$theta_global, c(2, 3, 4), 
                                 function(x) stats::median(x, na.rm = TRUE))
  # Posterior modal exposure categories for each exposure item and reduced class
  theta_global_modes <- apply(theta_global_med_temp, c(1, 2), which.max)
  # Identify unique classes
  unique_global_classes <- which(!duplicated(theta_global_modes, MARGIN = 2))
  # Number of unique classes
  K_red <- length(unique_global_classes)
  
  
  ### Use new classes to adjust and re-normalize posterior samples 
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
  
  ## Get posterior parameter samples for unique classes for theta 
  # Global theta: MxJxKxR
  theta_global_red <- post_MCMC_out$theta_global[, , unique_global_classes, , 
                                                 drop = FALSE]
  theta_global_red <- plyr::aaply(theta_global_red, c(1, 2, 3), function(x) 
    x / sum(x, na.rm = TRUE), .drop = FALSE) # Re-normalize
  # Local theta: MxJxHxR
  theta_local_red <- plyr::aaply(MCMC_out$theta_local_MCMC, c(1, 2, 3), function(x) 
    x / sum(x, na.rm = TRUE), .drop = FALSE) # Re-normalize
  
  
  ### Posterior median estimates 
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
  # Placeholder c_all
  c_all <- MCMC_out$c_all_MCMC[M, ]
  # Global-local assignments to use for obtaining c_all 
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
  
  
  ### Update g using posterior estimates
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
  
  ### Return variables
  estimates <- list(K_red = K_red, pi_red = pi_red, 
                    theta_global_red = theta_global_red,
                    theta_local_red = theta_local_red, nu_red = nu_red, 
                    pi_med = pi_med, theta_global_med = theta_global_med, 
                    theta_local_med = theta_local_med, nu_med = nu_med, 
                    theta_global_modes = theta_global_modes,
                    theta_local_modes = theta_local_modes, c_all = c_all, 
                    g_mat = g_mat, 
                    pred_global_class_probs = pred_global_class_probs)
  return(estimates)
}



