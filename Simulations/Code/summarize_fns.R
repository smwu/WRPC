#=================================================
# Functions for summarizing model output
# Author: Stephanie Wu
# Date created: 2024/05/29
# Date updated: 2024/05/29
#=================================================



save_scen_metrics <- function(scenario, samp_i_seq, save_path, wd, data_dir, 
                              WRPC, UNWT, res_dir, subset = FALSE, 
                              dist_type = "mean_abs") {
  
  # Get metrics for models
  metrics_all <- list()
  if (WRPC) {
    print("Getting WRPC results...")
    model <- "wrpc"
    metrics_wrpc <- get_metrics_wrpc(wd = wd, data_dir = data_dir, 
                                         res_dir = res_dir, scenario = scenario, 
                                         model = model, samp_i_seq = samp_i_seq,
                                         subset = subset, dist_type = dist_type,
                                         save_path = save_path)
    metrics_all$metrics_wrpc <- metrics_wrpc
  } 
  if (UNWT) {
    print("Getting unweighted RPC results...")
    model <- "unwt_wrpc"
    metrics_unwt_wrpc <- get_metrics_wrpc(wd = wd, data_dir = data_dir, 
                                        res_dir = res_dir, scenario = scenario, 
                                        model = model, samp_i_seq = samp_i_seq,
                                        subset = subset, dist_type = dist_type,
                                        save_path = save_path)
    metrics_all$metrics_unwt_wrpc <- metrics_unwt_wrpc
  } 

  # Save summary metrics
  save(metrics_all, 
       file = paste0(save_path, "summary.RData"))
  
}


get_metrics_wrpc <- function(wd, data_dir, res_dir, sum_dir, 
                               scenario, model, samp_i_seq, 
                               dist_type = "mean_abs", subset = FALSE,
                               save_path) {
  
  #============== Load data and initialize variables ===========================
  # Load simulated population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop.RData")
  load(pop_data_path)
  
  # Obtain true observed population parameters
  true_params <- get_true_params_wrpc(sim_pop = sim_pop)
  true_K <- as.vector(sim_pop$true_K)
  
  # Initialize variables
  L <- length(samp_i_seq)  # Number of samples
  J <- dim(sim_pop$true_global_thetas)[1]        # Number of exposure items
  runtime_all <- numeric(L)
  # Bias squared using posterior median
  K_all <- K_dist <- pi_dist <- nu_dist <- theta_global_dist <- 
    theta_local_dist <- theta_global_mode_dist <- theta_local_mode_dist <- rep(NA, L) 
  # Posterior variance
  pi_var_all <- nu_var_all <- theta_global_var_all <- theta_local_var_all <- rep(NA, L) 
  # Coverage variables
  pi_cover_all <- matrix(0, nrow=L, ncol=length(true_params$true_pi))
  nu_cover_all <- array(0, dim = c(L, dim(true_params$true_nu)))
  theta_global_cover_all <- array(0, c(L, dim(true_params$true_global_theta)[c(1,2)]))
  theta_local_cover_all <- array(0, c(L, dim(true_params$true_local_theta)[c(1,2)]))
  # MSE for all iterations
  pi_mse_all <- nu_mse_all <- theta_global_mse_all <- theta_local_mse_all <- rep(NA, L)
  
  # Initialize plotting structures
  pi_all <- matrix(NA, nrow=L, ncol=true_K)
  nu_all <- array(NA, dim = c(L, dim(sim_pop$true_nu)))
  theta_global_mode_all <- array(NA, dim=c(L, dim(sim_pop$true_global_patterns)))
  theta_local_mode_all <- array(NA, dim=c(L, dim(sim_pop$true_local_patterns)))
  mode_global_mis_all <- mode_local_mis_all <- rep(NA, L)
  
  
  #============== Get performance metrics for each iteration ===================

  # Get performance metrics for each sample iteration
  for (l in 1:L) { 
    samp_i <- samp_i_seq[l]
    summ_l <- get_metrics_wrpc_i(samp_i = samp_i,  sim_pop = sim_pop, 
                                   wd = wd, data_dir = data_dir, 
                                   res_dir = res_dir, scenario = scenario, 
                                   model = model, dist_type = dist_type, 
                                   subset = subset, save_path = save_path, 
                                   true_params = true_params)
    
    if (!is.null(summ_l)) {
      if (!is.null(summ_l$runtime)) {
        runtime_all[l] <- summ_l$runtime
      }
      K_all[l] <- summ_l$K
      K_dist[l] <- summ_l$K_dist
      pi_dist[l] <- summ_l$pi_dist
      nu_dist[l] <- summ_l$nu_dist
      theta_global_dist[l] <- summ_l$theta_global_dist
      theta_global_mode_dist[l] <- summ_l$theta_global_mode_dist  # for modal probs
      theta_local_dist[l] <- summ_l$theta_local_dist
      theta_local_mode_dist[l] <- summ_l$theta_local_mode_dist  # for modal probs
      
      pi_cover_all[l, ] <- summ_l$pi_cover
      nu_cover_all[l, , ] <- summ_l$nu_cover
      theta_global_cover_all[l, , 1:(dim(summ_l$theta_global_cover)[2])] <- summ_l$theta_global_cover
      theta_local_cover_all[l, , 1:(dim(summ_l$theta_local_cover)[2])] <- summ_l$theta_local_cover
      
      pi_var_all[l] <- summ_l$pi_var
      nu_var_all[l] <- summ_l$nu_var
      theta_global_var_all[l] <- summ_l$theta_global_var
      theta_local_var_all[l] <- summ_l$theta_local_var
      
      pi_mse_all[l] <- summ_l$pi_mse
      nu_mse_all[l] <- summ_l$nu_mse
      theta_global_mse_all[l] <- summ_l$theta_global_mse
      theta_local_mse_all[l] <- summ_l$theta_local_mse
      mode_global_mis_all[l] <- summ_l$mode_global_mis
      mode_local_mis_all[l] <- summ_l$mode_local_mis
      
      # Handle extra estimated classes if necessary
      K_l <- length(summ_l$pi)
      K_storage <- dim(pi_all)[2]
      if (K_storage < K_l) { 
        # Expand estimated matrix and array sizes
        extra <- K_l - true_K
        # Expand pi_all
        filler_pi <- array(NA, dim=c(L, extra))
        pi_all <- abind::abind(pi_all, filler_pi, along = 2)
        # Expand theta_global_mode_all
        filler_theta <- array(NA, dim=c(L, J, extra))
        theta_global_mode_all <- abind::abind(theta_global_mode_all, 
                                              filler_theta, along = 3)
      }
      
      theta_global_mode_all[l, , 1:(dim(summ_l$theta_global_mode)[2])] <- 
        summ_l$theta_global_mode
      theta_local_mode_all[l, , ] <- summ_l$theta_local_mode
      pi_all[l, 1:length(summ_l$pi)] <- summ_l$pi  
      nu_all[l, , ] <- summ_l$nu
    }
  }
  
  
  #============== Calculate bias^2 averaged over sample iterations =============
  K_bias <- mean(K_dist, na.rm = TRUE)
  pi_bias <- mean(pi_dist, na.rm = TRUE)
  nu_bias <- mean(nu_dist, na.rm = TRUE)
  theta_global_bias <- mean(theta_global_dist, na.rm = TRUE)
  theta_local_bias <- mean(theta_local_dist, na.rm = TRUE)
  theta_global_mode_bias <- mean(theta_global_mode_dist, na.rm = TRUE) # for modal probs
  theta_local_mode_bias <- mean(theta_local_mode_dist, na.rm = TRUE) # for modal probs
  
  # Calculated CI width, averaged across iterations
  pi_var <- mean(pi_var_all, na.rm = TRUE)
  nu_var <- mean(nu_var_all, na.rm = TRUE)
  theta_global_var <- mean(theta_global_var_all, na.rm = TRUE)
  theta_local_var <- mean(theta_local_var_all, na.rm = TRUE)
  
  # Calculate class-specific coverage, averaged across iterations
  # Coverage for pi
  pi_cover_avg <- colMeans(pi_cover_all, na.rm = TRUE)
  # Coverage for nu: average over food items
  nu_cover_avg <- rowMeans(colMeans(nu_cover_all, na.rm = TRUE), na.rm = TRUE)
  # Coverage for theta: average over food items
  theta_global_cover_avg <- colMeans(colMeans(theta_global_cover_all, 
                                              na.rm = TRUE), na.rm = TRUE)
  theta_local_cover_avg <- colMeans(colMeans(theta_local_cover_all, 
                                              na.rm = TRUE), na.rm = TRUE)
  
  runtime_avg <- mean(runtime_all, na.rm = TRUE)
  
  #============== Return results ===============================================
  ret_list <- list(K_bias = K_bias, pi_bias = pi_bias, nu_bias = nu_bias,
                   theta_global_bias = theta_global_bias, 
                   theta_local_bias = theta_local_bias,
                   theta_global_mode_bias = theta_global_mode_bias,
                   theta_local_mode_bias = theta_local_mode_bias,
                   pi_var = pi_var, nu_var = nu_var, 
                   theta_global_var = theta_global_var,
                   theta_local_var = theta_local_var,
                   pi_cover_avg = pi_cover_avg, nu_cover_avg = nu_cover_avg,
                   theta_global_cover_avg = theta_global_cover_avg,
                   theta_local_cover_avg = theta_local_cover_avg,
                   runtime_avg = runtime_avg, 
                   # results over all sample iterations
                   K_dist = K_dist, pi_dist = pi_dist, nu_dist = nu_dist,
                   theta_global_dist = theta_global_dist, 
                   theta_global_mode_dist = theta_global_mode_dist,
                   theta_local_mode_dist = theta_local_mode_dist,
                   pi_cover_all = rowMeans(pi_cover_all, na.rm = TRUE), # avg over k
                   nu_cover_all = apply(nu_cover_all, 1, mean), # avg over h,j
                   theta_global_cover_all = apply(theta_global_cover_all, 1, mean), # avg over j,k
                   theta_local_cover_all = apply(theta_local_cover_all, 1, mean), # avg over j,h
                   pi_var_all = pi_var_all, nu_var_all = nu_var_all,
                   theta_global_var_all = theta_global_var_all, 
                   theta_local_var_all = theta_local_var_all,
                   pi_mse_all = pi_mse_all, nu_mse_all = nu_mse_all,
                   theta_global_mse_all = theta_global_mse_all,
                   theta_local_mse_all = theta_local_mse_all)
  
  ret_list[["pi_all"]] <- pi_all
  ret_list[["nu_all"]] <- nu_all
  theta_global_mode <- apply(theta_global_mode_all, c(2,3), function(x) 
    mean(x, na.rm = TRUE))
  mode_global_mis <- mean(mode_global_mis_all, na.rm = TRUE)
  ret_list[["theta_global_mode"]] <- theta_global_mode
  ret_list[["mode_global_mis"]] <- mode_global_mis
  ret_list[["mode_global_mis_all"]] <- mode_global_mis_all
  theta_local_mode <- apply(theta_local_mode_all, c(2,3), function(x) 
    mean(x, na.rm = TRUE))
  mode_local_mis <- mean(mode_local_mis_all, na.rm = TRUE)
  ret_list[["theta_local_mode"]] <- theta_local_mode
  ret_list[["mode_local_mis"]] <- mode_local_mis
  ret_list[["mode_local_mis_all"]] <- mode_local_mis_all
  ret_list[["K_all"]] <- K_all
  
  return(ret_list)
}


# Function to get performance metrics for one iteration
# Inputs:
#   samp_i: sample iteration index
#   model: must be "wrpc" or "unwt_wrpc"
get_metrics_wrpc_i <- function(samp_i, sim_pop, wd, data_dir, res_dir, 
                                 scenario, model, dist_type, subset, save_path,
                                 true_params) {
  
  # Initialize return
  summ_i <- NULL
  
  # Check if summary already exists
  if (file.exists(paste0(save_path, "samp_", samp_i, "_", model, ".RData"))) {
    # Read in summary if it already exists
    print(paste0("Summary for samp ", samp_i, " already exists."))
    load(paste0(save_path, "samp_", samp_i, "_", model, ".RData"))
    
  } else {
    # Obtain true observed population parameters
    # # Need to reload so that filler dimensions do not keep adding over iterations
    true_K <- as.vector(sim_pop$true_K)
    
    # Check that sample data file exists
    sim_data_path <- paste0(wd, data_dir, "scen_", scenario, "/", samp_i, 
                            "sim_samp.RData")
    # Check that results file exists
    sim_res_path <- paste0(wd, res_dir, "scen_", scenario, "/samp_", samp_i, 
                           "_", model, "_results.RData")
    if (!file.exists(sim_data_path)) {
      print(paste0("File does not exist: ", sim_data_path))
    } else if (!file.exists(sim_res_path)) {
      print(paste0("File does not exist: ", sim_res_path))
    } else {
      print(samp_i)
      
      # Read in sample data
      load(sim_data_path)
      
      # Read in results data
      load(sim_res_path)
      if (!is.null(res$runtime)) {
        runtime <- res$runtime
      } else {
        runtime <- NULL
      }
      
      estimates <- res$estimates
      
      M <- dim(estimates$theta_global_red)[1]        # Number of MCMC iterations 
      J <- dim(estimates$theta_global_red)[2]        # Number of exposure items
      R <- dim(estimates$theta_global_red)[4]        # Number of exposure levels
      K <- length(estimates$pi_med)                  # Number of classes
      H <- dim(estimates$theta_local_red)[3]         # Number of subpopulations
      
      # If number of classes is incorrect, fill remaining components with 0's
      if (K > true_K) {
        # If there are extra estimated classes, add 0s to true parameters
        extra <- K - true_K
        true_params$true_pi <- c(true_params$true_pi, rep(0, extra))
        filler <- array(0, dim=c(dim(estimates$theta_global_med)[1], extra, 
                                 dim(estimates$theta_global_med)[3]))
        true_params$true_global_theta <- abind::abind(true_params$true_global_theta, 
                                               filler, along = 2)
      } else if (K < true_K) {
        # If there are missing estimated classes, add 0s to estimated parameters
        missing <- true_K - K
        estimates$pi_med <- c(estimates$pi_med, rep(0, missing))
        filler <- array(0, dim=c(dim(estimates$theta_global_med)[1], missing, 
                                 dim(estimates$theta_global_med)[3]))
        estimates$theta_global_med <- abind::abind(estimates$theta_global_med, 
                                                   filler, along = 2)  
        
        # Add 0's to full MCMC outputs for the missing classes
        estimates$pi_red <- abind::abind(estimates$pi_red, 
                                         array(0, dim=c(M, missing)), along=2)
        estimates$theta_global_red <- abind::abind(estimates$theta_global_red, 
                                            array(0, dim=c(M, J, missing, R)), 
                                            along = 3)
      }
      
      #============== Calculated mean absolute distance (abs bias) ===============
      
      ##### Number of classes, K
      K_dist <- get_dist_wrpc(K, true_K, dist_type = dist_type)
      
      ##### global theta: get dist (Eucl norm) and optimal ordering
      theta_global_perm <- get_theta_dist_wrpc(est_theta = estimates$theta_global_med, 
                                          true_theta = true_params$true_global_theta, 
                                          est_K = K, true_K = true_K, subset = subset,
                                          dist_type = dist_type)
      theta_global_dist <- theta_global_perm$theta_dist
      order <- theta_global_perm$order
      if (subset) {
        order_sub_est <- theta_global_perm$order_sub_est
        order_sub_true <- theta_global_perm$order_sub_true
        K_min <- length(order_sub_est)
      } else {
        order_sub_est <- order
        order_sub_true <- 1:length(order)
        K_min <- true_K
      }
      
      ##### local theta: get dist (Eucl norm). No reordering required
      theta_local_perm <- get_theta_dist_wrpc(est_theta = estimates$theta_local_med, 
                                               true_theta = true_params$true_local_theta, 
                                               est_K = H, true_K = H, subset = subset,
                                               dist_type = dist_type)
      theta_local_dist <- theta_local_perm$theta_dist
      
      ##### pi 
      pi_perm <- get_pi_dist_wrpc(est_pi = estimates$pi_med, 
                                    true_pi = true_params$true_pi, order = order, 
                                    est_K = K, true_K = true_K, subset = subset,
                                    order_sub_est = order_sub_est, 
                                    order_sub_true = order_sub_true,
                                    dist_type = dist_type)
      pi_dist <- pi_perm$pi_dist
      
      ##### nu 
      nu_dist <- get_dist_wrpc(estimates$nu_med, sim_pop$true_nu, 
                                 dist_type = dist_type)
      
      #============== Calculate coverage and CI widths ===========================
      ##### pi
      # Obtain credible intervals for each of the K true clusters
      pi_CI <- apply(estimates$pi_red[, order_sub_est], 2, 
                     function(x) quantile(x, c(0.025, 0.975)))
      # Assign 1 if interval covers true value, 0 if not
      # If a class is missing, defaults to 0 (not covered)
      pi_cover <- numeric(length(order_sub_true))
      pi_cover[order_sub_true] <- ifelse(
        (true_params$true_pi[order_sub_true] >= pi_CI[1,]) & 
          (true_params$true_pi[order_sub_true] <= pi_CI[2,]), 1, 0)
      # Subset to the true number of classes for pi_cover
      pi_cover <- pi_cover[1:true_K]
      # CI width averaged over the components
      pi_var <- mean(apply(pi_CI, 2, diff))
      # MSE
      pi_mse <- mean(apply(estimates$pi_red, 1, function(x) 
        get_dist_wrpc(x[order_sub_est], true_params$true_pi[order_sub_true], 
                        "mean_sq")))
      
      ##### nu
      # Get variance and coverage
      nu_cover <- matrix(NA, nrow = H, ncol = J)
      nu_var_temp <- numeric(H)
      for (h in 1:H) {
        # Obtain credible intervals for each item 
        # Margins of apply are the dimensions that should be preserved
        nu_h_CI <- apply(estimates$nu_red[, h, ], 2, 
                         function(x) quantile(x, c(0.025, 0.975)))
        # Assign 1 if interval covers true value, 0 if not
        nu_cover[h, ] <- ifelse((true_params$true_nu[h, ] >= nu_h_CI[1, ]) &
                               (true_params$true_nu[h, ] <= nu_h_CI[2, ]), 1, 0)
        # CI width measures variation in estimating the modes for each k,
        # averaged over the items
        nu_var_temp[h] <- mean(apply(nu_h_CI, 2, diff))
      }
      # CI width averaged over the subpopulations
      nu_var <- mean(nu_var_temp)
      # MSE
      nu_mse <- mean(apply(estimates$nu_red, 1, function(x) 
        get_dist_wrpc(x, true_params$true_nu, "mean_sq")))
      
      
      ##### global theta
      # Theta mode consumption levels for each item and class (pxK)
      est_global_modes <- apply(estimates$theta_global_med[, order_sub_est, ], c(1,2), which.max)
      true_global_modes <- apply(true_params$true_global_theta[, order_sub_true, ], c(1,2), 
                          which.max)
      # True modal probabilities for each item and class (pxK_true)
      true_theta_global_modal <- apply(true_params$true_global_theta[ , order_sub_true, ], 
                                c(1,2), max) 
      # Estimated modal probabilities for each item and class (pxK)
      est_theta_global_modal <- apply(estimates$theta_global_med[, order_sub_est, ], 
                               c(1, 2), max)
      # Get distance for modal probabilities
      theta_global_mode_dist <- get_dist_wrpc(est_theta_global_modal, 
                                              true_theta_global_modal, 
                                              dist_type = dist_type) 
      # Get variance and coverage
      # Initialize theta_cover  
      theta_global_var_temp <- numeric(K_min)
      theta_global_cover <- array(NA, dim = dim(true_params$true_global_theta)[c(1,2)])
      for (k in 1:K_min) {
        # Subset theta for cluster k
        est_theta_global_k <- estimates$theta_global_red[, , order_sub_est[k], ]
        # Each row provides the indices for one row of modal probabilities
        modal_idx <- cbind(rep(1:M, each = J), rep(1:J, times = M), 
                           rep(est_global_modes[, k], times = M))
        # estimated probabilities for the mode for cluster k (Mxp)
        est_theta_global_k_modal <- matrix(est_theta_global_k[modal_idx], 
                                           ncol = J, byrow = TRUE)
        # Obtain credible intervals for each item 
        # Margins of apply are the dimensions that should be preserved
        theta_global_CI <- apply(est_theta_global_k_modal, 2, 
                          function(x) quantile(x, c(0.025, 0.975)))
        theta_global_cover[, order_sub_true[k]] <- ifelse(
          (true_theta_global_modal[, k] >= theta_global_CI[1, ]) &
            (true_theta_global_modal[, k] <= theta_global_CI[2, ]), 1, 0)
        # CI width measures variation in estimating the modes for each k,
        # averaged over the items
        theta_global_var_temp[k] <- mean(apply(theta_global_CI, 2, diff))
      }
      # Subset to the true number of classes for theta_cover
      theta_global_cover <- theta_global_cover[, 1:true_K]
      # CI width averaged over the classes
      theta_global_var <- mean(theta_global_var_temp)
      # MSE
      theta_global_mse <- mean(apply(estimates$theta_global_red, 1, function(x) 
        get_dist_wrpc(x[, order_sub_est, ], 
                        true_params$true_global_theta[, order_sub_true, ], "mean_sq")))
      
      
      ##### local theta
      # Theta mode consumption levels for each item and class (pxK)
      est_local_modes <- apply(estimates$theta_local_med, c(1,2), which.max)
      true_local_modes <- apply(true_params$true_local_theta, c(1,2), which.max)
      # True modal probabilities for each item and class (pxK_true)
      true_theta_local_modal <- apply(true_params$true_local_theta, c(1,2), max) 
      # Estimated modal probabilities for each item and class (pxK)
      est_theta_local_modal <- apply(estimates$theta_local_med, c(1, 2), max)
      # Get distance for modal probabilities
      theta_local_mode_dist <- get_dist_wrpc(est_theta_local_modal, 
                                             true_theta_local_modal, 
                                             dist_type = dist_type) 
      # Get variance and coverage
      # Initialize theta_cover  
      theta_local_var_temp <- numeric(H)
      theta_local_cover <- array(NA, dim = dim(true_params$true_local_theta)[c(1,2)])
      for (h in 1:H) {
        # Subset theta for subpop h
        est_theta_local_h <- estimates$theta_local_red[, , h, ]
        # Each row provides the indices for one row of modal probabilities
        modal_idx <- cbind(rep(1:M, each = J), rep(1:J, times = M), 
                           rep(est_local_modes[, h], times = M))
        # estimated probabilities for the mode for subpop h (MxJ)
        est_theta_local_h_modal <- matrix(est_theta_local_h[modal_idx], 
                                           ncol = J, byrow = TRUE)
        # Obtain credible intervals for each item 
        # Margins of apply are the dimensions that should be preserved
        theta_local_CI <- apply(est_theta_local_h_modal, 2, 
                                 function(x) quantile(x, c(0.025, 0.975)))
        theta_local_cover[, h] <- ifelse(
          (true_theta_local_modal[, h] >= theta_local_CI[1, ]) &
            (true_theta_local_modal[, h] <= theta_local_CI[2, ]), 1, 0)
        # CI width measures variation in estimating the modes for each k,
        # averaged over the items
        theta_local_var_temp[h] <- mean(apply(theta_local_CI, 2, diff))
      }
      # CI width averaged over the classes
      theta_local_var <- mean(theta_local_var_temp)
      # MSE
      theta_local_mse <- mean(apply(estimates$theta_local_red, 1, function(x) 
        get_dist_wrpc(x, true_params$true_local_theta, "mean_sq")))
      
      
      #============== Parameter estimate plot outputs ============================
      ##### global theta
      # Global theta mode consumption levels for each item and class (pxK)
      est_global_modes <- apply(estimates$theta_global_med[, order, ], c(1,2), 
                                which.max)
      # Get theta mode
      theta_global_mode <- array(NA, dim = c(J, max(K, true_K)))
      theta_global_mode[, 1:length(order)] <- est_global_modes
      # Mode mismatches
      mode_global_mis <- sum(abs(est_global_modes[, 1:true_K] - 
                                   sim_pop$true_global_patterns))
      
      ##### local theta
      # Local theta mode consumption levels for each item and class (pxK)
      theta_local_mode <- apply(estimates$theta_local_med, c(1,2), which.max)
      # Mode mismatches
      mode_local_mis <- sum(abs(theta_local_mode - sim_pop$true_local_patterns))
      
      ##### pi
      pi <- numeric(length(order))
      pi[1:length(order)] <- estimates$pi_med[order]
      
      ##### nu
      nu <- estimates$nu_med
      
      # Return performance metrics for the iteration
      summ_i <- list(runtime = runtime, 
                     # bias using dist_type
                     K_dist = K_dist, pi_dist = pi_dist, 
                     nu_dist = nu_dist, theta_global_dist = theta_global_dist, 
                     theta_local_dist = theta_local_dist, 
                     # coverage
                     pi_cover = pi_cover, nu_cover = nu_cover, 
                     theta_global_cover = theta_global_cover, 
                     theta_local_cover = theta_local_cover,
                     # CI width
                     pi_var = pi_var, nu_var = nu_var,
                     theta_global_var = theta_global_var, 
                     theta_local_var = theta_local_var,
                     # MSE
                     pi_mse = pi_mse, nu_mse = nu_mse, 
                     theta_global_mse = theta_global_mse, 
                     theta_local_mse = theta_local_mse,
                     # Components for plotting 
                     theta_global_mode = theta_global_mode, 
                     theta_global_mode_dist = theta_global_mode_dist, 
                     theta_local_mode = theta_local_mode,
                     theta_local_mode_dist = theta_local_mode_dist,
                     mode_global_mis = mode_global_mis,
                     mode_local_mis = mode_local_mis,
                     pi = pi, nu = nu, K = K)
      
      # Save summary metrics
      save(summ_i, file = paste0(save_path, "samp_", samp_i, "_", model, ".RData"))
    }
  }
  
  # Return results
  return(summ_i)
}


get_true_params_wrpc <- function(sim_pop) {
  # Get true pi using population data
  true_pi <- tabulate(sim_pop$true_Ci) / length(sim_pop$true_Ci)
  
  # Get true nu
  true_nu_dim <- dim(sim_pop$true_nu)
  true_nu <- matrix(NA, nrow = true_nu_dim[1], ncol = true_nu_dim[2])
  for (h in 1:true_nu_dim[1]) {
    for (j in 1:true_nu_dim[2]) {
      true_nu[h, j] <- sum((sim_pop$true_Gij[, j] == 1) & (sim_pop$true_Hi == h)) /
        sum(sim_pop$true_Hi == h)
    }
  }
  
  # Get true global theta
  theta_global_dim <- dim(sim_pop$true_global_thetas)
  true_global_theta <- array(NA, dim=theta_global_dim)
  for (j in 1:theta_global_dim[1]) {
    for (k in 1:theta_global_dim[2]) {
      for (r in 1:theta_global_dim[3]) {
        true_global_theta[j,k,r] <- 
          sum((sim_pop$X_data[,j]==r) & (sim_pop$true_Ci==k) & (sim_pop$true_Gij[,j]==1)) / 
          sum((sim_pop$true_Ci==k) & (sim_pop$true_Gij[,j]==1))
      }
    }
  }
  # Get true local theta
  theta_local_dim <- dim(sim_pop$true_local_thetas)
  true_local_theta <- array(NA, dim=theta_local_dim)
  for (j in 1:theta_local_dim[1]) {
    for (h in 1:theta_local_dim[2]) {
      for (r in 1:theta_local_dim[3]) {
        true_local_theta[j,h,r] <- 
          sum((sim_pop$X_data[,j]==r) & (sim_pop$true_Hi==h) & (sim_pop$true_Gij[,j]==0)) / 
          sum((sim_pop$true_Hi==h) & (sim_pop$true_Gij[,j]==0))
      }
    }
  }
  
  return(list(true_pi = true_pi, true_global_theta = true_global_theta,
              true_local_theta = true_local_theta,
              true_nu = true_nu))
}


# Calculate distance between par1 and par2 using a specified distance metric
get_dist_wrpc <- function(par1, par2, dist_type = "mean_abs") {
  if (dist_type == "mean_abs") {  # Mean absolute error
    dist <- mean(abs(par1 - par2))
  } else if (dist_type == "sum_sq") {  # Frobenius norm / squared Euclidean norm
    dist <- sum((par1 - par2)^2)
  } else if (dist_type == "mean_sq") {  # MSE
    dist <- mean((par1 - par2)^2)
  } else {
    stop("Error: dist_type must be 'mean_abs', 'sum_sq', or 'mean_sq' ")
  }
  return(dist)
}


get_theta_dist_wrpc <- function(est_theta, true_theta, est_K, true_K, 
                                  subset, dist_type) {
  
  ### First get minimum distance using full vectors (with additional 0's)
  ### to get optimal ordering
  # Find all permutations of est_theta and true_theta with filler 0's
  all_perms <- gtools::permutations(n = dim(est_theta)[2], r = dim(true_theta)[2])
  # Obtain vector of mean absolute distance between est and true theta, 
  # calculated for each permutation
  dist_all_perms <- numeric(nrow(all_perms))
  for (i in 1:nrow(all_perms)) {
    est_theta_perm <- est_theta[,all_perms[i, ],]
    dist_all_perms[i] <- get_dist_wrpc(est_theta_perm, true_theta, 
                                         dist_type = dist_type) 
  }
  # Obtain optimal ordering of classes
  order <- all_perms[which.min(dist_all_perms), ]
  est_theta_perm <- est_theta[ , order, ]
  
  # Initialize subset ordering of classes
  order_sub_est <- order
  order_sub_true <- 1:true_K
  # Lowest dist out of all permutations
  theta_dist <- min(dist_all_perms)
  
  ### Option to use this minimum distance, or calculate minimum distance after
  ### subsetting (subset == TRUE)
  if (subset) {   # Calculate distance after subsetting
    if (est_K < true_K) {  # If missing a true class
      theta_sub <- get_subset_dist_wrpc(large_par = true_theta[, 1:true_K, ], 
                                          small_par = est_theta[, 1:est_K, ], 
                                          param_name = "theta", dist_type = dist_type)
      theta_dist <- theta_sub$par_dist         # Distance
      order_sub_true <- theta_sub$order_large  # Subset order for true_theta
      order_sub_est <- 1:est_K                 
    } else if (true_K < est_K) {  # If extra class
      theta_sub <- get_subset_dist_wrpc(large_par = est_theta[, 1:est_K, ], 
                                          small_par = true_theta[, 1:true_K, ],
                                          param_name = "theta", dist_type = dist_type)
      theta_dist <- theta_sub$par_dist
      order_sub_est <- theta_sub$order_large  # Subset order for est_theta
      order_sub_true <- 1:true_K
    } else {  # true_K == est_K
      # Lowest dist out of all permutations
      theta_dist <- min(dist_all_perms)
      order_sub_est <- order
      order_sub_true <- 1:true_K
    }
  }
  
  ### Return dist, ordering, reordered estimate, and subsetted orderings
  return(list("theta_dist" = theta_dist, "order" = order, 
              "est_theta_perm" = est_theta_perm, 
              "order_sub_est" = order_sub_est,
              "order_sub_true" = order_sub_true))
}



get_pi_dist_wrpc <- function(est_pi, true_pi, order, est_K, true_K, subset = TRUE,
                             order_sub_est, order_sub_true, dist_type) {
  ### Use input optimal ordering
  ### Option to subset when calculating minimum distance (subset == TRUE)
  if (!subset) {  # No subsetting
    pi_dist <- get_dist_wrpc(est_pi[order], true_pi, dist_type = dist_type)
  } else {  # Calculate distance after subsetting
    pi_dist <- get_dist_wrpc(est_pi[order_sub_est], true_pi[order_sub_true], 
                             dist_type = dist_type)
  }
  
  ### Return dist, ordering, and reordered estimate
  return(list("pi_dist" = pi_dist, "order" = order, "est_pi_perm" = est_pi[order]))
}


# Get distance between `large_par` and `small_par`, subsetting to the dimensions 
# of `small_par`. Used in situation when theta distance is to be calculated 
# after subsetting to the smaller number of classes in the case of a class 
# mismatch
get_subset_dist_wrpc <- function(large_par, small_par, param_name, dist_type) {
  if (param_name == "theta") {
    large_K <- dim(large_par)[2]   ## Change to handle 0's
    sum(large_par[1, , 1] != 0)
    small_K <- dim(small_par)[2]
  } else if (param_name == "pi") {
    large_K <- length(large_par)
    small_K <- length(small_par)
  } else {
    stop("Error: 'param_name' must be either 'theta' or 'pi'")
  }
  # Find all subsets of large_K with size equal to small_K
  sub_perms <- gtools::permutations(n = large_K, r = small_K)
  # Obtain dist (Frobenius norm) between large_par and small_par per permutation
  dist_sub_perms <- numeric(nrow(sub_perms))
  for (i in 1:nrow(sub_perms)) {
    if (param_name == "theta") {
      large_par_sub <- large_par[ , sub_perms[i, ], ]
    } else if (param_name == "pi") {
      large_par_sub <- large_par[sub_perms[i, ]]
    } 
    dist_sub_perms[i] <- get_dist_wrpc(small_par, large_par_sub, "mean_abs")
  }
  # Lowest dist out of all permutations
  par_dist <- min(dist_sub_perms)
  # Ordering corresponding to lowest dist
  order_large <- sub_perms[which.min(dist_sub_perms), ]
  
  return(list(par_dist = par_dist, order_large = order_large))
}



#==================== Tables ===================================================
create_app_tables_wrpc <- function(save_paths, scenarios, scen_names, 
                                     overall_name, format = "latex", 
                                     digits = 3, WRPC = TRUE, UNWT = TRUE,
                                     modal = FALSE) {
  num_scen <- length(scenarios)
  # models depending on WRPC and UNWT true/false
  model <- list()
  if (UNWT) {
    model <- c(model, "Unweighted RPC")
  }
  if (WRPC) {
    model <- c(model, "WRPC")
  }
  model <- unlist(model)
  # multiplier for number of rows depending on models
  mult <- length(model)
  
  metrics_wrpc_df <- as.data.frame(matrix(NA, nrow = mult*length(scenarios), 
                                            ncol = 15))
  colnames(metrics_wrpc_df) <- c(overall_name, "Model", 
                                   "$K$ |Bias|", "$\\pi$ |Bias|",  
                                   "$\\theta_G$ |Bias|", 
                                   "$\\theta_L$ |Bias|", "$\\nu$ |Bias|",
                                   "$\\pi$ Var", "$\\theta_G$ Var", 
                                   "$\\theta_L$ Var", "$\\nu$ Var",
                                   "$\\pi$ Cov","$\\theta_G$ Cov",
                                   "$\\theta_L$ Cov", "$\\nu$ Cov")
  metrics_wrpc_df[, 1] <- rep(scen_names, each = mult)
  metrics_wrpc_df[, 2] <- rep(model, num_scen)  
  # output_inds <- 1:7
  if (modal) {  # use modal theta for bias 
    output_inds <- c(1, 2, 6, 7, 3,  # bias
                     8, 10, 11, 9) # ci-width
  } else {  # use all theta for bias
    output_inds <- c(1, 2, 4, 5, 3,  # bias
                     8, 10, 11, 9) # ci-width
  }
  row_ind <- 1
  for (i in 1:num_scen) {
    save_path <- save_paths[i]
    load(paste0(save_path, "summary.RData"))
    if (UNWT) {
      metrics_wrpc_df[row_ind, -c(1,2)] <- 
        c(metrics_all$metrics_unwt_wrpc[output_inds], 
          mean(metrics_all$metrics_unwt_wrpc$pi_cover_avg), 
          mean(metrics_all$metrics_unwt_wrpc$theta_global_cover_avg), 
          mean(metrics_all$metrics_unwt_wrpc$theta_local_cover_avg),
          mean(metrics_all$metrics_unwt_wrpc$nu_cover_avg))
      row_ind <- row_ind + 1
    }
    if (WRPC) {
      metrics_wrpc_df[row_ind, -c(1,2)] <- 
        c(metrics_all$metrics_wrpc[output_inds], 
          mean(metrics_all$metrics_wrpc$pi_cover_avg), 
          mean(metrics_all$metrics_wrpc$theta_global_cover_avg), 
          mean(metrics_all$metrics_wrpc$theta_local_cover_avg),
          mean(metrics_all$metrics_wrpc$nu_cover_avg))
      row_ind <- row_ind + 1
    }
  }
  
  metrics_wrpc_df %>% 
    kbl(digits = digits, align = "rrrrrrrrrrrr", booktabs = TRUE, format = format,
        caption = "Summary of mean absolute bias, 95% credible interval width, and coverage for simulations based on posterior samples.") %>%
    kable_classic() %>%
    kable_styling(full_width = FALSE)
}


#========================= Plot theta ==========================================

plot_theta_patterns_wolcan <- function(wd, data_dir, scenario, save_path) {
  # Load summary
  load(paste0(save_path, "summary.RData"))
  
  # Load true population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData")
  load(pop_data_path)
  # Obtain true observed population parameters
  true_params <- get_true_params_wolcan(sim_pop = sim_pop)
  
  p_true <- theta_mode_plot(sim_pop$true_global_patterns, "True Classes") + 
    guides(fill = guide_legend(reverse = FALSE)) +
    labs(fill = "Modal θ Level")
  p_unwt <- theta_mode_plot(metrics_all$metrics_unwt_wrpc$theta_mode, "Unweighted Classes")
  p_WOLCAN <- theta_mode_plot(metrics_all$metrics_wrpc$theta_mode, "WOLCAN Classes")
  p_comb <- ggarrange(p_true, 
                      p_unwt + theme(axis.title.y = element_blank()), 
                      p_WOLCAN + theme(axis.title.y = element_blank()), 
                      nrow = 1, common.legend = TRUE, legend = "top")
  return(p_comb)
}

theta_mode_plot <- function(theta_plot_data, x_label) {
  p <- dim(theta_plot_data)[1]
  K <- dim(theta_plot_data)[2]
  Item <- factor(as.character(1:p), levels = as.character(p:1))
  theta_plot <- data.frame(theta_plot_data, Item)
  colnames(theta_plot) <- c(1:K, "Item")
  theta_plot <- theta_plot %>% gather("Class", "Level", 1:K) 
  patterns <- ggplot(theta_plot, aes(x=Class, y=Item, fill=Level)) + 
    theme_classic() +
    xlab(x_label) +
    geom_tile(color="gray") + 
    geom_text(aes(label = round(Level,2)), col="white", cex=2.5) +
    scale_fill_gradient(trans = "reverse")
  return(patterns)
}

theta_mode_plot_wolcan <- function(theta_plot_data, x_label) {
  p <- dim(theta_plot_data)[1]
  K <- dim(theta_plot_data)[2]
  Item <- factor(as.character(1:p), levels = as.character(p:1))
  theta_plot <- data.frame(theta_plot_data, Item)
  colnames(theta_plot) <- c(1:K, "Item")
  theta_plot <- theta_plot %>% gather("Class", "Level", 1:K) 
  patterns <- ggplot(theta_plot, aes(x=Class, y=Item, fill=as.factor(Level))) + 
    theme_classic() +
    xlab(x_label) +
    geom_tile(color="black", linewidth = 0.1) + 
    # geom_text(aes(label = round(Level,2)), col="white", cex=2.5) +
    scale_fill_brewer(type = "seq", palette = "RdYlBu", direction = -1,
                      name = "Level")
  # scale_fill_gradient(trans = "reverse")
  return(patterns)
}

#===================== Plot pi =================================================

plot_pi_patterns_wolcan <- function(wd, data_dir, scenario, samp_i_seq, 
                                    save_path, y_lim = c(0,1)) {
  # Load true population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData")
  load(pop_data_path)
  # Obtain true observed population parameters
  true_params <- get_true_params_wolcan(sim_pop = sim_pop)
  true_pi <- true_params$true_pi
  true_K <- as.vector(sim_pop$K)
  samp_pi <- get_avg_over_samps_wolcan(wd = wd, data_dir = data_dir, 
                                       scenario = scenario, 
                                       samp_i_seq = samp_i_seq)$avg_samp_pi
  
  # Load summary
  load(paste0(save_path, "summary.RData"))
  
  # Load simulated sample data
  L <- length(metrics_all$metrics_wrpc$K_dist)
  
  pi_plot_data <- as.data.frame(rbind(metrics_all$metrics_unwt_wrpc$pi_all[, 1:true_K],
                                      metrics_all$metrics_wrpc$pi_all[, 1:true_K]))
  colnames(pi_plot_data) <- 1:ncol(pi_plot_data)
  pi_plot_data$Model <- c(rep("Unweighted", times=L),
                          rep("WOLCAN", times=L))
  pi_plot_data <- pi_plot_data %>% gather("pi_component", "value", -Model)
  ggplot(pi_plot_data, aes(x=pi_component, y=value, fill=Model)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() +
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_params$true_pi[1],
                             yend=true_params$true_pi[1]),color="red") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_params$true_pi[2],
                             yend=true_params$true_pi[2]),color="red") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_params$true_pi[3],
                             yend=true_params$true_pi[3]),color="red") +
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=samp_pi[1], yend=samp_pi[1]),
                 color="red", linetype = "dashed") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=samp_pi[2], yend=samp_pi[2]),
                 color="red", linetype = "dashed") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=samp_pi[3], yend=samp_pi[3]),
                 color="red", linetype = "dashed") +
    # ggtitle(paste0("Parameter estimation for π across samples")) + 
    xlab("Latent Class") + ylab("π Value") + 
    theme(legend.position="top") 
}


get_avg_over_samps_wolcan <- function(wd, data_dir, scenario, samp_i_seq) {
  # Load true population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData")
  load(pop_data_path)
  # Obtain true observed population parameters
  true_params <- get_true_params_wolcan(sim_pop = sim_pop)
  
  L <- length(samp_i_seq)
  samp_pis <- array(NA, dim = c(length(samp_i_seq), length(sim_data$true_pi)))
  for (l in 1:L) { # For each sample iteration
    samp_i = samp_i_seq[l]
    sim_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", 
                            samp_i, "_B_wolcan.RData")
    load(sim_data_path)
    samp_pis[l, ] <- tabulate(sim_samp_B$c_all) / length(sim_samp_B$c_all)
  }
  avg_samp_pi <- apply(samp_pis, 2, function(x) mean(x, na.rm = TRUE))
  
  return(list(avg_samp_pi = avg_samp_pi))
}


# # Compare distributions
# hist(sim_samp_B$true_pi_B, breaks = 30, col = rgb(0,0,1,1/4), ylim = c(0, 100),
#      xlab = "true_pi_B (blue) vs. hat_pi_B (red)")
# hist(hat_pi_B, breaks = 30, col = rgb(1,0,0,1/4), add = TRUE)
# hist(hat_pi_z, breaks = 20)
# mean(abs(hat_pi_B - sim_samp_B$true_pi_B))
# 
# hist(hat_pi_B_dist, breaks = 30)
# hist(apply(hat_pi_B_dist, 1, max), breaks = 30)
# hist(apply(hat_pi_B_dist, 1, function(x) sum(x > 1)), breaks = 30)
# 
# # Compare distributions
# hist(sim_samp_B$true_pi_B, breaks = 30, col = rgb(0,0,1,1/4), ylim = c(0, 100),
#      xlab = "true_pi_B (blue) vs. hat_pi_B (red)")
# hist(hat_pi_B, breaks = 30, col = rgb(1,0,0,1/4), add = TRUE)
# 
# 
# 
# get_summ_stats <- function(res, res_true = NULL) {
#   if (!is.null(res$estimates_adjust)) {
#     est <- res$estimates_adjust
#   } else {
#     est <- res$estimates
#   }
#   # This is a rough fix right now. Should be updated to proper version using theta
#   est_pi <- sort(est$pi_med)
#   #true_pi <- sort(res_true$estimates$pi_med)
#   true_pi <- sort(prop.table(table(sim_pop$c_all)))
#   ci_est <- apply(est$pi_red, 2, function(x) quantile(x, c(0.025, 0.975)))
#   ci_est <- t(apply(ci_est, 1, sort))
#   
#   abs_bias <- mean(abs(est_pi - true_pi))
#   ci_width <- mean(ci_est[2, ] - ci_est[1, ])
#   cover <- mean(sapply(1:ncol(ci_est), function(x) 
#     ci_est[1, x] <= true_pi[x] & true_pi[x] <= ci_est[2, x]))
#   
#   return(list(abs_bias = abs_bias, ci_width = ci_width, cover = cover))
# }