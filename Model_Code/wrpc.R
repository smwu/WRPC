#' Run the WRPC model

library(Rcpp)
library(RcppArmadillo)

wrpc <- function(x_mat, h_all, sampling_wt = NULL, cluster_id = NULL, 
                 stratum_id = NULL, run_sampler = "both", K_max = 30, 
                 adapt_seed = NULL, K_fixed = NULL, fixed_seed = NULL, 
                 class_cutoff_global = 0.05, n_runs = 20000, burn = 10000, 
                 thin = 5, update = 10, switch = 50, save_res = TRUE, 
                 save_path = NULL, alpha_adapt = NULL, eta_global_adapt = NULL, 
                 eta_local_adapt = NULL, a_adapt = NULL, b_adapt = NULL,
                 alpha_fixed = NULL, eta_global_fixed = NULL, 
                 eta_local_fixed = NULL, a_fixed = NULL, b_fixed = NULL) {

  # Begin runtime tracker
  start_time <- Sys.time()
  #================= Read in data ==============================================
  
  print("Read in data")
  # Obtain dimensions
  n <- dim(x_mat)[1]        # Number of individuals
  J <- dim(x_mat)[2]        # Number of exposure items
  R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
               function(x) length(unique(x)))  
  R <- max(R_j)             # Maximum number of exposure categories across items
  H <- length(unique(h_all)) # Number of subpopulations
  
  # If no sampling weights, set all weights to 1 
  if (is.null(sampling_wt)) {
    w_all <- rep(1, n)
  } else {  # Otherwise, obtain normalized weights
    # Obtain normalized weights
    kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
    w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
  }
  
  
  #================= ADAPTIVE SAMPLER ==========================================
  
  ### Run adaptive sampler to obtain number of latent classes, K_fixed
  if (run_sampler %in% c("both", "adapt")) { 
    print("Running adaptive sampler...")
    
    # Set seed
    if(!is.null(adapt_seed)) {
      set.seed(adapt_seed)
    }
    
    ### Initialize hyperparameters 
    # Default hyperparameters for nu
    if (is.null(a_adapt)) {
      a_adapt <- matrix(1, nrow = H, ncol = J)
    }
    if (is.null(b_adapt)) {
      b_adapt <- matrix(1, nrow = H, ncol = J)
    }
    # Default hyperparameters for pi 
    if (is.null(alpha_adapt)) {
      # Hyperparameter for prior for global pi
      alpha_adapt <- rep(1, K_max) / K_max   
    }
    # Default hyperparameters for theta
    if (is.null(eta_global_adapt)) {
      # Hyperparameter for prior for global theta
      # Unviable categories have value 0.01 to prevent rank deficiency issues
      eta_global_adapt <- matrix(0.01, nrow = J, ncol = R) 
      for (j in 1:J) {
        eta_global_adapt[j, 1:R_j[j]] <- rep(1, R_j[j]) 
      }
    }
    if (is.null(eta_local_adapt)) {
      # Hyperparameter for prior for local theta
      # Unviable categories have value 0.01 to prevent rank deficiency issues
      eta_local_adapt <- array(0.01, dim = c(J, H, R)) 
      for (j in 1:J) {
        for (h in 1:H) {
          eta_local_adapt[j, h, 1:R_j[j]] <- rep(1, R_j[j]) 
        }        
      }
    }
    
    ### Initialize RPC model
    # Obtain pi, theta_global, theta_local, nu, c_all, g_mat
    RPC_params <- init_WRPC(alpha = alpha_adapt, a = a_adapt, b = b_adapt, 
                           eta_global = eta_global_adapt, 
                           eta_local = eta_local_adapt, n = n, K = K_max, 
                           H = H, J = J, R = R, h_all = h_all)  
    
    ### Run MCMC adaptive sampler 
    # Obtain nu_MCMC, pi_MCMC, theta_global_MCMC, theta_local_MCMC, c_all_MCMC, 
    # g_mat_MCMC
    MCMC_out <- run_MCMC_Rcpp_WRPC(RPC_params = RPC_params, n = n, K = K_max, 
                                   H = H, J = J, R = R, w_all = w_all, 
                                   h_all = h_all, x_mat = x_mat, g_mat = g_mat,
                                   alpha = alpha_adapt, a = a_adapt, b = b_adapt, 
                                   eta_global = eta_global_adapt, 
                                   eta_local = eta_local_adapt, 
                                   n_runs = n_runs, burn = burn, thin = thin, 
                                   update = update, switch = switch)
    
    ### Post-processing for adaptive sampler 
    # Get median number of global classes with >= cutoff% of individuals, over 
    # all iterations
    M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
    K_MCMC <- rowSums(MCMC_out$pi_MCMC >= class_cutoff_global)
    # Number of classes for the fixed sampler
    K_fixed <- round(stats::median(K_MCMC, na.rm = TRUE))
    print(paste0("K_fixed: ", K_fixed))
    
    ### Create adaptive output list (fixed sampler replaces this if run)
    res <- list(MCMC_out = MCMC_out, K_fixed = K_fixed, K_MCMC = K_MCMC)
    # Save output
    if (save_res) {
      save(res, file = paste0(save_path, "_wrpc_adapt.RData"))
    }
    
    # Reduce memory burden
    rm(RPC_params, MCMC_out)
  }
  
  
  #================= FIXED SAMPLER =============================================
  
  ### Run fixed sampler to obtain posteriors
  if (run_sampler %in% c("both", "fixed")) {
    print("Running fixed sampler...")
    
    # Set seed
    if (!is.null(fixed_seed)) {
      set.seed(fixed_seed)
    }
    
    ### Initialize hyperparameters 
    # Default parameters for beta hyperprior for nu
    if (is.null(a_fixed)) {
      a_fixed <- matrix(1, nrow = H, ncol = J)
    }
    if (is.null(b_fixed)) {
      b_fixed <- matrix(1, nrow = H, ncol = J)
    }
    # Default hyperparameters for pi 
    if (is.null(alpha_fixed)) {
      # Hyperparameter for prior for global pi
      alpha_fixed <- rep(1, K_max) / K_max   
    }
    # Default hyperparameters for theta
    if (is.null(eta_global_fixed)) {
      # Hyperparameter for prior for global theta
      # Unviable categories have value 0.01 to prevent rank deficiency issues
      eta_global_fixed <- matrix(0.01, nrow = J, ncol = R) 
      for (j in 1:J) {
        eta_global_fixed[j, 1:R_j[j]] <- rep(1, R_j[j]) 
      }
    }
    if (is.null(eta_local_fixed)) {
      # Hyperparameter for prior for local theta
      # Unviable categories have value 0.01 to prevent rank deficiency issues
      eta_local_fixed <- array(0.01, dim = c(J, H, R)) 
      for (j in 1:J) {
        for (h in 1:H) {
          eta_local_fixed[j, h, 1:R_j[j]] <- rep(1, R_j[j]) 
        }        
      }
    }
    
    ### Initialize RPC model
    # Obtain pi, theta_global, theta_local, nu, c_all, g_mat
    RPC_params <- init_WRPC(alpha = alpha_fixed, a = a_fixed, b = b_fixed, 
                           eta_global = eta_global_fixed, 
                           eta_local = eta_local_fixed, n = n, K = K_max, 
                           H = H, J = J, R = R, h_all = h_all)  
    ### Run MCMC fixed sampler
    # Obtain nu_MCMC, pi_MCMC, theta_global_MCMC, theta_local_MCMC, c_all_MCMC, 
    # g_mat_MCMC
    MCMC_out <- run_MCMC_Rcpp_WRPC(RPC_params = RPC_params, n = n, K = K_max, 
                                   H = H, J = J, R = R, w_all = w_all, 
                                   h_all = h_all, x_mat = x_mat, g_mat = g_mat,
                                   alpha = alpha_fixed, a = a_fixed, b = b_fixed, 
                                   eta_global = eta_global_fixed, 
                                   eta_local = eta_local_fixed, 
                                   n_runs = n_runs, burn = burn, thin = thin, 
                                   update = update, switch = switch)
    
    # Save output
    if (save_res) {
      save(MCMC_out, file = paste0(save_path, "_wrpc_MCMC_out.RData"))
    }
    ### Post-processing to recalibrate labels and remove extraneous empty classes
    # Obtain K_med, pi, theta_global, theta_local, dendrogram_global
    post_MCMC_out <- post_process_WRPC(MCMC_out = MCMC_out, J = J, R = R, H = H, 
                                       class_cutoff_global = class_cutoff_global)
    
    ### Obtain posterior estimates, reduce number of classes, analyze results
    # Obtain K_red, pi_red, theta_global_red, theta_local_red, nu_red, pi_med
    # theta_global_med, theta_local_med, nu_med, c_all, g_mat, 
    # pred_global_class_probs
    estimates <- get_estimates_WRPC(MCMC_out = MCMC_out, 
                                    post_MCMC_out = post_MCMC_out, n = n, J = J,
                                    H = H, R = R, x_mat = x_mat)
    
    ### Create output list. Replaces adaptive sampler output list
    res <- list(estimates = estimates, MCMC_out = MCMC_out,
                post_MCMC_out = post_MCMC_out, K_fixed = K_fixed)
    
    # Store data variables used
    data_vars <- list(n = n, J = J, R_j = R_j, R = R, H = H, w_all = w_all, 
                      sampling_wt = sampling_wt, x_mat = x_mat,
                      stratum_id = stratum_id, cluster_id = cluster_id)
    res$data_vars <- data_vars
    
    class(res) <- "wrpc"
    
    # Save output
    if (save_res) {
      save(res, file = paste0(save_path, "_wrpc_results.RData"))
    }
  }
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  res$runtime <- runtime
  
  # Return output
  return(res)
}



