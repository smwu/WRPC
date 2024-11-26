#' Run the WRPC model

library(Rcpp)
library(RcppArmadillo)

#' Run the WRPC model
#'
#' @description
#' `wrpc` runs a weighted robust profile clustering model, then saves and 
#' returns the results.
#'
#' @param x_mat Matrix of multivariate categorical exposures. nxJ
#' @param h_all Vector of subpopulation indicators. nx1. 
#' @param sampling_wt Vector of survey sampling weights. nx1. Default is `NULL`, 
#' indicating no sampling weights and setting all weights to 1. 
#' @param cluster_id Vector of individual cluster IDs. nx1. Default is `NULL`,
#' indicating each individual is their own cluster.
#' @param stratum_id Vector of individual stratum IDs. nx1. Default is `NULL`,
#' indicating no stratification.
#' @param run_sampler String specifying which sampler(s) should be run. Must be 
#' one of `"both"` (default), `"fixed"`, or `"adapt"`. See Details.
#' @param K_max Upper limit for number of classes. Default is 30.
#' @param adapt_seed Numeric seed for adaptive sampler. Default is `NULL`.
#' @param K_fixed True number of classes, if known. Default is `NULL`, as it is
#' not necessary for the adaptive sampler. If bypassing the adaptive sampler and
#' running the fixed sampler directly, need to specify a value here. See Details.
#' @param fixed_seed Numeric seed for fixed sampler. Default is `NULL`.
#' @param class_cutoff_global Minimum class size proportion when determining 
#' number of global classes. Default is 0.05.
#' @param n_runs Number of MCMC iterations. Default is 20000.
#' @param burn Number of MCMC iterations to drop as a burn-in period. Default is 10000.
#' @param thin Thinning factor for MCMC iterations. Default is 5.
#' @param switch Number specifying how often label switching via a random 
#' permutation sampler should be applied to encourage mixing. Default is 50.
#' @param update Number specifying that MCMC progress updates should be printed 
#' every `update` iterations. Default is 10.
#' @param save_res Boolean specifying if results should be saved. Default = `TRUE`.
#' @param save_path String specifying directory and file name to save results, 
#' e.g., "~/Documents/run". Default is `NULL`.
#' @param alpha_adapt Adaptive sampler hyperparameter for prior for global class 
#' membership probabilities \eqn{\pi}. Default is `NULL` and default values are 
#' used (see Details). If specified, must be (`K_max`)x1 vector. 
#' @param eta_global_adapt Adaptive sampler hyperparameter for prior for global item 
#' level probabilities \eqn{\theta_{G, jkr}} for each item \eqn{j}, global class 
#' \eqn{k}, and category \eqn{r}, assumed to be the same across classes. 
#' Default is `NULL` and default values are used (see Details). 
#' If specified, must be JxR matrix, where J is the number of exposure items and 
#' R is the maximum number of categories for the exposure.
#' @param eta_local_adapt Adaptive sampler hyperparameter for prior for local item 
#' level probabilities \eqn{\theta_{L, jh\cdot}} for each item \eqn{j}, 
#' subpopulation \eqn{h}, and category \eqn{r}. Default is `NULL` and default 
#' values are used (see Details). If specified, must be JxHxR array, where J is 
#' the number of exposure items, H is the number of subpopulations, and R is the 
#' maximum number of categories for the exposure.
#' @param a_adapt Adaptive sampler shape hyperparameter for \eqn{\nu}. Default 
#' is `NULL` and default values are used (see Details). If specified, must be HxJ matrix.
#' @param b_adapt Adaptive sampler rate hyperparameter for \eqn{\nu}. Default 
#' is `NULL` and default values are used (see Details). If specified, must be HxJ matrix.
#' @param alpha_fixed Fixed sampler hyperparameter for prior for \eqn{\pi}. Default is 
#' `NULL` and default values are used (see Details). If specified, must be (`K_fixed`)x1 vector. 
#' @param eta_global_fixed Fixed sampler hyperparameter for prior for global item 
#' level probabilities \eqn{\theta_{G, jkr}} for each item \eqn{j}, global class 
#' \eqn{k}, and category \eqn{r}, assumed to be the same across classes. 
#' Default is `NULL` and default values are used (see Details). If specified, 
#' must be JxR matrix.
#' @param eta_local_fixed Fixed sampler hyperparameter for prior for local item 
#' level probabilities \eqn{\theta_{L, jh\cdot}} for each item \eqn{j}, 
#' subpopulation \eqn{h}, and category \eqn{r}. Default is `NULL` and default 
#' values are used (see Details). If specified, must be JxHxR array.
#' @param a_fixed Fixed sampler shape hyperparameter for \eqn{\nu}. Default 
#' is `NULL` and default values are used (see Details). If specified, must be HxJ matrix.
#' @param b_fixed Fixed sampler rate hyperparameter for \eqn{\nu}. Default 
#' is `NULL` and default values are used (see Details). If specified, must be HxJ matrix.
#' 
#' @details 
#' If no survey sample adjustments are desired, leave `sampling_wt`, `stratum_id`, 
#' and `cluster_id` to their default `NULL` values.
#' 
#' By default, the function will run two samplers: the adaptive sampler followed 
#' by the fixed sampler. The adaptive sampler determines the number of latent global 
#' classes, which is then used in the fixed sampler for parameter estimation. 
#' If the number of latent global classes is already known and only the fixed sampler 
#' needs to be run, specify `"fixed"` for the `run_sampler` argument and specify a 
#' number for `K_fixed`. If only the adaptive sampler is to be run, specify 
#' `"adapt"` for the `run_sampler` argument. Use `adapt_seed` (default is `NULL`) 
#' to specify a seed for the adaptive sampler, and use `fixed_seed` (default is 
#' `NULL`) to specify a separate seed for the fixed sampler. 
#' 
#' `x_mat` is an nxJ matrix with each row corresponding to the J-dimensional 
#' categorical exposure for an individual. `K_max` is the maximum number of 
#' global latent classes allowable, to be used if the adaptive sampler is run. 
#' `class_cutoff_global` is the minimum size of each global class as a proportion 
#' of the population, used when determining the number of latent classes.  
#' 
#' To save results, set `save_res = TRUE` (default) and `save_path` to a string
#' that specifies both the location and the beginning of the file name 
#' (e.g., "~/Documents/run"). The file name will have "_wrpc_adapt.RData" or 
#' "_wrpc_results.RData" appended to it.
#'
#' If hyperparameters for the adaptive or fixed sampler are left as `NULL` 
#' (default), the following default values are used. Let \eqn{K} refer to 
#' `K_max` for the adaptive sampler and `K_fixed` for the fixed sampler. 
#' For \eqn{\pi}, a Dirichlet prior with hyperparameter \eqn{\alpha = 1/K} for 
#' each component. For \eqn{\theta_{G, jk\cdot}} and \eqn{\theta_{L, jh\cdot}}, 
#' a Dirichlet prior with hyperparameter \eqn{\eta_j} equal to `rep(1, R_j)` 
#' where `R_j` is the number of categories for exposure item j. If `R_j < R`, 
#' the remaining categories have hyperparameter set to 0.01. This is done 
#' independently for each exposure item j and is assumed to be the same across 
#' latent classes. For \eqn{\nu_{hj}}, a Beta prior with shape \eqn{a = 1} and 
#' rate \eqn{b = 1}. Note that hyperparameters for the fixed sampler should 
#' probably only be specified if running the fixed sampler directly, bypassing 
#' the adaptive sampler. 
#'
#' @return
#' If the fixed sampler is run, returns an object `res` of class `"wrpc"`; a 
#' list containing the following:
#' \describe{
#'   \item{\code{estimates}}{List of posterior model results, resulting from a 
#'   call to [get_estimates_wrpc()]}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used, including:
#'   `n`: Sample size.
#'   `J`: Number of exposure items.
#'   `R_j`: Number vector of number of exposure categories for each item; Jx1.
#'   `R`: Maximum number of exposure categories across items.
#'   `H`: Number of subpopulations.
#'   `w_all`: Vector of sampling weights normalized to sum to n; nx1.
#'   `sampling_wt`: Vector of survey sampling weights; nx1.
#'   `x_mat`: Matrix of multivariate categorical exposures; nxJ.
#'   `h_all`: Vector of subpopulation indicators; nx1.
#'   `stratum_id`: Vector of individual stratum IDs; nx1 or NULL.
#'   `cluster_id`: Vector of individual cluster IDs; nx1 or NULL. 
#'   }
#'   \item{\code{MCMC_out}}{List of full MCMC output, resulting from a call to 
#'   [run_MCMC_Rcpp()]}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling, resulting 
#'   from a call to [post_process()]}
#'   \item{\code{K_fixed}}{Number of global classes used for the fixed sampler}
#' }
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_wrpc_results.RData`. 
#' 
#' If only the adaptive sampler is run (i.e., `run_sampler` = `"adapt"`), returns
#' list `res` containing:
#' \describe{
#'   \item{\code{MCMC_out}}{List of full MCMC output, resulting from a call to 
#'   [run_MCMC_Rcpp()]}
#'   \item{\code{K_fixed}}{Number of global classes used for the fixed sampler, 
#' obtained from the adaptive sampler}
#'   \item{\code{K_MCMC}}{Adaptive sampler MCMC output for K; Mx1, where M is 
#'   the number of MCMC iterations after burn-in and thinning.}
#' }
#' If `save_res = TRUE` (default), also saves `res` as 
#' `[save_path]_wrpc_adapt.RData`. 
#' 
#' @importFrom RcppTN rtn
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats rnorm median confint
#' @importFrom survey svydesign svyglm
#' @export
#'
#' @examples
#' # Load libraries
#' library(baysc)  # Install at https://github.com/smwu/baysc
#' library(tidyverse)
#' 
#' # Source functions
#' source("Model_Code/wrpc_mcmc_fns.R")         # MCMC functions
#' Rcpp::sourceCpp("Model_Code/wrpc_mcmc.cpp")  # MCMC functions Rcpp code
#' source("Model_Code/wrpc_utilities.R")        # Additional utility functions
#' source("Model_Code/wrpc_plotting_fns.R")     # Plotting functions
#' 
#' ## Load and prepare data
#' # Load data example NHANES dataset (more information in `baysc` package)
#' load(file = "Model_Code/data_nhanes.rda")
#' # Categorical exposure matrix of food groups, nxJ
#' x_mat <- as.matrix(data_nhanes[, c(11:38)])
#' # Subpopulation indicators using age, nx1
#' h_all <- as.numeric(data_nhanes$age_cat)
#' # Binary outcome of hypertension
#' y_all <- data_nhanes$BP_flag
#' # Survey stratum indicators
#' stratum_id <- data_nhanes$stratum_id
#' # Survey cluster indicatos
#' cluster_id <- data_nhanes$cluster_id
#' # Survey sampling weights
#' sampling_wt <- data_nhanes$sample_wt
#' 
#' ## Run wrpc
#' res <- wrpc(x_mat = x_mat, h_all = h_all, sampling_wt = sampling_wt, 
#'              cluster_id = cluster_id, stratum_id = stratum_id, 
#'              run_sampler = "both", adapt_seed = 1, class_cutoff_global = 0.05, 
#'              n_runs = 50, burn = 25, switch = 10, thin = 1, save_res = FALSE)
#' 
#' ## Plot results
#' # Global pattern profiles              
#' plot_wrpc_global_pattern_profiles(res = res)
#' # Global pattern item consumption level probabilities
#' plot_wrpc_global_pattern_probs(res = res) 
#' # Local pattern deviations for items allocated to local (nu < 0.5)
#' plot_wrpc_local_profiles_allocation(res = res)
#' # Full local pattern profiles
#' plot_wrpc_local_pattern_profiles(res = res)
#' # Global-local deviation
#' plot_wrpc_allocation(res = res)        
#' # Distribution of global classes by subpopulation            
#' plot_wrpc_class_subgroup_dist(res = res)
#' 
#' 
#' ## Example for unweighted model
#' # Leave stratum_id, cluster_id, and sampling_wt as NULL
#' res_unwt <- wrpc(x_mat = x_mat, h_all = h_all, run_sampler = "both", 
#'                  adapt_seed = 1, class_cutoff_global = 0.05, n_runs = 50, 
#'                  burn = 25, switch = 10, thin = 1, save_res = FALSE)
#'
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
  # Initialize sampler name to save file
  sampler_name <- ""
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
  
  #================= Catch errors ==============================================
  catch_errors_wrpc(x_mat = x_mat, h_all = h_all, sampling_wt = sampling_wt, 
                     cluster_id = cluster_id, stratum_id = stratum_id, 
                     run_sampler = run_sampler, K_max = K_max, 
                     adapt_seed = adapt_seed, K_fixed = K_fixed, 
                     fixed_seed = fixed_seed, class_cutoff_global = class_cutoff_global, 
                     n_runs = n_runs, burn = burn, thin = thin, 
                    switch = switch, update = update,
                     save_res = save_res, save_path = save_path,
                     alpha_adapt = alpha_adapt, 
                     eta_global_adapt = eta_global_adapt, 
                     eta_local_adapt = eta_local_adapt, 
                     a_adapt = a_adapt, b_adapt = b_adapt,  
                     alpha_fixed = alpha_fixed, 
                     eta_global_fixed = eta_global_fixed, 
                     eta_local_fixed = eta_local_fixed, 
                     a_fixed = a_fixed, b_fixed = b_fixed)
  
  
  #================= ADAPTIVE SAMPLER ==========================================
  
  ### Run adaptive sampler to obtain number of latent classes, K_fixed
  if (run_sampler %in% c("both", "adapt")) { 
    sampler_name <- "_adapt"
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
      save(res, file = paste0(save_path, "_wrpc", sampler_name, ".RData"))
    }
    
    # Reduce memory burden
    rm(RPC_params, MCMC_out)
  }
  
  
  #================= FIXED SAMPLER =============================================
  
  ### Run fixed sampler to obtain posteriors
  if (run_sampler %in% c("both", "fixed")) {
    sampler_name <- "_results"
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
    
    # # Save output
    # if (save_res) {
    #   save(MCMC_out, file = paste0(save_path, "_wrpc_MCMC_out.RData"))
    # }
    
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
                      sampling_wt = sampling_wt, x_mat = x_mat, h_all = h_all,
                      stratum_id = stratum_id, cluster_id = cluster_id)
    res$data_vars <- data_vars
    
    class(res) <- "wrpc"
    
  }
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  res$runtime <- runtime
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_wrpc", sampler_name, ".RData"))
  }
  
  # Return output
  return(res)
}



