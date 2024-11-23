#' Catch errors in input variables for external functions
#' 
#' @description
#' Catch input errors in variables necessary for package functions. All 
#' parameters are set to `NULL` by default so that error checks 
#' are only performed on relevant variables.
#' 
#' @inheritParams swolca
#' @param model String specifying which model is used. Must be one of `swolca` 
#' (default) or `wolca`
#' @return The function stops and an error message is displayed if the input 
#' variables are not acceptable
#' @details All parameters are set to `NULL` by default so that error checks 
#' are only performed on relevant variables.
#' 
#' @importFrom stats terms as.formula
#' @importFrom stringr str_detect
#' @keywords internal
#' @export

catch_errors_wrpc <- function(x_mat = NULL, h_all = NULL, sampling_wt = NULL, 
                              cluster_id = NULL, stratum_id = NULL, 
                              run_sampler = NULL, K_max = NULL, L_max = NULL,
                              adapt_seed = NULL, K_fixed = NULL, 
                              fixed_seed = NULL, class_cutoff_global = NULL, 
                              n_runs = NULL, burn = NULL, thin = NULL, update = NULL, 
                              num_reps = NULL, save_res = NULL, save_path = NULL,
                              alpha_adapt = NULL, eta_global_adapt = NULL, 
                              eta_local_adapt = NULL, a_adapt = NULL, b_adapt = NULL,
                              alpha_fixed = NULL, eta_global_fixed = NULL, 
                              eta_local_fixed = NULL, a_fixed = NULL, b_fixed = NULL) {
  if (is.null(x_mat) | is.null(h_all)) {
    stop("need to specify exposure matrix, x_mat, and subpopulation vector, h_all")
  } else {
    
    # Check type for x_mat
    if (!is.matrix(x_mat)) {
      stop("x_mat must be a numeric matrix")
    }
    # Check type for h_all
    if (!is.vector(h_all)) {
      stop("h_all must be a numeric vector")
    }
    
    # Obtain dimensions
    n <- dim(x_mat)[1]        # Number of individuals
    J <- dim(x_mat)[2]        # Number of exposure items
    R <- max(apply(x_mat, 2,  # Number of exposure categories
                   function(x) length(unique(x))))  
    H <- length(unique(h_all)) # Number of subpopulations
    
    # Check sampler specification
    if (!is.null(run_sampler)) {
      if (!(run_sampler %in% c("both", "adapt", "fixed"))) {
        stop("run_sampler must be one of `both`, `adapt`, or `fixed`")
      }
      if (run_sampler == "fixed") {
        if (is.null(K_fixed)) {
          stop("K_fixed must be specified")
        }
      }
    }
    
    # Check same number of individuals for x and h
    if (n != length(h_all)) {
      stop("number of rows in x_mat must match length of h_all")
    }
    # Check h_all coding
    if (!all(unique(h_all) %in% 1:H)) {
      stop(paste0("please recode h_all so that its values range from 1 to ", H))
    }
    
    # Check same number of individuals for x and sampling weights
    if (!is.null(sampling_wt)) {
      if (!is.vector(sampling_wt)) {
        stop("sampling_wt must be a numeric vector")
      }
      if (n != length(sampling_wt)) {
        stop("number of rows in x_mat must match length of sampling_wt")
      }
    }
    
    # If no clustering, assign each individual to their own cluster. Else, check
    # same number of individuals for x and clusters
    if (is.null(cluster_id)) {
      cluster_id <- 1:n
    } else {
      if (!is.vector(cluster_id)) {
        stop("cluster_id must be a numeric vector")
      }
      if (n != length(cluster_id)) {
        stop("number of rows in x_mat must match length of cluster_id")
      }
    }
    
    # Check same number of individuals for x and strata
    if (!is.null(stratum_id)) {
      if (!is.vector(stratum_id)) {
        stop("stratum_id must be a numeric vector")
      }
      if (n != length(stratum_id)) {
        stop("number of rows in x_mat must match length of stratum_id")
      }
    }
    
    # Check cutoff is between 0 and 1
    if (!is.null(class_cutoff_global)) {
      if (class_cutoff_global <= 0 | class_cutoff_global >= 1) {
        stop("class_cutoff_global must be a proportion in (0,1)")
      }
    }
    
    # Check seeds
    if (!is.null(adapt_seed)) {
      if (!is.numeric(adapt_seed)) {
        stop("adapt_seed must be numeric")
      }
    }
    if (!is.null(fixed_seed)) {
      if (!is.numeric(fixed_seed)) {
        stop("fixed_seed must be numeric")
      }
    }
    
    # Check hyperparameter dimensions for adaptive sampler
    if (!is.null(alpha_adapt)) {
      if (length(alpha_adapt) != K_max) {
        stop("length of alpha_adapt must be the same as K_max")
      }
    }
    if (!is.null(eta_global_adapt)) {
      if ((nrow(eta_global_adapt) != J) | (ncol(eta_global_adapt) != R)) {
        stop("eta_global_adapt must be a matrix with J rows and R columns, where J is
             the number of exposure items and R is the maximum number of 
             exposure categories")
      }
      if (any(eta_global_adapt == 0)) {
        warning("eta_global_adapt has 0 values and may result in rank-difficiency 
                issues during the Hessian calculation in the var_adjust() function")
      }
    }
    if (!is.null(eta_local_adapt)) {
      if ((dim(eta_local_adapt)[1] != J) | (dim(eta_local_adapt)[2] != H) | 
          (dim(eta_local_adapt)[3] != R)) {
        stop(paste0("eta_local_adapt must be an array of dimension JxHxR, ",
                    "where J is the number of exposure items, H is the number of ",
                    "subpopulations, and R is the maximum number of exposure categories"))
      }
      if (any(eta_local_adapt == 0)) {
        warning("eta_local_adapt has 0 values and may result in rank-difficiency 
                issues during the Hessian calculation in the var_adjust() function")
      }
    }
    if (!is.null(a_adapt)) {
      if (!is.numeric(a_adapt) | dim(a_adapt)[1] != H | dim(a_adapt)[2] != J) {
        stop(paste0("a_adapt must be a numeric matrix of dimension HxJ, ", 
                    "where H is the number of subpopulations and J is the ",
                    "number of exposure item"))
      }
    }
    if (!is.null(b_adapt)) {
      if (!is.numeric(b_adapt) | dim(b_adapt)[1] != H | dim(b_adapt)[2] != J) {
        stop(paste0("b_adapt must be a numeric matrix of dimension HxJ, ", 
                    "where H is the number of subpopulations and J is the ",
                    "number of exposure item"))
      }
    }
    
    # Check number of classes
    if (!is.null(K_max)) {
      if (K_max < 1) {
        stop("Maximum number of global classes in K_max must be at least 1")
      }
    }
    if (!is.null(K_fixed)) {
      if (K_fixed < 1) {
        stop("Maximum number of global classes in K_fixed must be at least 1")
      }
      # Check hyperparameter dimensions for fixed sampler
      if (!is.null(alpha_fixed)) {
        if (length(alpha_fixed) != K_fixed) {
          stop("length of alpha_fixed must be the same as K_fixed")
        }
      }
      if (!is.null(eta_global_fixed)) {
        if ((nrow(eta_global_fixed) != J) | (ncol(eta_global_fixed) != R)) {
          stop("eta_global_fixed must be a matrix with J rows and R columns, 
          where J is the number of exposure items and R is the maximum number of 
          exposure categories")
        }
        if (any(eta_global_fixed == 0)) {
          warning("eta_global_fixed has 0 values and may result in rank-difficiency issues
                during the Hessian calculation in the var_adjust() function")
        }
      }
    } else {
      if (any(!is.null(c(alpha_fixed, eta_global_fixed)))) {
        stop("K_fixed must be specified alongside alpha_fixed and eta_global_fixed")
      }
    }
    if (!is.null(eta_local_fixed)) {
      if ((dim(eta_local_fixed)[1] != J) | (dim(eta_local_fixed)[2] != H) | 
          (dim(eta_local_fixed)[3] != R)) {
        stop("eta_local_fixed must be am array with dimensions JxHxR")
      }
      if (any(eta_local_fixed == 0)) {
        warning("eta_local_fixed has 0 values and may result in rank-difficiency issues
              during the Hessian calculation in the var_adjust() function")
      }
    }
    if (!is.null(a_fixed)) {
      if (!is.numeric(a_fixed) | dim(a_fixed)[1] != H | dim(a_fixed)[2] != J) {
        stop(paste0("a_fixed must be a numeric matrix of dimension HxJ, ", 
                    "where H is the number of subpopulations and J is the ",
                    "number of exposure item"))
      }
    }
    if (!is.null(b_fixed)) {
      if (!is.numeric(b_fixed) | dim(b_fixed)[1] != H | dim(b_fixed)[2] != J) {
        stop(paste0("b_fixed must be a numeric matrix of dimension HxJ, ", 
                    "where H is the number of subpopulations and J is the ",
                    "number of exposure item"))
      }
    }
    
    # Check MCMC parameters
    if (!all(is.null(c(n_runs, burn, thin, switch, update)))) {
      if (!all(c(n_runs, burn, thin, switch, update) %% 1 == 0) | 
          !all(c(n_runs, burn, thin, switch, update) >= 0)) {
        stop("n_runs, burn, thin, and update must be whole numbers")
      }
      if (burn > n_runs) {
        stop("n_runs must be larger than burn")
      }
    }
    
    # Check saving parameters
    if (!is.null(save_res)) {
      if (!is.logical(save_res)) {
        stop("save_res must be a boolean specifying if results should be saved")
      }
      if (save_res) {
        if (is.null(save_path) | !is.character(save_path)) {
          stop("save_path must be a string specifying a path and file name, such as '~/Documents/run'")
        } else {
          last_slash_ind <- regexpr("\\/[^\\/]*$", save_path)
          if (!dir.exists(substr(save_path, start = 1, stop = last_slash_ind))) {
            stop("directory specified in save_path does not exist")
          }
          if (last_slash_ind == length(save_path)) {
            stop("please append the start of a file name to the end of save_path. 
            For example, '~/Documents/run' can produce a saved file named 
            'run_wrpc_results.RData'")
          }
        }
      }
    }
  }
}
