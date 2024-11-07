#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


//' Draw from a Dirichlet distribution
//' 
//' @param alpha Vector of positive concentration parameters. Length determines
//' the number of dimensions for the Dirichlet distribution.
//' @return Vector result of a single draw from the specified Dirichlet distribution
//' @importFrom stats rgamma
//' @keywords internal
// [[Rcpp::export]]
arma::vec rdirichlet_cpp(arma::vec alpha) {
 int distribution_size = alpha.size();
 // Draw from a Dirichlet
 arma::vec distribution(distribution_size);
 
 double sum_term = 0;
 // Loop through the distribution and draw Gamma variables
 for (int i = 0; i < distribution_size; i++) {
   // Use R instead of Rcpp b/c faster for a scalar
   // For the R functions, omit the "n" argument
   // Source: https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html
   double cur = R::rgamma(alpha[i], 1.0);
   distribution(i) = cur;
   sum_term += cur;
 }
 // Now normalize
 for (int i = 0; i < distribution_size; i++) {
   distribution(i) = distribution(i)/sum_term;
 }
 return distribution;
}
 
 
//' Draw from a Categorical distribution
//' 
//' @param probs Row vector of category event probabilities. Length determines
//' the number of categories for the Categorical distribution.
//' @return Integer specifying the category (1-based indexing) of a single draw 
//' from the specified Categorical distribution
//' @importFrom stats rmultinom
//' @keywords internal
// [[Rcpp::export]]
int rcat_cpp(arma::rowvec probs) {
 int num_categs = probs.size();
 IntegerVector draw(num_categs);
 R::rmultinom(1, probs.begin(), num_categs, draw.begin());
 int categ = which_max(draw);
 return categ; 
}


//' Apply log-sum-exp trick 
//' 
//' @description
//' `logSumExp_cpp` computes the logarithm of the sum of exponentials, adapted from 
//' https://github.com/helske/seqHMM/blob/main/src/logSumExp.cpp.
//'
//' @param x Row vector of input values
//' @return Result of computing the logarithm of the sum of exponentials of the
//' input values.
//' @keywords internal
// [[Rcpp::export]]
double logSumExp_cpp(const arma::rowvec& x) {
 int maxi = x.index_max();
 double maxv = x(maxi);
 if (!(maxv > -arma::datum::inf)) {
   return -arma::datum::inf;
 }
 double cumsum = 0.0;
 for (unsigned int i = 0; i < x.n_elem; i++) {
   if ((i != maxi) && (x(i) > -arma::datum::inf)) {
     cumsum += exp(x(i) - maxv);
   }
 }
 return maxv + log1p(cumsum);
}


//' Update g for WRPC
//' 
//' `update_g_WRPC` updates the global-local assignment for each individual and 
//' item.
//' 
//' @inheritParams update_c_WRPC
//' 
//' @return Updated `g_mat` matrix of global-local assignments by drawing from 
//' Bernoulli distributions with updated global assignment probabilities
//' @keywords internal
// [[Rcpp::export]]
arma::mat update_g_WRPC(arma::mat& g_mat, const int& J, const int& n, 
                       const arma::vec& h_all, const arma::vec& c_all, 
                       const arma::mat& x_mat, 
                       const arma::mat& nu, const arma::cube& theta_global, 
                       const arma::cube& theta_local) {
 
 // Calculate posterior global-local assignment, for each individual and item
 for (int i = 0; i < n; i++) {
   // Get subpopulation index. Careful about 0-based indexing
   int h_i_ind = h_all(i) - 1;  
   // Posterior global-local assignment for each item
   for (int j = 0; j < J; j++) {
     double global = nu(h_i_ind, j) * 
       theta_global(j, c_all(i) - 1, x_mat(i, j) - 1);
     double local = (1 - nu(h_i_ind, j)) * 
       theta_local(j, h_i_ind, x_mat(i, j) - 1);
     double p_ij = global / (global + local);
     // Faster in R for scalar
     g_mat(i, j) = R::rbinom(1, p_ij);
   }
 }
 return g_mat;
}


//' Update pi for WRPC
//' 
//' `update_pi_WRPC` updates the vector parameter of global class membership 
//' probabilities by drawing from its posterior.
//' 
//' @inheritParams run_MCMC_Rcpp_WRPC
//' @param pi Vector parameter of global class membership probabilities. Kx1
//' 
//' @return Updated `pi` vector after drawing from its posterior distribution
//' @keywords internal
// [[Rcpp::export]]
arma::vec update_pi_WRPC(arma::vec& pi, const arma::vec& w_all, 
                          const arma::vec& c_all, const int& K, 
                          const arma::vec& alpha) {
 NumericVector w_all_copy = as<NumericVector>(wrap(w_all));
 // Posterior parameters for pi
 arma::vec alpha_post(K);  
 for (int k = 0; k < K; k++) {
   // Add sum of normalized weights for those assigned to class k, equiv. to
   // weighted number of individuals assigned to each class
   // Careful with 0-based indexing
   LogicalVector indiv_k = (as<IntegerVector>(wrap(c_all)) == (k + 1));
   NumericVector weights_k = w_all_copy[(indiv_k)];
   alpha_post[k] = alpha[k] + sum(weights_k);
 }

 // Draw pi from posterior
 arma::vec out = rdirichlet_cpp(alpha_post);
 pi = out;
 
 return pi;
}


//' Update c for WRPC
//' 
//' `update_c_WRPC` updates the vector of global individual class assignments 
//' by drawing from a Categorical distribution with updated category probabilities.
//' 
//' @inheritParams run_MCMC_Rcpp
//' @param pi_global Vector parameter of global class membership probabilities. Kx1
//' @param theta_global Array parameter of global item level probabilities. JxKxR
//' @param g_mat Matrix of global-local assignments for each individual and item. nxJ
//' 
//' @return Updated `c_all` vector after drawing from a Categorical distribution
//' with updated category event probabilities.
//' @keywords internal
// [[Rcpp::export]]
arma::vec update_c_WRPC(arma::vec& c_all, const int& n, const int& K, const int& J, 
                        const arma::cube& theta_global, const arma::mat& x_mat, 
                        const arma::vec& pi, const arma::mat& g_mat) {
 arma::mat log_cond_c(n, K);        // Individual log-likelihood for each class
 arma::mat pred_class_probs(n, K);  // Posterior class membership probabilities
 
 // Calculate posterior class membership, p(c_i=k|-), for each class k and
 // update class assignments
 for (int i = 0; i < n; i++) {
   for (int k = 0; k < K; k++) {
     // Calculate theta component of individual log-likelihood for class k
     double log_theta_comp_k = 0.0;
     for (int j = 0; j < J; j++) {
       // Add in global thetas for items with global assignment
       if (g_mat(i, j) == 1) {
         // Subtract 1 from exposure value due to 0-based indexing
         log_theta_comp_k += log(theta_global(j, k, x_mat(i, j) - 1));
       }
     }
     // Individual log-likelihood for class k
     log_cond_c(i, k) = log(pi(k)) + log_theta_comp_k;
   }
   // Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
   pred_class_probs.row(i) = exp(log_cond_c.row(i) - logSumExp_cpp(log_cond_c.row(i)));
   // Update class assignment using the posterior probabilities
   // Be careful of 0-based indexing
   c_all(i) = rcat_cpp(pred_class_probs.row(i)) + 1;
 }
 return c_all;
}


//' Update global theta for WRPC
//' 
//' `update_theta_global_WRPC` updates the array of global item level 
//' probabilities by drawing from its posterior.
//' 
//' @inheritParams update_c_WRPC
//' 
//' @return Updated `theta_global` array after drawing from its posterior distribution.
//' @keywords internal
// [[Rcpp::export]]
arma::cube update_theta_global_WRPC(arma::cube& theta_global, const int& J, 
                                    const int& K, const int& R, 
                                    const arma::mat& eta_global, 
                                    const arma::vec& w_all, const arma::vec& c_all, 
                                    const arma::mat& x_mat, const arma::mat& g_mat) {
 NumericVector w_all_copy = as<NumericVector>(wrap(w_all));
 
 for (int j = 0; j < J; j++) {
   NumericVector eta_global_post(R);  // Posterior parameters for theta for item j
   // Select those with item j assigned global
   LogicalVector indiv_global = (as<NumericVector>(wrap(g_mat.col(j))) == 1);
   for (int k = 0; k < K; k++) {
     LogicalVector indiv_k = (as<IntegerVector>(wrap(c_all)) == (k + 1));
     for (int r = 0; r < R; r++) {
       // Add sum of normalized weights for those assigned to class k with 
       // x_ij = r and global assignment for item j
       // Careful with 0-based indexing
       LogicalVector indiv_r = (as<NumericVector>(wrap(x_mat.col(j))) == (r + 1));
       NumericVector weights_r = w_all_copy[(indiv_global & indiv_k & indiv_r)];
       eta_global_post(r) = eta_global(j, r) + sum(weights_r);
     }
     // Draw theta from posterior
     theta_global.tube(j, k) = rdirichlet_cpp(eta_global_post);
   }
 }
 return theta_global;
}


//' Update local theta for WRPC
//' 
//' `update_theta_local_WRPC` updates the array of local item level 
//' probabilities by drawing from its posterior.
//' 
//' @inheritParams update_c_WRPC
//' @param eta_local JxHxR array of hyperparameters for theta_local
//' 
//' @return Updated `theta_local` array after drawing from its posterior distribution.
//' @keywords internal
// [[Rcpp::export]]
arma::cube update_theta_local_WRPC(arma::cube& theta_local, const int& H, 
                                   const int& J, const int& R, 
                                   const arma::cube& eta_local, 
                                   const arma::vec& w_all, const arma::vec& h_all,
                                   const arma::mat& x_mat, const arma::mat& g_mat) {
 NumericVector w_all_copy = as<NumericVector>(wrap(w_all));
 
 for (int j = 0; j < J; j++) {
   NumericVector eta_local_post(R);  // Posterior parameters for theta for item j
   // Select those with item j assigned local
   LogicalVector indiv_local = (as<NumericVector>(wrap(g_mat.col(j))) == 0);
   for (int h = 0; h < H; h++) {
     // Select those in subpopulation h 
     LogicalVector indiv_h = (as<IntegerVector>(wrap(h_all)) == (h + 1));
     for (int r = 0; r < R; r++) {
       // Add sum of normalized weights for those assigned to class k with 
       // x_ij = r in subpop h. Careful with 0-based indexing
       LogicalVector indiv_r = (as<NumericVector>(wrap(x_mat.col(j))) == (r + 1));
       NumericVector weights_r = w_all_copy[(indiv_h & indiv_local & indiv_r)];
       eta_local_post(r) = eta_local(j, h, r) + sum(weights_r);
     }
     // Draw theta from posterior
     theta_local.tube(j, h) = rdirichlet_cpp(eta_local_post);
   }
 }
 return theta_local;
}


//' Update nu for WRPC
//' 
//' `update_nu_WRPC` updates the global-local assignment probabilities, nu, by
//' drawing from its posterior.
//' 
//' @inheritParams update_c
//' 
//' @return Updated `nu` array after drawing from its posterior distribution.
//' @keywords internal
// [[Rcpp::export]]
arma::mat update_nu_WRPC(arma::mat& nu, const int&H, const int& J, const int& n, 
                         const arma::mat& a, const arma::mat& b, 
                         const arma::vec& w_all, const arma::vec& h_all, 
                         const arma::mat& g_mat) {
 
 NumericVector w_all_copy = as<NumericVector>(wrap(w_all));
 double a_post;
 double b_post;
 
 // Update for each subgroup h and item j
 for (int h = 0; h < H; h++) {
   for (int j = 0; j < J; j++) {
     // Get individuals in subgroup h, separating out global vs local for j
     LogicalVector indiv_h = (as<IntegerVector>(wrap(h_all)) == (h + 1));
     LogicalVector indiv_global = (as<IntegerVector>(wrap(g_mat.col(j))) == 1);
     LogicalVector indiv_local = (as<IntegerVector>(wrap(g_mat.col(j))) == 0);
     // Weights for subsetted individuals
     NumericVector weights_h_global = w_all_copy[(indiv_h & indiv_global)];
     NumericVector weights_h_local = w_all_copy[(indiv_h & indiv_local)];
     
     // Posterior hyperparameters
     a_post = a(h, j) + sum(weights_h_global);
     b_post = b(h, j) + sum(weights_h_local);
     
     // Draw nu from posterior. Faster in R for scalar 
     nu(h, j) = R::rbeta(a_post, b_post);
   }
 }
 return nu;
}



/*** R
# ### TESTING CODE
# 
# setwd(wd)
# code_dir <- "Code/Nimble/"
# data_dir <- "Code/Nimble/"
# res_dir <- "Code/Nimble/"
# # Read in data
# load(paste0(wd, data_dir, "w1_sim_samp.RData"))
# data_vars <- sim_samp
# 
# # Fixed K
# K <- 3
# 
# # Run WRPC
# x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
# stratum_id <- NULL      # Stratum indicators, nx1
# cluster_id <- NULL                   # No clustering
# sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
# h_all <- data_vars$true_Hi           # Subpopulation indicators, nx1
# 
# print("Read in data")
# # Obtain dimensions
# n <- dim(x_mat)[1]        # Number of individuals
# J <- dim(x_mat)[2]        # Number of exposure items
# R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
#              function(x) length(unique(x)))  
# R <- max(R_j)             # Maximum number of exposure categories across items
# H <- length(unique(h_all)) # Number of subpopulations
# # Obtain normalized weights
# kappa <- sum(sampling_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
# w_all <- c(sampling_wt / kappa) # Weights normalized to sum to n, nx1
# 
# ### Set hyperparameters
# a <- b <- matrix(1, nrow = H, ncol = J)
# alpha <- rep(1, K) / K
# eta_global <- matrix(0.01, nrow = J, ncol = R)
# for (j in 1:J) {
#   eta_global[j, 1:R_j[j]] <- rep(1, R_j[j])
# }
# eta_local <- array(0.01, dim = c(J, H, R))
# for (j in 1:J) {
#   for (h in 1:H) {
#     eta_local[j, h, 1:R_j[j]] <- rep(1, R_j[j])
#   }
# }
# 
# ### Initialize model
# # Initialize at different values
# nu <- matrix(0.5, nrow = H, ncol = J)
# pi <- rep(1, times = K) / K
# c_all <- rep(1:K, length.out = nrow(x_mat))
# theta_global <- array(1 / R, dim = c(J, K, R))
# theta_local <- array(1 / R, dim = c(J, H, R))
# g_mat <- matrix(rbinom(n = nrow(x_mat), size = 1, prob = 0.5), 
#                 nrow = nrow(x_mat), ncol = J)
#
#
#
# ### Test g_mat
# g_mat[1:5, 1:5]
# ## Test: Rcpp version
# set.seed(1)
# g_mat <- update_g_WRPC(g_mat = g_mat, J = J, n = n, h_all = h_all, c_all = c_all, 
#                         nu = nu, x_mat = x_mat, theta_global = theta_global, 
#                         theta_local = theta_local)
# g_mat[1:5, 1:5]
# 
# g_mat <- matrix(rbinom(n = nrow(x_mat), size = 1, prob = 0.5), 
#                 nrow = nrow(x_mat), ncol = J)  # Reset
# ## Update g_mat: R version 
# set.seed(1)
# for (i in 1:n) {
#   h_i <- h_all[i]
#   p_i <- numeric(J)
#   # Posterior global-local assignment probability for each item
#   for (j in 1:J) {
#     global_prob <- nu[h_i, j] * theta_global[j, c_all[i], x_mat[i, j]]
#     local_prob <- (1 - nu[h_i, j]) * theta_local[j, h_i, x_mat[i, j]]
#     p_i[j] <- global_prob / (global_prob + local_prob)
#   }
#   # For each individual, assign items to global or local
#   g_mat[i, ] <- rbinom(n = J, size = 1, prob = p_i)
# }
# g_mat[1:5, 1:5]
# 
# 
# ### Test pi
# pi
# ## Test: Rcpp version
# set.seed(1)
# pi <- update_pi_WRPC(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
# pi
# 
# pi <- rep(1, times = K) / K
# ## Update pi: R version
# set.seed(1)
# alpha_post <- numeric(K)
# for (k in 1:K) {
#   # Add sum of normalized weights for those assigned to class k, equiv. to
#   # weighted number of individuals assigned to each class
#   alpha_post[k] <- alpha[k] + sum(w_all[c_all == k])
# }
# pi <- c(gtools::rdirichlet(n = 1, alpha = alpha_post))
# pi
# 
# 
# ### Test c_all
# c_all[1:10]
# ## Test: Rcpp version
# set.seed(1)
# c_all <- update_c_WRPC(c_all = c_all, n = n, K = K, J = J, 
#                        theta_global = theta_global, x_mat = x_mat, 
#                        pi = pi, g_mat = g_mat)
# c_all[1:10]
# 
# c_all <- rep(1:K, length.out = nrow(x_mat))
# ## Update c_all: R version
# set.seed(1)
# # Posterior class membership probabilities
# pred_global_class_probs <- matrix(NA, nrow = n, ncol = K) 
# # Individual log-likelihood for each class
# log_cond_c <- matrix(NA, nrow = n, ncol = K)       
# # Calculate posterior class membership, p(c_i=k|-), for each class k
# for (i in 1:n) {
#   for (k in 1:K) {
#     # Calculate theta component of individual log-likelihood assuming class k
#     log_global_theta_comp_k <- 0
#     for (j in 1:J) {
#       # Add in global thetas for items with global assignment
#       if (g_mat[i, j] == 1) {
#         log_global_theta_comp_k <- log_global_theta_comp_k + 
#           log(theta_global[j, k, x_mat[i, j]])
#       }
#     }
#     # Individual log-likelihood for class k
#     log_cond_c[i, k] <- log(pi[k]) + log_global_theta_comp_k
#   }
#   # Calculate p(c_i=k|-) = p(x,c_i=k) / p(x)
#   pred_global_class_probs[i, ] <- 
#     exp(log_cond_c[i, ] - matrixStats::logSumExp(log_cond_c[i, ], na.rm = TRUE))
#   # Update class assignment using the posterior probabilities
#   c_all[i] <- which(stats::rmultinom(n = 1, size = 1, 
#                                      prob = pred_global_class_probs[i, ]) == 1)
# }
# c_all[1:10]
# 
# 
# ### Test theta_global
# theta_global[1:5, 1:3, 1]
# ## Test: Rcpp version
# set.seed(1)
# theta_global <- update_theta_global_WRPC(theta_global = theta_global, J = J, 
#                                          K = K, R = R, eta_global = eta_global, 
#                                          w_all = w_all, c_all = c_all, 
#                                          x_mat = x_mat, g_mat = g_mat)
# theta_global[1:5, 1:3, 1]
# 
# theta_global <- array(1 / R, dim = c(J, K, R))
# ## Update theta_global: R version
# set.seed(1)
# eta_global_post <- numeric(R)
# for (j in 1:J) {
#   for (k in 1:K) {
#     for (r in 1:R) {
#       # Add sum of normalized weights for those assigned to class k with x_ij = r
#       eta_global_post[r] <- eta_global[r] + 
#         sum(w_all[(g_mat[, j] == 1) & (c_all == k) & (x_mat[, j] == r)])
#     }
#     theta_global[j, k, ] <- c(gtools::rdirichlet(n = 1, alpha = eta_global_post))
#   }
# }
# theta_global[1:5, 1:3, 1]
# 
# 
# ### Test theta_local
# theta_local[1:5, 1:2, 1]
# ## Test: Rcpp version
# set.seed(1)
# theta_local <- update_theta_local_WRPC(theta_local = theta_local, H = H, 
#                                         J = J, R = R, eta_local = eta_local, 
#                                         w_all = w_all, h_all = h_all, 
#                                         x_mat = x_mat, g_mat = g_mat)
# theta_local[1:5, 1:2, 1]
# 
# theta_local <- array(1 / R, dim = c(J, H, R))
# ## Update theta_local: R version
# set.seed(1)
# eta_local_post <- numeric(R)
# for (j in 1:J) {
#   for (h in 1:H) {
#     for (r in 1:R) {
#       # Add sum of normalized weights for those assigned to class k with x_ij = r
#       eta_local_post[r] <- eta_local[r] + 
#         sum(w_all[(g_mat[, j] == 0) & (h_all == h) & (x_mat[, j] == r)])
#     }
#     theta_local[j, h, ] <- c(gtools::rdirichlet(n = 1, alpha = eta_local_post))
#   }
# }
# theta_local[1:5, 1:2, 1]
# 
# 
# ### Test nu
# nu[1:2, 1:5]
# ## Test: Rcpp version
# set.seed(1)
# nu <- update_nu_WRPC(nu = nu, H = H, J = J, n = n, a = a, b = b, w_all = w_all, 
#                      h_all = h_all, g_mat = g_mat)
# nu[1:2, 1:5]
# 
# nu <- matrix(0.5, nrow = H, ncol = J)
# ## Update nu: R version
# set.seed(1)
# for (h in 1:H) {
#   for (j in 1:J) {
#     a_post <- a[h, j] + sum(w_all[(g_mat[, j] == 1) & (h_all == h)])
#     b_post <- b[h, j] + sum(w_all[(g_mat[, j] == 0) & (h_all == h)])
#     nu[h, j] <- stats::rbeta(n = 1, shape1 = a_post, shape2 = b_post)
#   }
# }
# nu[1:2, 1:5]

*/



