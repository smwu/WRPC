#=================================================
# Run WRPC model on simulated data
# Author: Stephanie Wu
# Date updated: 2024/11/25
#=================================================

# Read in two command line arguments
args <- commandArgs(trailingOnly = TRUE)
scenario <- args[[1]]   # Simulation scenario
samp_i <- args[[2]]     # Sample number
samp_i <- as.numeric(samp_i)


# Load libraries
library(baysc)
library(BART)  # BART
library(abind)  # binding arrays
library(parallel)  # parallel processing
library(plyr)
library(dplyr)
library(LaplacesDemon)  # distributions
library(truncnorm)
library(fastDummies)
library(matrixStats)
library(Matrix)
library(gtools)
library(e1071)
library(rstan)  # stan
library(survey)  # survey functions
library(Rcpp)  # Rcpp interface
library(RcppArmadillo)
library(RcppTN)

# Set directories
wd <- "~/Documents/GitHub/WRPC/"  # Working directory
wd <- "/n/netscratch/stephenson_lab/Lab/stephwu18/WRPC/"
data_dir <- "Simulations/Data/"    # Data directory
res_dir <- "Simulations/Results/"  # Results directory
code_dir <- "Model_Code/"  # Model code directory

# Create scenario results folder if it doesn't exist
dir_path <- paste0(wd, res_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

# Define path to save results
save_path <- paste0(wd, res_dir, "scen_", scenario, "/samp_", samp_i)

# Check if results already exist
already_done <- file.exists(paste0(save_path, "_wrpc_results.RData"))
if (already_done) {
  print(paste0('Scenario ', scenario, ' samp ', samp_i, ' already exists.'))
} else {
  
  print(paste0('Running scenario ', scenario, ' samp ', samp_i, '...'))
  
  # Source R model functions
  source(paste0(wd, code_dir, "wrpc.R"))                  # Main WRPC function
  source(paste0(wd, code_dir, "wrpc_mcmc_fns.R"))         # MCMC functions
  Rcpp::sourceCpp(paste0(wd, code_dir, "wrpc_mcmc.cpp"))  # MCMC functions Rcpp code
  source(paste0(wd, code_dir, "wrpc_utilities.R"))        # Additional utility functions
  
  ### Read in population and sample data
  load(paste0(wd, data_dir, "scen_", scenario, "/sim_pop.RData"))
  load(paste0(wd, data_dir, "scen_", scenario, "/", samp_i, "sim_samp.RData"))
  
  
  ### Define data and parameters
  # Estimating pseudo-weights
  x_mat <- sim_samp$X_data  # Multivariate categorical variables
  h_all <- sim_samp$true_Hi # Subpopulation indicator
  sampling_wt <- sim_samp$sample_wt # Sampling weight
  
  K_max <- 20          # Max number of latent classes for adaptive sampler
  adapt_seed <- samp_i # Seed for adaptive sampler
  fixed_seed <- 1    # Seed for fixed sampler
  n_runs <- 20000    # Number of MCMC iterations
  burn <- 10000      # Burn-in
  thin <- 5          # Thinning
  switch <- 100
  update <- 5000     # Display update frequency
  class_cutoff_global <- 0.05  # Cutoff proportion for global classes
  
  # # For tests
  # n_runs <- 1000
  # burn <- 500
  # thin <- 5
  # update <- 100
  
  
  ### Run weighted model
  res <- wrpc(x_mat = x_mat, h_all = h_all, sampling_wt = sampling_wt, 
              K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
              fixed_seed = fixed_seed, run_sampler = "both", 
              class_cutoff_global = 0.05, switch = switch, 
              n_runs = n_runs, burn = burn, thin = thin, update = update,
              save_res = TRUE, save_path = save_path)
  
  ### Run unweighted model
  res_unwt <- wrpc(x_mat = x_mat, h_all = h_all, sampling_wt = rep(1, nrow(x_mat)), 
                   K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
                   fixed_seed = fixed_seed, run_sampler = "both", 
                   class_cutoff_global = 0.05, switch = switch, 
                   n_runs = n_runs, burn = burn, thin = thin, update = update,
                   save_res = TRUE, save_path = paste0(save_path, "_unwt"))
}

