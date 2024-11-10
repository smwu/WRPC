#============
# WRPC Tests
#============

rm(list=ls())

library(tidyverse)
library(ggpubr)
library(gtools)
library(testthat)
library(sfsmisc)  # mult.fig

wd <- "~/Documents/Github/WRPC/"
code_dir <- "Code/"
data_dir <- "Code/Nimble/"
res_dir <- "Results/"

# Source functions
source(paste0(wd, code_dir, "wrpc.R"))
source(paste0(wd, code_dir, "wrpc_mcmc_fns.R"))
Rcpp::sourceCpp(paste0(wd, code_dir, "wrpc_mcmc.cpp"))

# Read in data
load(paste0(wd, data_dir, "w1_sim_pop.RData"))
load(paste0(wd, data_dir, "w1_sim_samp.RData"))
data_vars <- sim_samp

#=============== Test adaptive sampler =========================================
# Overfitted K
K_max <- 30

# Run WRPC
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
stratum_id <- NULL      # Stratum indicators, nx1
cluster_id <- NULL                   # No clustering
sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
h_all <- data_vars$true_Hi           # Subpopulation indicators, nx1

n_runs <- 2000
burn <- 0
thin <- 1
update <- n_runs / 10
switch <- 50
class_cutoff_global <- 0.05

set.seed(20241102)
res_adapt <- wrpc(K_max = K_max, x_mat = x_mat, h_all = h_all, 
                  sampling_wt = sampling_wt, run_sampler = "adapt", 
                  adapt_seed = 1, n_runs = n_runs, burn = burn, thin = thin, 
                  update = update, save_res = TRUE, 
                  save_path = paste0(wd, res_dir, "test_adapt"))

test_that("Adaptive sampler works", {
  expect_equal(res_adapt$K_fixed, 3) 
})


#=============== Test adaptive and fixed sampler ===============================
# Overfitted K
K_max <- 20

# Run WRPC
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxJ
stratum_id <- NULL      # Stratum indicators, nx1
cluster_id <- NULL                   # No clustering
sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
h_all <- data_vars$true_Hi           # Subpopulation indicators, nx1

n_runs <- 2000
burn <- 0
thin <- 1
update <- n_runs / 10
switch <- 50
class_cutoff_global <- 0.05

set.seed(20241102)
res_both <- wrpc(K_max = K_max, x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, run_sampler = "both", 
                 adapt_seed = 1, n_runs = n_runs, burn = burn, thin = thin, 
                 update = update, save_res = TRUE,
                 save_path = paste0(wd, res_dir, "test_Rcpp_adapt_fixed"))

test_that("Adaptive and fixed sampler work", {
  expect_equal(res_both$K_fixed, 3) 
  expect_equal(res_both$estimates$K_red, 3)
})


#================ Examine results ==============================================
# Define true parameters
parameter_true <- list(
  c_all = sim_pop$true_Ci,
  g_mat = sim_pop$true_Gij,
  pi = as.vector(sim_pop$true_pi_global),
  theta_global = sim_pop$true_global_thetas,
  theta_local = sim_pop$true_local_thetas,
  nu = sim_pop$true_nu
)

# Reorder to match true class order
ordering <- c(2, 1, 3)
theta_global_plot <- res_both$estimates$theta_global_red[, , ordering, ]
theta_local_plot <- res_both$estimates$theta_local_red
pi_plot <- res_both$estimates$pi_red[, ordering]
nu_plot <- res_both$estimates$nu_red
c_plot <- ifelse(res_both$estimates$c_all == ordering[1], 1, 
                 ifelse(res_both$estimates$c_all == ordering[2], 2, 3))
g_mat <- res_both$estimates$g_mat
K <- res_both$estimates$K_red
theta_global_modes <- res_both$estimates$theta_global_modes
theta_local_modes <- res_both$estimates$theta_local_modes
nu_red <- res_both$estimates$nu_red
H <- res_both$data_vars$H
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
# Check misclassification
table(sim_samp$true_Ci, c_plot)
prop.table(table(sim_samp$true_Ci, c_plot))

## look at g_mat --------
# Check misclassification
table(c(sim_samp$true_Gij), c(g_mat))
prop.table(table(c(sim_samp$true_Gij), c(g_mat)))
