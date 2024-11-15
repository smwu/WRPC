#====================================================
# Analyzing HCHS data considering condition awareness
# including those with self-report
# Author: Stephanie Wu
# Date Updated: 2024/11/13
#====================================================

#==================== Read in HCHS data ========================================

rm(list = ls())

library(tidyverse)  # data wrangling
library(baysc)      # bayesian survey clustering
library(haven)      # read sas file
library(naniar)     # missingness
library(sfsmisc)    # mult.fig

wd <- "~/Documents/Github/WRPC/"
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WRPC/"
code_dir <- "Model_Code/"
data_dir <- "Application/HCHS_Data/"
res_dir <- "Results/"

# Read in raw derived variables
raw_derv_vars <- read_sas(paste0(wd, data_dir, "part_derv_inv4.sas7bdat"))
# Read in raw self-report variables
raw_mhea_vars <- read_sas(paste0(wd, data_dir, "mhea_inv4.sas7bdat"))
# Read in raw FPQ data
raw_fpq_vars <- read_sas(paste0(wd, data_dir, "fpe_inv4.sas7bdat"))
# Read in cleaned FPQ data
sol_data <- read.csv(paste0(wd, data_dir, "FPQ129_CVD07FEB2023.csv"))

# Source WRPC functions
source(paste0(wd, code_dir, "wrpc.R"))
source(paste0(wd, code_dir, "wrpc_mcmc_fns.R"))
Rcpp::sourceCpp(paste0(wd, code_dir, "wrpc_mcmc.cpp"))

# Read in cleaned data
study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_unaware_full.csv"))

### Subset to Carribean backgrounds (Dominican, Cuban, and Puerto Rican)
# n = 5084
study_data_dropna <- study_data_dropna %>%
  filter(BKGRD1_C6 %in% c(0, 2, 4))

#============== Derive confounder risk groups supervised by CVD ================

# Confounder risk factors
x_mat <- as.matrix(study_data_dropna %>% 
                     select(AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL,
                            DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                            EVER_MI_RELATIVE, EVER_CABG_RELATIVE))
summary(x_mat)
# Binary outcome
y_all <- c(as.numeric(study_data_dropna$CVD_FRAME))  
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
glm_form <- "~ 1"  # no additional confounders
seed <- 1
res_conf_sup <- baysc::swolca(x_mat = x_mat, y_all = y_all, glm_form = glm_form, 
                              sampling_wt = sampling_wt, cluster_id = cluster_id, 
                              stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                              adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                              burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                              save_path = paste0(wd, res_dir, "conf_sup_carrib"))
res_conf_sup$runtime  # 70 mins, 4 classes
item_labels <- c("Age", "Gender", "Smoking", "Alcohol", "PA", "FH_Diab", 
                 "FH_CHD", "FH_Stroke", "Rel_Angina", "Rel_MI", "Rel_CABG")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes (USE NEW REORDER FUNCTION!)
res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(3, 4, 1, 2))
# Use some of the changed functions
plot_pattern_probs(res_conf_sup, item_labels = item_labels, 
                   categ_title = categ_title, y_title = y_title, num_rows = 2)
baysc::plot_pattern_profiles(res_conf_sup, item_labels = item_labels,
                             categ_title = categ_title, y_title = y_title, )
plot_class_dist(res_conf_sup)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_conf_sup$post_MCMC_out$dendrogram)
regr_coefs <- get_regr_coefs(res = res_conf_sup)
plot_regr_coefs(regr_coefs = regr_coefs, res = res_conf_sup)
plot_outcome_probs(res_conf_sup)
# Average posterior probability: 0.92
mean(apply(res_conf_sup$estimates$pred_class_probs, 1, max))


# Framingham score among those with no CVD
study_temp <- study_data_dropna %>% 
  left_join(raw_derv_vars %>% select(ID, FRAME_CVD_RISK_10YR) %>%
              mutate(ID = as.integer(ID)), by = join_by(ID))
summary(study_temp$FRAME_CVD_RISK_10YR[study_temp$CVD_FRAME == 0])
summary(study_temp$FRAME_CVD_RISK_10YR[study_temp$CVD_FRAME == 1])
study_temp$conf_class <- res_conf_sup$estimates$c_all
# Framingham score among those with no CVD, for each class
round(sapply(1:4, function(x) 
  mean(study_temp$FRAME_CVD_RISK_10YR[study_temp$conf_class == x], na.rm = TRUE)), 2)
# svy_design <- survey::svydesign(id = ~cluster_id, weights = ~sampling_wt, 
#                                 strata = ~stratum_id, data = data.frame(cov_df, Class = c_all))


table(study_data_dropna$AGE_C3, res_conf_sup$estimates$c_all)
# Class distribution before and after removing those with CVD
res_conf_sup$estimates$pi_med
prop.table(table(res_conf_sup$estimates$c_all))
prop.table(table(res_conf_sup$estimates$c_all[res_conf_sup$data_vars$y_all == 0]))


### Trace plots
ordering <- 1:res_conf_sup$estimates$K_red
theta_plot <- res_conf_sup$estimates$theta_red[, , ordering, ]
pi_plot <- res_conf_sup$estimates$pi_red[, ordering]
## plot trace of theta
K <- res_conf_sup$estimates$K_red
j <- 1
r <- 1
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Item probability: item ", j, ",theta"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  plot(theta_plot[, j, k, r], type = "l", ylim = c(0,1), xlab = "", ylab = "")
}


## plot trace of pi
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for pi"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  plot(pi_plot[, k], type = "l", ylim = c(0,1), xlab = "", ylab = "")
}


#============== Derive cardiometabolic risk groups =============================

### Supervised by CVD 
# Cardiometabolic indicators risk factors
x_mat <- as.matrix(study_data_dropna %>% 
                     select(bmi_c4, hyp_c3, diab_c4, dyslip_c3))
# , ckd_lab))
summary(x_mat)
# Binary outcome
y_all <- c(as.numeric(study_data_dropna$CVD_FRAME))  
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
glm_form <- "~ 1"  # no additional confounders
seed <- 1
res_cmd_sup <- baysc::swolca(x_mat = x_mat, y_all = y_all, glm_form = glm_form, 
                             sampling_wt = sampling_wt, cluster_id = cluster_id, 
                             stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                             adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                             burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                             save_path = paste0(wd, res_dir, "cmd_sup_carrib"))
res_cmd_sup$runtime  # 45 mins, 5 classes
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(3, 4, 5, 2, 1))
# Use some of the changed functions
plot_pattern_probs(res_cmd_sup, item_labels = item_labels, 
                   categ_title = categ_title, y_title = y_title, num_rows = 1)
baysc::plot_pattern_profiles(res_cmd_sup, item_labels = item_labels, 
                             categ_title = categ_title, y_title = y_title)
plot_class_dist(res_cmd_sup)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_cmd_sup$post_MCMC_out$dendrogram)
regr_coefs <- get_regr_coefs(res = res_cmd_sup)
plot_regr_coefs(regr_coefs = regr_coefs, res = res_cmd_sup)
plot_outcome_probs(res_cmd_sup)
# Average posterior probability: 63
mean(apply(res_cmd_sup$estimates$pred_class_probs, 1, max))

table(study_data_dropna$AGE_C3, res_cmd_sup$estimates$c_all)


# Check framingham risk score
# Framingham score among those with no CVD
study_temp <- study_data_dropna %>% 
  left_join(raw_derv_vars %>% select(ID, FRAME_CVD_RISK_10YR) %>%
              mutate(ID = as.integer(ID)), by = join_by(ID))
summary(study_temp$FRAME_CVD_RISK_10YR[study_temp$CVD_FRAME == 0])
summary(study_temp$FRAME_CVD_RISK_10YR[study_temp$CVD_FRAME == 1])
study_temp$cmd_class <- res_cmd_sup$estimates$c_all
# Framingham score among those with no CVD, for each class
round(sapply(1:3, function(x) 
  mean(study_temp$FRAME_CVD_RISK_10YR[study_temp$cmd_class == x], na.rm = TRUE)), 2)
# svy_design <- survey::svydesign(id = ~cluster_id, weights = ~sampling_wt, 
#                                 strata = ~stratum_id, data = data.frame(cov_df, Class = c_all))


# Class distribution before and after removing those with CVD
res_cmd_sup$estimates$pi_med
prop.table(table(res_cmd_sup$estimates$c_all))
prop.table(table(res_cmd_sup$estimates$c_all[res_cmd_sup$data_vars$y_all == 0]))


### Trace plots
ordering <- 1:res_cmd_sup$estimates$K_red
theta_plot <- res_cmd_sup$estimates$theta_red[, , ordering, ]
pi_plot <- res_cmd_sup$estimates$pi_red[, ordering]
## plot trace of theta
K <- res_cmd_sup$estimates$K_red
j <- 1
r <- 1
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Item probability: item ", j, ",theta"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  plot(theta_plot[, j, k, r], type = "l", ylim = c(0,1), xlab = "", ylab = "")
}


## plot trace of pi
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for pi"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  plot(pi_plot[, k], type = "l", ylim = c(0,1), xlab = "", ylab = "")
}


#=============== Comparing confounder and cardiometabolic risk groups ==========
### NOTE: USE LOCAL VERSION OF reorder_classes
# load(paste0(wd, res_dir, "conf_sup_unaware_carrib_swolca_results.RData"))
# res_conf_sup <- res
# res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(4, 3, 2, 1))
# load(paste0(wd, res_dir, "cmd_sup_unaware_carrib_swolca_results.RData"))
# res_cmd_sup <- res
# res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(2, 3, 1))
# study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_unaware_carrib.csv"))

# Look at grid of confounder and cardiometabolic risk factors
table(res_conf_sup$estimates$c_all, res_cmd_sup$estimates$c_all, useNA = "always")
# By levels of the confounder risk group
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 1])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 2])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 3])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 4])), 2)

## Restrict to those without CVD
# Look at grid of confounder and cardiometabolic risk factors
table(res_conf_sup$estimates$c_all[res_conf_sup$data_vars$y_all == 0], 
      res_cmd_sup$estimates$c_all[res_conf_sup$data_vars$y_all == 0], useNA = "always")
# By levels of the confounder risk group
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$data_vars$y_all == 0 & 
                                      res_conf_sup$estimates$c_all == 1])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$data_vars$y_all == 0 & 
                                      res_conf_sup$estimates$c_all == 2])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$data_vars$y_all == 0 & 
                                      res_conf_sup$estimates$c_all == 3])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$data_vars$y_all == 0 & 
                                      res_conf_sup$estimates$c_all == 4])), 2)

# Create subgroup
study_data_dropna <- study_data_dropna %>%
  mutate(conf_risk = res_conf_sup$estimates$c_all, 
         cmd_risk = res_cmd_sup$estimates$c_all) %>%
  mutate(subgroup = case_when(
    cmd_risk == 1 & conf_risk == 1 ~ 1,  # low-lowF
    cmd_risk == 1 & conf_risk == 2 ~ 2,  # low-lowM
    cmd_risk == 1 & conf_risk == 3 ~ 3,  # low-med
    cmd_risk == 1 & conf_risk == 4 ~ 4,  # low-high
    cmd_risk == 2 & conf_risk == 1 ~ 5,  # med-lowF
    cmd_risk == 2 & conf_risk == 2 ~ 6,  # med-lowM
    cmd_risk == 2 & conf_risk == 3 ~ 7,  # med-med
    cmd_risk == 2 & conf_risk == 4 ~ 8,  # med-high
    cmd_risk == 3 & conf_risk == 1 ~ 9,  # high-lowF
    cmd_risk == 3 & conf_risk == 2 ~ 10,  # high-lowM
    cmd_risk == 3 & conf_risk == 3 ~ 11,  # high-med
    cmd_risk == 3 & conf_risk == 4 ~ 12,  # high-high
    .default = NA
  ), 
  # Subgroup version 2, combining the two confounder low-risk groups
  subgroup2 = case_when(
    cmd_risk == 1 & conf_risk %in% c(1,2) ~ 1,  # low-low
    cmd_risk == 1 & conf_risk == 3 ~ 2,  # low-med
    cmd_risk == 1 & conf_risk == 4 ~ 3,  # low-high
    cmd_risk == 2 & conf_risk %in% c(1,2) ~ 4,  # med-low
    cmd_risk == 2 & conf_risk == 3 ~ 5,  # med-med
    cmd_risk == 2 & conf_risk == 4 ~ 6,  # med-high
    cmd_risk == 3 & conf_risk %in% c(1,2) ~ 7,  # high-low
    cmd_risk == 3 & conf_risk == 3 ~ 8,  # high-med
    cmd_risk == 3 & conf_risk == 4 ~ 9,  # high-high
    .default = NA
  ))
table(study_data_dropna$subgroup)
table(study_data_dropna$subgroup2)


### Restrict to those without CVD
study_data_dropna <- study_data_dropna %>%
  filter(CVD_FRAME == 0)
table(study_data_dropna$subgroup)
table(study_data_dropna$subgroup2)


#================= WRPC model with CVD risk subgroups ==========================

### Run WRPC with subgroup version 1

# Categorical exposure matrix, nxJ, 2879x129
x_mat <- as.matrix(study_data_dropna %>% select(fp2:fp1)) 
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 5000, 
                 burn = 2000, thin = 3, update = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_unaware_carrib"))
res_wrpc$runtime

load(paste0(wd, res_dir, "HCHS_unaware_carrib_adapt.RData"))
res$K_fixed  # 8
hist(res$K_MCMC)
# # Try running fixed sampler with K=8
# seed <- 1
# res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
#                     sampling_wt = sampling_wt, cluster_id = cluster_id, 
#                     stratum_id = stratum_id, run_sampler = "fixed", K_fixed = 8, 
#                     fixed_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
#                     burn = 5000, thin = 5, update = 1000, save_res = TRUE, 
#                     save_path = paste0(wd, res_dir, "HCHS_PR_DR_CVD"))
res_wrpc$runtime  # 3.5 hours
res_wrpc$estimates$pi_med  # 3 classes
res_wrpc$post_MCMC_out$K_med
plot(res_wrpc$post_MCMC_out$dendrogram_global)
mean(apply(res_wrpc$estimates$pred_global_class_probs, 1, max))
plot(res_wrpc$MCMC_out$pi_MCMC[, 1], type = "l")
plot(res_wrpc$post_MCMC_out$pi[, 1], type = "l")
boxplot(res_wrpc$post_MCMC_out$pi[, 1], ylim=c(0,1))
plot(res_wrpc$MCMC_out$theta_local_MCMC[, 1, 1, 4], type = "l")
plot(res_wrpc$MCMC_out$nu_MCMC[, 1, 1], type = "l")
hist(round(res_wrpc$estimates$nu_med, 2), breaks=30)
# Which foods are local
which(res_wrpc$estimates$nu_med < 0.3, arr.ind = TRUE)
plot(res_wrpc$MCMC_out$nu_MCMC[, 12, 6], type = "l")
plot(res_wrpc$MCMC_out$nu_MCMC[, 10, 16], type = "l")
res_wrpc$estimates$theta_local_modes[13, ]
res_wrpc$estimates$theta_global_modes[13, ]
res_wrpc$estimates$theta_local_modes 



### Run WRPC with subgroup version 2

# Categorical exposure matrix, nxJ, 2879x129
x_mat <- as.matrix(study_data_dropna %>% select(fp2:fp1)) 
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup2))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 5000, 
                 burn = 2000, thin = 3, update = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_unaware_carrib_v2"))
res_wrpc$runtime




#=== Try unsupervised risk factor clustering, excluding those aware or w/ CVD ==

# Read in cleaned data
study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_unaware_full.csv"))

### Subset to Carribean backgrounds (Dominican, Cuban, and Puerto Rican)
# n = 5084
study_data_dropna <- study_data_dropna %>%
  filter(BKGRD1_C6 %in% c(0, 2, 4))

### Restrict to those without CVD
# n = 3242
study_data_dropna <- study_data_dropna %>%
  filter(CVD_FRAME == 0)

### Sample size restricted to exclude those with any self-reported conditions 
# n = 2169
study_data_dropna <- study_data_dropna %>%
  filter((hyp_c3 != 3) & (diab_c4 != 3) & (dyslip_c3 != 3))



### Confounder risk factors
x_mat <- as.matrix(study_data_dropna %>% 
                     select(AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL,
                            DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                            EVER_MI_RELATIVE, EVER_CABG_RELATIVE))
summary(x_mat)
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
seed <- 1
res_conf_sup <- baysc::wolca(x_mat = x_mat, 
                              sampling_wt = sampling_wt, cluster_id = cluster_id, 
                              stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                              adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                              burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                              save_path = paste0(wd, res_dir, "conf_carrib_unaware_nocvd"))
res_conf_sup$runtime  # 20 mins, 4 classes
item_labels <- c("Age", "Gender", "Smoking", "Alcohol", "PA", "FH_Diab", 
                 "FH_CHD", "FH_Stroke", "Rel_Angina", "Rel_MI", "Rel_CABG")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes (USE NEW REORDER FUNCTION!)
res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(4, 3, 1, 2))
# Use some of the changed functions
plot_pattern_probs(res_conf_sup, item_labels = item_labels, 
                   categ_title = categ_title, y_title = y_title, num_rows = 2)
baysc::plot_pattern_profiles(res_conf_sup, item_labels = item_labels,
                             categ_title = categ_title, y_title = y_title, )
plot_class_dist(res_conf_sup)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_conf_sup$post_MCMC_out$dendrogram)
# Average posterior probability: 0.79
mean(apply(res_conf_sup$estimates$pred_class_probs, 1, max))
study_temp <- study_data_dropna %>% 
  left_join(raw_derv_vars %>% select(ID, FRAME_CVD_RISK_10YR) %>%
              mutate(ID = as.integer(ID)), by = join_by(ID))
study_temp$conf_class <- res_conf_sup$estimates$c_all
# Framingham score among those with no CVD, for each class
round(sapply(1:4, function(x) 
  mean(study_temp$FRAME_CVD_RISK_10YR[study_temp$conf_class == x], na.rm = TRUE)), 2)



### Cardiometabolic indicators risk factors
x_mat <- as.matrix(study_data_dropna %>% 
                     select(bmi_c4, hyp_c3, diab_c4, dyslip_c3))
summary(x_mat)
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
seed <- 1
res_cmd_sup <- baysc::wolca(x_mat = x_mat, 
                             sampling_wt = sampling_wt, cluster_id = cluster_id, 
                             stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                             adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                             burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                             save_path = paste0(wd, res_dir, "cmd_carrib_unaware_nocvd"))
res_cmd_sup$runtime  # 11 mins, 5 classes
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(3, 1, 2))
# Use some of the changed functions
plot_pattern_probs(res_cmd_sup, item_labels = item_labels, 
                   categ_title = categ_title, y_title = y_title, num_rows = 1)
baysc::plot_pattern_profiles(res_cmd_sup, item_labels = item_labels, 
                             categ_title = categ_title, y_title = y_title)
plot_class_dist(res_cmd_sup)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_cmd_sup$post_MCMC_out$dendrogram)
# Average posterior probability: 69
mean(apply(res_cmd_sup$estimates$pred_class_probs, 1, max))
study_temp <- study_data_dropna %>% 
  left_join(raw_derv_vars %>% select(ID, FRAME_CVD_RISK_10YR) %>%
              mutate(ID = as.integer(ID)), by = join_by(ID))
study_temp$cmd_class <- res_cmd_sup$estimates$c_all
# Framingham score among those with no CVD, for each class
round(sapply(1:3, function(x) 
  mean(study_temp$FRAME_CVD_RISK_10YR[study_temp$cmd_class == x], na.rm = TRUE)), 2)


# Look at grid of confounder and cardiometabolic risk factors
table(res_conf_sup$estimates$c_all, res_cmd_sup$estimates$c_all, useNA = "always")
# By levels of the confounder risk group
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 1])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 2])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 3])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 4])), 2)


# Create subgroup
study_data_dropna <- study_data_dropna %>%
  mutate(conf_risk = res_conf_sup$estimates$c_all, 
         cmd_risk = res_cmd_sup$estimates$c_all) %>%
  mutate(subgroup = case_when(
    cmd_risk == 1 & conf_risk == 1 ~ 1,  # low-lowF
    cmd_risk == 1 & conf_risk == 2 ~ 2,  # low-lowM
    cmd_risk == 1 & conf_risk == 3 ~ 3,  # low-med
    cmd_risk == 1 & conf_risk == 4 ~ 4,  # low-high
    cmd_risk == 2 & conf_risk == 1 ~ 5,  # med-lowF
    cmd_risk == 2 & conf_risk == 2 ~ 6,  # med-lowM
    cmd_risk == 2 & conf_risk == 3 ~ 7,  # med-med
    cmd_risk == 2 & conf_risk == 4 ~ 8,  # med-high
    cmd_risk == 3 & conf_risk == 1 ~ 9,  # high-lowF
    cmd_risk == 3 & conf_risk == 2 ~ 10,  # high-lowM
    cmd_risk == 3 & conf_risk == 3 ~ 11,  # high-med
    cmd_risk == 3 & conf_risk == 4 ~ 12,  # high-high
    .default = NA
  ), 
  # Subgroup version 2, combining the two confounder low-risk groups
  subgroup2 = case_when(
    cmd_risk == 1 & conf_risk %in% c(1,2) ~ 1,  # low-low
    cmd_risk == 1 & conf_risk == 3 ~ 2,  # low-med
    cmd_risk == 1 & conf_risk == 4 ~ 3,  # low-high
    cmd_risk == 2 & conf_risk %in% c(1,2) ~ 4,  # med-low
    cmd_risk == 2 & conf_risk == 3 ~ 5,  # med-med
    cmd_risk == 2 & conf_risk == 4 ~ 6,  # med-high
    cmd_risk == 3 & conf_risk %in% c(1,2) ~ 7,  # high-low
    cmd_risk == 3 & conf_risk == 3 ~ 8,  # high-med
    cmd_risk == 3 & conf_risk == 4 ~ 9,  # high-high
    .default = NA
  ))
table(study_data_dropna$subgroup)
table(study_data_dropna$subgroup2)
