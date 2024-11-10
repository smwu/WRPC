#=========================
# Analyzing HCHS data
# Author: Stephanie Wu
# Date Updated: 2024/11/06
#=========================

rm(list = ls())

library(tidyverse)  # data wrangling
library(baysc)      # bayesian survey clustering
library(haven)      # read sas file
library(naniar)     # missingness
library(sfsmisc)    # mult.fig

wd <- "~/Documents/Github/WRPC/"
code_dir <- "Model_Code/"
data_dir <- "Application/HCHS_Data/"
res_dir <- "Results/"

# Read in cleaned data, subsetted to PR and complete data
study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data.csv"))

# Source WRPC functions
source(paste0(wd, code_dir, "wrpc.R"))
source(paste0(wd, code_dir, "wrpc_mcmc_fns.R"))
Rcpp::sourceCpp(paste0(wd, code_dir, "wrpc_mcmc.cpp"))


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
                     save_path = paste0(wd, res_dir, "conf_sup"))
res_conf_sup$runtime  # 17 mins
item_labels <- c("Age", "Gender", "Smoking", "Alcohol", "PA", "FH_Diab", 
                 "FH_CHD", "FH_Stroke", "Rel_Angina", "Rel_MI", "Rel_CABG")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(2, 1, 3))
# Use some of the changed functions
plot_pattern_probs(res_conf_sup, item_labels = item_labels, 
                   categ_title = categ_title, y_title = y_title, num_rows = 2)
baysc::plot_pattern_profiles(res_conf_sup, item_labels = item_labels,
                             categ_title = categ_title, y_title = y_title, )
baysc::plot_class_dist(res_conf_sup)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_conf_sup$post_MCMC_out$dendrogram)
regr_coefs <- get_regr_coefs(res = res_conf_sup)
plot_regr_coefs(regr_coefs = regr_coefs, res = res_conf_sup)
plot_outcome_probs(res_conf_sup)
# Average posterior probability 
mean(apply(res_conf_sup$estimates$pred_class_probs, 1, max))

# # For no_recode version
# res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(3, 2, 1, 4))


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
                     select(bmi_c3, hyp_lab_med, diab_lab_med, dyslip_lab_med))
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
                          save_path = paste0(wd, res_dir, "cmd_sup_no_ckd"))
res_cmd_sup$runtime  # 10 mins
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia", "CKD")
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(2, 3, 1, 4))
# Use some of the changed functions
plot_pattern_probs(res_cmd_sup, item_labels = item_labels, 
                   categ_title = categ_title, y_title = y_title, num_rows = 1)
baysc::plot_pattern_profiles(res_cmd_sup, item_labels = item_labels, 
                             categ_title = categ_title, y_title = y_title)
baysc::plot_class_dist(res_cmd_sup)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_cmd_sup$post_MCMC_out$dendrogram)
regr_coefs <- get_regr_coefs(res = res_cmd_sup)
plot_regr_coefs(regr_coefs = regr_coefs, res = res_cmd_sup)
plot_outcome_probs(res_cmd_sup)
# Average posterior probability 
mean(apply(res_cmd_sup$estimates$pred_class_probs, 1, max))

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
load(paste0(wd, res_dir, "conf_sup_swolca_results.RData"))
res_conf_sup <- res
res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(2, 1, 3))
load(paste0(wd, res_dir, "cmd_sup_no_ckd_swolca_results.RData"))
res_cmd_sup <- res
res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(2, 3, 1, 4))

# Look at grid of confounder and cardiometabolic risk factors
table(res_conf_sup$estimates$c_all, res_cmd_sup$estimates$c_all, useNA = "always")
# By levels of the confounder risk group
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 1])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 2])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 3])), 2)

# Create subgroup
study_data_dropna <- study_data_dropna %>%
  mutate(conf_risk = res_conf_sup$estimates$c_all, 
         cmd_risk = res_cmd_sup$estimates$c_all) %>%
  mutate(subgroup = case_when(
    cmd_risk == 1 & conf_risk == 1 ~ 1,  # low-low
    cmd_risk == 1 & conf_risk == 2 ~ 2,  # low-med
    cmd_risk == 1 & conf_risk == 3 ~ 3,  # low-high
    cmd_risk == 2 & conf_risk == 1 ~ 4,  # low-dys-low
    cmd_risk == 2 & conf_risk == 2 ~ 5,  # low-dys-med
    cmd_risk == 2 & conf_risk == 3 ~ 6,  # low-dys-high
    cmd_risk == 3 & conf_risk == 1 ~ 7,  # med-low
    cmd_risk == 3 & conf_risk == 2 ~ 8,  # med-med
    cmd_risk == 3 & conf_risk == 3 ~ 9,  # med-high
    cmd_risk == 4 & conf_risk == 1 ~ 10,  # high-low
    cmd_risk == 4 & conf_risk == 2 ~ 11,  # high-med
    cmd_risk == 4 & conf_risk == 3 ~ 12,  # high-high
    .default = NA
  ))
table(study_data_dropna$subgroup)


### Try version with 4 confounder groups
res_conf_sup2 <- load(paste0(wd, res_dir, "conf_sup_no_recode_swolca_results.RData"))
res_conf_sup2 <- reorder_classes(res = res_conf_sup2, new_order = c(3, 2, 1, 4))
study_data_dropna <- study_data_dropna %>%
  mutate(conf_risk2 = res_conf_sup2$estimates$c_all)
mutate(subgroup2 = case_when(
  cmd_risk == 1 & conf_risk2 == 1 ~ 1,  # low-low
  cmd_risk == 1 & conf_risk2 == 2 ~ 2,  # low-med
  cmd_risk == 1 & conf_risk2 == 3 ~ 3,  # low-high
  cmd_risk == 1 & conf_risk2 == 4 ~ 4,  # low-high_fh
  cmd_risk == 2 & conf_risk2 == 1 ~ 5,  # low-dys-low
  cmd_risk == 2 & conf_risk2 == 2 ~ 6,  # low-dys-med
  cmd_risk == 2 & conf_risk2 == 3 ~ 7,  # low-dys-high
  cmd_risk == 2 & conf_risk2 == 4 ~ 8,  # low-dys-high_fh
  cmd_risk == 3 & conf_risk2 == 1 ~ 9,  # med-low
  cmd_risk == 3 & conf_risk2 == 2 ~ 10,  # med-med
  cmd_risk == 3 & conf_risk2 == 3 ~ 11,  # med-high
  cmd_risk == 3 & conf_risk2 == 4 ~ 12,  # med-high_fh
  cmd_risk == 4 & conf_risk2 == 1 ~ 13,  # high-low
  cmd_risk == 4 & conf_risk2 == 2 ~ 14,  # high-med
  cmd_risk == 4 & conf_risk2 == 3 ~ 15,  # high-high
  cmd_risk == 4 & conf_risk2 == 4 ~ 16,  # high-high_fh
  .default = NA
))

#================= WRPC model with CVD risk subgroups ==========================

### Run WRPC
# Categorical exposure matrix, nxJ, 1988x129
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
res_cmd_sup <- wrpc(x_mat = x_mat, h_all = h_all, 
                    sampling_wt = sampling_wt, cluster_id = cluster_id, 
                    stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                    adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                    burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                    save_path = paste0(wd, res_dir, "HCHS_PR_CVD"))


# # Try version with 4 confounder groups
# # Subgroup variable, nx1
# h_all <- c(as.numerc(study_data_dropna$subgroup2))
# 
# seed <- 1
# res_cmd_sup <- wrpc(x_mat = x_mat, y_all = y_all, glm_form = glm_form, 
#                     sampling_wt = sampling_wt, cluster_id = cluster_id, 
#                     stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
#                     adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
#                     burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
#                     save_path = paste0(wd, res_dir, "HCHS_PR_CVD_4conf"))



#================= OLD CODE ====================================================
#================= Derive cardiometabolic risk groups: unsupervised ============
### Unsupervised
# Cardiometabolic indicators risk factors
x_mat <- as.matrix(study_data_dropna %>% 
                     select(bmi_c3, hyp_lab_med, diab_lab_med, dyslip_lab_med, 
                            ckd_lab))
summary(x_mat)
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
seed <- 20241107
res_cmd <- baysc::wolca(x_mat = x_mat, sampling_wt = sampling_wt, 
                        cluster_id = cluster_id, stratum_id = stratum_id, 
                        run_sampler = "both", K_max = 30, adapt_seed = seed, 
                        class_cutoff = 0.05, n_runs = 20000, burn = 10000, 
                        thin = 5, update = 1000, save_res = TRUE, 
                        save_path = paste0(wd, res_dir, "cmd"))
res_cmd$runtime  # 10 mins
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia", "CKD")
# Use some of the changed functions
plot_pattern_probs(res_cmd, item_labels = item_labels, num_rows = 1)
baysc::plot_pattern_profiles(res_cmd, item_labels = item_labels)
baysc::plot_class_dist(res_cmd)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_cmd$post_MCMC_out$dendrogram)
# Average posterior probability 
mean(apply(res_cmd$estimates$pred_class_probs, 1, max))


### Trace plots
ordering <- 1:res_cmd$estimates$K_red
theta_plot <- res_cmd$estimates$theta_red[, , ordering, ]
pi_plot <- res_cmd$estimates$pi_red[, ordering]
## plot trace of theta
K <- res_cmd$estimates$K_red
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


#================= WRPC model with diabetes ====================================

# Source functions
source(paste0(wd, code_dir, "wrpc.R"))
source(paste0(wd, code_dir, "wrpc_mcmc_fns.R"))
Rcpp::sourceCpp(paste0(wd, code_dir, "wrpc_mcmc.cpp"))


# Overfitted K
K_max <- 20

# Run WRPC
x_mat <- as.matrix(study_data %>% select(fp2:fp1)) # Categorical exposure matrix, nxJ
stratum_id <- study_data$STRAT      # Stratum indicators, nx1
cluster_id <- study_data$PSU_ID                   # Cluster indicators, nx1
sampling_wt <- study_data$WEIGHT_FINAL_EXPANDED   # Survey sampling weights, nx1
h_all <- study_data$subgroup           # Subpopulation indicators, nx1

### Age and gender
n_runs <- 2000
burn <- 1000
thin <- 5
update <- 100
switch <- 50
class_cutoff_global <- 0.05
seed <- 20241105
res_both <- wrpc(K_max = K_max, x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, run_sampler = "both", 
                 adapt_seed = seed, n_runs = n_runs, burn = burn, thin = thin, 
                 update = update, save_res = TRUE,
                 save_path = paste0(wd, res_dir, "HCHS_age_gender"))
### Age and gender
n_runs <- 10000
burn <- 5000
thin <- 5
update <- 100
switch <- 50
class_cutoff_global <- 0.05
seed <- 20241105
res_both <- wrpc(K_max = K_max, x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, run_sampler = "both", 
                 adapt_seed = seed, n_runs = n_runs, burn = burn, thin = thin, 
                 update = update, save_res = TRUE,
                 save_path = paste0(wd, res_dir, "HCHS_age_gender10000"))


### Test version
n_runs <- 200
burn <- 0
thin <- 1
update <- n_runs / 10
update <- 10
switch <- 50
class_cutoff_global <- 0.05

seed <- 20241102
res_both <- wrpc(K_max = K_max, x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, run_sampler = "both", 
                 adapt_seed = seed, n_runs = n_runs, burn = burn, thin = thin, 
                 update = update, save_res = TRUE,
                 save_path = paste0(wd, res_dir, "test_HCHS_diabetes"))
res_fixed <- wrpc(K_max = K_max, x_mat = x_mat, h_all = h_all, 
                  sampling_wt = sampling_wt, run_sampler = "fixed", K_fixed = 7, 
                  adapt_seed = seed, n_runs = n_runs, burn = burn, thin = thin, 
                  update = update, save_res = TRUE,
                  save_path = paste0(wd, res_dir, "test_HCHS_diabetes"))


K <- res_both$estimates$K_red
theta_global_plot <- res_both$estimates$theta_global_red[, , , ]
theta_local_plot <- res_both$estimates$theta_local_red
# pi_plot <- res_both$estimates$pi_red[, ordering]
# nu_plot <- res_both$estimates$nu_red
# c_plot <- ifelse(res_both$estimates$c_all == ordering[1], 1, 
#                  ifelse(res_both$estimates$c_all == ordering[2], 2, 3))
# g_mat <- res_both$estimates$g_mat
# K <- res_both$estimates$K_red
# theta_global_modes <- res_both$estimates$theta_global_modes
# theta_local_modes <- res_both$estimates$theta_local_modes
# nu_red <- res_both$estimates$nu_red
# H <- res_both$data_vars$H

## plot trace of theta_global
plot(res_both$estimates$theta_global_red[, 1,1,1], type = "l")
plot(res_both$estimates$theta_local_red[, 1,1,1], type = "l")

theta_global_modes

theta_local_modes

round(res_both$estimates$nu_med, 3)