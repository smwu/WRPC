#=================================
# Analyzing HCHS data w/ PR and DR
# Author: Stephanie Wu
# Date Updated: 2024/11/11
#=================================

#==================== Preparing HCHS data w/ PR and DR =========================

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
# Read in raw FPQ data
raw_fpq_vars <- read_sas(paste0(wd, data_dir, "fpe_inv4.sas7bdat"))
# Read in cleaned FPQ data
sol_data <- read.csv(paste0(wd, data_dir, "FPQ129_CVD07FEB2023.csv"))

# Source WRPC functions
source(paste0(wd, code_dir, "wrpc.R"))
source(paste0(wd, code_dir, "wrpc_mcmc_fns.R"))
Rcpp::sourceCpp(paste0(wd, code_dir, "wrpc_mcmc.cpp"))

#================= FPQ cleaning ================================================

# Select diet variables
fpq_data <- sol_data %>% 
  select(ID, STRAT, PSU_ID, CENTER, BKGRD1_C6, fp2:fp1)

### Check site x background numbers
# Site: B=Bronx, C=Chicago, M=Miami, S=San Diego
# Background: 0=Dominican, 1=Central or South American, 2=Cuban, 3=Mexican, 
# 4=Puerto Rican, 5=More than one/Other heritage 
table(fpq_data$CENTER, fpq_data$BKGRD1_C6)


### Subset to Puerto Rican background and complete data: total 3218
fpq_PR_DR <- fpq_data %>% 
  filter(BKGRD1_C6 %in% c(0, 4)) %>%
  drop_na()

# Merge in derived variables, which include demographic and outcome variables
# Combined study data: 3218 x 393
study_data <- raw_derv_vars %>% 
  mutate(ID = as.integer(ID)) %>%
  right_join(fpq_PR_DR %>% 
               select(ID, fp2:fp1) %>% 
               mutate(across(fp2:fp1, as.numeric)), by = join_by(ID))

# Find 5-level foods
which(apply(fpq_PR_DR[, -c(1:5)], 2, 
            function(x) max(x, na.rm = TRUE) == 5))
# Convert 5-level foods to 4 levels
study_data <- study_data %>%
  mutate_at(c("fp2", "fp116", "fp7", "fp10", "fp11", "fp12", "fp13", "fp17", 
              "fp53", "fp101"), ~ifelse(.x == 5, 4, .x))
# Check no 5-level foods
apply(study_data %>% select(fp2:fp1), 2, function(x) max(x, na.rm = TRUE))


#==================== Confounders cleaning ========================================

# Look at distribution of various variables
vars_summ <- study_data %>% 
  select(AGE, GENDER, WAIST_HIP,
         CIGARETTE_PACK_YEARS_C3, # 1: 0, 2: 0-10, 3: >10. 37 NAs
         CIGARETTE_USE, # 1: never, 2: former, 3: current. 5 NAs
         ALCOHOL_USE_LEVEL, # 1: no current, 2: low level, 3: high level. 4 NAs
         ALCOHOL_USE, # 1: never, 2: former, 3: current. 3 NAs
         GPAQ_LEVEL, # Physical activity. 1:high, 2: moderate, 3: low. 11 NAs
         DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE, # 0: No, 1: Yes. ~44 NAs
         EVER_MI_RELATIVE, EVER_CABG_RELATIVE,  # 0: No, 1: Yes. ~24 NAs
         CKD, # 1:normal, 2:mild, 3:moderate, 4:severe, 5:end-stage. 15 NAs
         AHI_GE15, # Sleep apnea. 1=yes. # 412 NAs
         COPD_BY_BD,  # 1=yes. 2944 NAs
         MED_LLD, # lipid lowering drugs/antihyperlipidemics. 1=yes. 64 NAs
         MED_ANTIDIAB, # Antidiabetics. 1=yes. 64 NAs
         MED_ANTIHYPERT # Antihypertensives. 1=yes. 64 NAs
  )
dim(vars_summ)  # n=3218
summary(vars_summ$AGE)
summary(vars_summ$WAIST_HIP)
apply(vars_summ[, -1], 2, function(x) 
  round(prop.table(table(x, useNA = "always")), 2))
apply(vars_summ[, -1], 2, function(x) table(x, useNA = "always"))

# CVD for supervising
# Composite CVD including CHD (coronary death, myocardial infarction, coronary 
# insufficiency, and angina), cerebrovascular events (including ischemic stroke, 
# hemorrhagic stroke, and transient ischemic attack), peripheral artery disease 
# (intermittent claudication), and heart failure.
table(study_data$CVD_FRAME, useNA = "always")


# Check missingness
# Note: FPQ is already subsetted to complete cases
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             AGE, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL, 
                             DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                             EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME))
naniar::vis_miss(study_data %>% 
                   select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                          AGE, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL, 
                          DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                          EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME))

# Variable conversions
study_data_conf <- study_data %>%
  mutate(GENDER = case_match(GENDER, 
                             "M" ~ 1, 
                             "F" ~ 3,
                             .default = NA),
         # Age categories
         AGE_C3 = case_when(
           AGE < 45 ~ 1, 
           (AGE >= 45) & (AGE < 65) ~ 2,
           AGE >= 65 ~ 3, 
           .default = NA)) %>%
  # Convert levels to 1, 2 for binary for swolca function
  mutate_at(c("DIAB_FAMHIST", "FH_CHD", "FH_STROKE", "EVER_ANGINA_RELATIVE", 
              "EVER_MI_RELATIVE", "EVER_CABG_RELATIVE"), function(x) 
                ifelse(x == 0, 1, ifelse(x == 1, 3, NA)))


#================= Cardiometabolic risk factor cleaning ========================


# Clinical cardiometabolic indicators
study_data_conf_cmd <- study_data_conf %>%
  mutate(
    bmi_c3 = case_when(
      BMIGRP_C4 %in% c(1, 2) ~ 1, # normal or underweight BMI<25
      BMIGRP_C4 == 3 ~ 2, # overweight 25<=BMI<30
      BMIGRP_C4 == 4 ~ 3, # obese BMI >= 30
      .default = NA),   
    hyp_lab_med = case_when(
      HYPERTENSION == 0 ~ 1, # none
      HYPERTENSION == 1 ~ 3, # 140/90 or scanned meds
      .default = NA),
    diab_lab_med = DIABETES2, # non, pre, diabetic. FBS>=>200 or post-OGTT>=200 or A1C>=6.5% or scanned meds
    dyslip_lab_med = case_when(
      DYSLIPIDEMIA == 0 ~ 1, # none
      DYSLIPIDEMIA == 1 | MED_LLD == 1 ~ 3, # LDL>=160 or HDL<=40 or trigl>=200 or scanned meds
      .default = NA),
    ckd_lab = case_when(
      CKD == 1 ~ 1, # normal, GFR>=90
      CKD == 2 ~ 2, # mild, 60<=GFR<90
      CKD %in% c(3, 4, 5) ~ 3, # moderate, severe, end-stage, GFR<60
    ))  

# Check missingness
naniar::gg_miss_var(study_data_conf_cmd %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL, 
                             DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                             EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
                             bmi_c3, hyp_lab_med, diab_lab_med, dyslip_lab_med, 
                             ckd_lab))
naniar::vis_miss(study_data_conf_cmd %>% 
                   select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                          AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL, 
                          DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                          EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
                          bmi_c3, hyp_lab_med, diab_lab_med, dyslip_lab_med, 
                          ckd_lab))

# Sample size after dropping all missingness: n = 3116
study_data_dropna <- study_data_conf_cmd %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, CENTER, BKGRD1_C6,
         AGE_C3, GENDER, CIGARETTE_USE, 
         ALCOHOL_USE, GPAQ_LEVEL, DIAB_FAMHIST, FH_CHD, FH_STROKE, 
         EVER_ANGINA_RELATIVE, EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
         bmi_c3, hyp_lab_med, diab_lab_med, dyslip_lab_med, ckd_lab, fp2:fp1) %>% 
  drop_na()

# Save cleaned data
write.csv(study_data_dropna, file = paste0(wd, data_dir, "cleaned_data_PR_DR.csv"),
          row.names = FALSE)



#============== Derive confounder risk groups supervised by CVD ================
# study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_PR_DR.csv"))

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
                              save_path = paste0(wd, res_dir, "conf_sup_PR_DR"))
res_conf_sup$runtime  # 40 mins, 3 classes
item_labels <- c("Age", "Gender", "Smoking", "Alcohol", "PA", "FH_Diab", 
                 "FH_CHD", "FH_Stroke", "Rel_Angina", "Rel_MI", "Rel_CABG")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(3, 2, 1))
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
# Average posterior probability: 0.96
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

# Try no thinning and Kmax=20
seed <- 1
res_conf_sup <- baysc::swolca(x_mat = x_mat, y_all = y_all, glm_form = glm_form, 
                              sampling_wt = sampling_wt, cluster_id = cluster_id, 
                              stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                              adapt_seed = seed, class_cutoff = 0.05, n_runs = 10000, 
                              burn = 5000, thin = 1, update = 1000, save_res = TRUE, 
                              save_path = paste0(wd, res_dir, "conf_sup_PR_DR_nothin"))

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
                             save_path = paste0(wd, res_dir, "cmd_sup_no_ckd_PR_DR"))
res_cmd_sup$runtime  # 23 mins, 4 classes
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(3, 4, 2, 1))
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
# Average posterior probability: 73
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


# Try no thinning and Kmax=20
seed <- 1
res_cmd_sup <- baysc::swolca(x_mat = x_mat, y_all = y_all, glm_form = glm_form, 
                              sampling_wt = sampling_wt, cluster_id = cluster_id, 
                              stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                              adapt_seed = seed, class_cutoff = 0.05, n_runs = 10000, 
                              burn = 5000, thin = 1, update = 1000, save_res = TRUE, 
                              save_path = paste0(wd, res_dir, "cmd_sup_PR_DR_nothin"))
res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(4, 3, 2, 1))



### Unsupervised 
# Cardiometabolic indicators risk factors
x_mat <- as.matrix(study_data_dropna %>% 
                     select(bmi_c3, hyp_lab_med, diab_lab_med, dyslip_lab_med))
# , ckd_lab))
summary(x_mat)
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
seed <- 1
res_cmd <- baysc::wolca(x_mat = x_mat, 
                             sampling_wt = sampling_wt, cluster_id = cluster_id, 
                             stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                             adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                             burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                             save_path = paste0(wd, res_dir, "cmd_no_ckd_PR_DR"))
res_cmd$runtime  # 17 mins, 5 classes
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
# res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(3, 4, 2, 1))
# Use some of the changed functions
plot_pattern_probs(res_cmd, item_labels = item_labels, 
                   categ_title = categ_title, y_title = y_title, num_rows = 1)
baysc::plot_pattern_profiles(res_cmd, item_labels = item_labels, 
                             categ_title = categ_title, y_title = y_title)
baysc::plot_class_dist(res_cmd)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_cmd$post_MCMC_out$dendrogram)
# Average posterior probability: 80
mean(apply(res_cmd$estimates$pred_class_probs, 1, max))

#=============== Comparing confounder and cardiometabolic risk groups ==========
### NOTE: USE LOCAL VERSION OF reorder_classes
# load(paste0(wd, res_dir, "conf_sup_PR_DR_swolca_results.RData"))
# res_conf_sup <- res
# res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(3, 2, 1))
# load(paste0(wd, res_dir, "cmd_sup_no_ckd_PR_DR_swolca_results.RData"))
# res_cmd_sup <- res
# res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(3, 4, 2, 1))
# study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_PR_DR.csv"))

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


#================= WRPC model with CVD risk subgroups ==========================

### Run WRPC
# Categorical exposure matrix, nxJ, 3116x129
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
                    adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                    burn = 5000, thin = 5, update = 1000, save_res = TRUE, 
                    save_path = paste0(wd, res_dir, "HCHS_PR_DR_CVD"))
load(paste0(wd, res_dir, "HCHS_PR_DR_CVD_wrpc_adapt.RData"))
res$K_fixed  # 8
hist(res$K_MCMC)
# Try running fixed sampler with K=8
seed <- 1
res_cmd_sup <- wrpc(x_mat = x_mat, h_all = h_all, 
                    sampling_wt = sampling_wt, cluster_id = cluster_id, 
                    stratum_id = stratum_id, run_sampler = "fixed", K_fixed = 8, 
                    fixed_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                    burn = 5000, thin = 5, update = 1000, save_res = TRUE, 
                    save_path = paste0(wd, res_dir, "HCHS_PR_DR_CVD"))
res_wrpc <- res_cmd_sup
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
  
