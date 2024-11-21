#====================================================
# Analyzing HCHS data for Carribean background
# Step 1: Deriving CVD risk groups using CVD and 
# including self-report in cardiometabolic definitions
# Author: Stephanie Wu
# Date Updated: 2024/11/17
#====================================================

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
# Read in cleaned FPQ data w/ 49 categories
fpq_49categ <- read.csv(paste0(wd, data_dir, "fpq_49categ_2024Nov19.csv"))

# Source WRPC functions
source(paste0(wd, code_dir, "wrpc.R"))
source(paste0(wd, code_dir, "wrpc_mcmc_fns.R"))
Rcpp::sourceCpp(paste0(wd, code_dir, "wrpc_mcmc.cpp"))

# Prevalent CVD 
prop.table(table(raw_derv_vars$PRECHD_ANGINA, useNA = "always"))
# Ever CVD
prop.table(table(raw_derv_vars$CVD_FRAME, useNA = "always"))

#================= FPQ cleaning ================================================

# Select diet variables
fpq_data <- fpq_49categ %>% 
  select(ID, citrus_juice:soup_oth)

# Drop NAs: n=12761
fpq_dropna <- fpq_data %>% drop_na()

# Merge in derived variables, which include demographic and outcome variables.
# Merge in self-report variables
# Combined study data: 12761 x 441
study_data <- raw_derv_vars %>% 
  mutate(ID = as.integer(ID)) %>%
  left_join(raw_mhea_vars %>% mutate(ID = as.integer(ID)), by = join_by(ID)) %>%
  right_join(fpq_dropna %>% 
               mutate(across(citrus_juice:soup_oth, as.numeric)), by = join_by(ID)) 

# Check all 5-level foods
apply(study_data %>% select(citrus_juice:soup_oth), 2, 
      function(x) max(x, na.rm = TRUE))

# 
# # Select diet variables
# fpq_data <- sol_data %>% 
#   select(ID, STRAT, PSU_ID, CENTER, BKGRD1_C6, fp2:fp1)
# 
# ### Check site x background numbers
# # Site: B=Bronx, C=Chicago, M=Miami, S=San Diego
# # Background: 0=Dominican, 1=Central or South American, 2=Cuban, 3=Mexican, 
# # 4=Puerto Rican, 5=More than one/Other heritage 
# table(fpq_data$CENTER, fpq_data$BKGRD1_C6)
# 
# # Drop NAs: n=12761
# fpq_dropna <- fpq_data %>% drop_na()
# 
# # Merge in derived variables, which include demographic and outcome variables.
# # Merge in self-report variables
# # Combined study data: 12761 x 521
# study_data <- raw_derv_vars %>% 
#   mutate(ID = as.integer(ID)) %>%
#   left_join(raw_mhea_vars %>% mutate(ID = as.integer(ID)), by = join_by(ID)) %>%
#   right_join(fpq_dropna %>% 
#                select(ID, fp2:fp1) %>% 
#                mutate(across(fp2:fp1, as.numeric)), by = join_by(ID)) 
# 
# # Find 5-level foods
# which(apply(fpq_dropna[, -c(1:5)], 2, 
#             function(x) max(x, na.rm = TRUE) == 5))
# # Convert 5-level foods to 4 levels
# study_data <- study_data %>%
#   mutate_at(c("fp2", "fp116", "fp7", "fp10", "fp11", "fp12", "fp13", "fp17", 
#               "fp53", "fp101"), ~ifelse(.x == 5, 4, .x))
# # Check no 5-level foods
# apply(study_data %>% select(fp2:fp1), 2, function(x) max(x, na.rm = TRUE))


#==================== Confounders cleaning ========================================

# Check missingness
# Note: FPQ is already subsetted to complete cases
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             AGE, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL, 
                             DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                             EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME))

# Variable conversions
study_data_conf <- study_data %>%
  mutate(GENDER = case_match(GENDER, 
                             "F" ~ 1, 
                             "M" ~ 3,
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

### Definitions for hypertension
# 3637 lab-based or scanned meds
table(study_data$HYPERTENSION, useNA = "always") 
# 4013 lab-based or scanned meds or self-report meds
table(study_data$HYPERTENSION2, useNA = "always") 
# 3044 of HYPERTENSION2 are aware
table(study_data$HYPERT_AWARENESS, useNA = "always") 
# 4 categories using HYPERTENSION2, MUEA33D, MED_ANTIHYPERT variables
table(study_data$HYPERTENSION_C4, useNA = "always") 
# 4018 self-report hypertension (some are in non- or pre- group for HYPERTENSION2)
table(study_data$MHEA1, useNA = "always") 
# Hypertension by awareness: unaware exists for both treated (274) and 
# untreated (686) hypertension
table(study_data$HYPERTENSION_C4, study_data$MHEA1, useNA = "always")
# WE DEFINE: 1) none, 2) pre-, 3) treated (scanned or self-report meds), 
# 4) untreated aware (no scanned or self-report meds, but self-report condition),
# 5) untreated unaware (lab-based or scanned meds, but no self-report meds or condition)

### Definitions for diabetes
# 2239 lab-based 
table(study_data$DIABETES_LAB, useNA = "always") 
# 2570 lab-based or scanned meds
table(study_data$DIABETES2, useNA = "always") 
# 2685 lab-based or scanned meds or self-report diabetes
table(study_data$DIABETES3, useNA = "always") 
# 4 categories using DIABETES3, MUEA33C (self-report meds), MED_ANTIDIAB variables
table(study_data$DIABETES_C4, useNA = "always") 
# 1982 self-report diabetes (some are in non- or pre- group for DIABETES3)
table(study_data$MHEA16, useNA = "always") 
# Diabetes by awareness: unaware exists for both treated (22) and 
# untreated (842) diabetes
table(study_data$DIABETES_C4, study_data$MHEA16, useNA = "always")
# WE DEFINE: 1) none, 2) pre-, 3) treated (scanned or self-report meds), 
# 4) untreated aware (no scanned or self-report meds, but self-report condition),
# 5) untreated unaware (lab-based or scanned meds, but no self-report meds or condition)

### Definitions for dyslipidemia
# 4946 lab-based 
table(study_data$DYSLIPIDEMIA, useNA = "always") 
# WE DON'T HAVE DYSLIPIDEMIA_C3
# # 3 categories using DYSLIPIDEMIA, MUEA33E (self-report meds), MED_LLD variables
# table(study_data$DYSLIPIDEMIA_C3, useNA = "always") 
# 1647 report antihyperlipidemic medication (some are in the none group for DYSLIPIDEMIA)
table(study_data$DYSLIPIDEMIA, study_data$MED_LLD, useNA = "always")
# 4703 self-report high cholesterol (some are in the none group for DYSLIPIDEMIA)
table(study_data$MHEA2, useNA = "always") 
# Dyslipidemia by awareness: unaware (2545) 
table(study_data$DYSLIPIDEMIA, study_data$MHEA2, useNA = "always")
# WE DEFINE: 1) none, 2) aware (lab-based or scanned med, & self-report condition), 
# 3) unaware (lab-based or scanned meds, but no self-report condition)


# Clinical cardiometabolic indicators
study_data_conf_cmd <- study_data_conf %>%
  mutate(
    bmi_c3 = case_when(
      BMIGRP_C6 %in% c(1, 2) ~ 1, # normal or underweight BMI<25
      BMIGRP_C6 == 3 ~ 2, # overweight 25<=BMI<30
      BMIGRP_C6 %in% c(4, 5, 6) ~ 3, # obese BMI >= 30
      .default = NA),  
    # hypertension self-report: 1 = yes
    hyp_sr = MHEA1,
    # hypertension including self-report, lab-based, and scanned meds
    hyp_lab_med_sr = case_when(
      HYPERTENSION == 0 ~ 1, # none
      HYPERTENSION == 1 | MHEA1 == 1 ~ 3, # hypertension (lab-based or scanned meds or self-report)
      .default = NA),
    # diabetes self-report: 1 = yes
    diab_sr = MHEA16,
    # diabetes including self-report, lab-based, and scanned meds
    diab_lab_med_sr = case_when(
      DIABETES2 == 1 ~ 1, # none 
      DIABETES2 == 2 ~ 2, # pre-diabetic (fasting time > 8 AND fasting glucose in range 100-125 mg/dL) or (post-OGTT glucose in range 140-199 mg/dL) or (5.7%≤A1C<6.5%)
      DIABETES2 == 3 | MHEA16 == 1 ~ 3, # diabetes (FBS>=>200 or post-OGTT>=200 or A1C>=6.5% or scanned meds or self-report)
      .default = NA),
    # dyslipidemia self-report: 1 = yes
    dyslip_sr = MHEA2,
    # dyslipidemia including self-report, lab-based, and scanned meds
    dyslip_lab_med_sr = case_when(
      DYSLIPIDEMIA == 0 ~ 1, # none
      DYSLIPIDEMIA == 1 | MED_LLD == 1 | MHEA2 == 1 ~ 3, # high chol (LDL>=160 or HDL<=40 or trigl>=200 or scanned meds, & self-report)
      .default = NA)
  )  

# Check missingness
naniar::gg_miss_var(study_data_conf_cmd %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL, 
                             DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                             EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
                             bmi_c3, hyp_sr, hyp_lab_med_sr, diab_sr, diab_lab_med_sr, 
                             dyslip_sr, dyslip_lab_med_sr))


#==================== Preparing HCHS data w/ PR, DR, and Cuban =================

# Sample size after dropping all missingness: n = 12489
study_data_dropna <- study_data_conf_cmd %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, CENTER, BKGRD1_C6,
         AGE_C3, GENDER, CIGARETTE_USE, 
         ALCOHOL_USE, GPAQ_LEVEL, DIAB_FAMHIST, FH_CHD, FH_STROKE, 
         EVER_ANGINA_RELATIVE, EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
         bmi_c3, hyp_sr, hyp_lab_med_sr, diab_sr, diab_lab_med_sr, 
         dyslip_sr, dyslip_lab_med_sr, citrus_juice:soup_oth) %>% 
  drop_na()

# Framingham CVD (ever) levels: 3616 w/ ever CVD (29%)
# Composite CVD including CHD (coronary death, myocardial infarction, coronary 
# insufficiency, and angina), cerebrovascular events (including ischemic stroke, 
# hemorrhagic stroke, and transient ischemic attack), peripheral artery disease 
# (intermittent claudication), and heart failure.
table(study_data_dropna$CVD_FRAME, useNA = "always")
prop.table(table(study_data_dropna$CVD_FRAME, useNA = "always"))

### Subset to Carribean backgrounds (Dominican, Cuban, and Puerto Rican)
# n = 5086
study_data_dropna <- study_data_dropna %>%
  filter(BKGRD1_C6 %in% c(0, 2, 4))

table(study_data_dropna$CVD_FRAME, useNA = "always")  # 1844
prop.table(table(study_data_dropna$CVD_FRAME, useNA = "always")) # 36%

# # Save cleaned data
# write.csv(study_data_dropna, file = paste0(wd, data_dir, "cleaned_data_carrib.csv"),
#           row.names = FALSE)


#============== Derive confounder risk groups supervised by CVD ================

# study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_carrib.csv"))

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
                              save_path = paste0(wd, res_dir, "conf_sup_carrib_full"))
res_conf_sup$runtime  # 70 mins, 4 classes
res_conf_sup_varadj <- baysc::swolca_var_adjust(res = res_conf_sup, save_res = TRUE, 
                                                save_path = paste0(wd, res_dir, "conf_sup_carrib_full_varadj"))
res_conf_sup <- res_conf_sup_varadj

# # Try with additional confounders
# # Additional confounders
# V_data <- as.data.frame(study_data_dropna %>% 
#   select(bmi_c3, hyp_lab_med_sr, diab_lab_med_sr, 
#          dyslip_lab_med_sr) %>%
#   mutate_all(as.factor))
# glm_form <- "~ bmi_c3 + hyp_lab_med_sr + diab_lab_med_sr + dyslip_lab_med_sr"
# set.seed(1)
# res_conf_sup_cov <- baysc::swolca(x_mat = x_mat, y_all = y_all, V_data = V_data,
#                                   glm_form = glm_form, 
#                               sampling_wt = sampling_wt, cluster_id = cluster_id, 
#                               stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
#                               adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
#                               burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
#                               save_path = paste0(wd, res_dir, "conf_sup_carrib_cov_full"))
# res_conf_sup_cov_varadj <- baysc::swolca_var_adjust(res = res_conf_sup_cov, save_res = FALSE)


### Check results
item_labels <- c("Age", "Gender", "Smoking", "Alcohol", "PA", "FH_Diab", 
                 "FH_CHD", "FH_Stroke", "Rel_Angina", "Rel_MI", "Rel_CABG")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes (USE NEW REORDER FUNCTION!)
# res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(4, 3, 2, 1))
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


# Framingham score 
study_temp <- study_data_dropna %>% 
  left_join(raw_derv_vars %>% select(ID, FRAME_CVD_RISK_10YR) %>%
              mutate(ID = as.integer(ID)), by = join_by(ID))
summary(study_temp$FRAME_CVD_RISK_10YR[study_temp$CVD_FRAME == 0])
summary(study_temp$FRAME_CVD_RISK_10YR[study_temp$CVD_FRAME == 1])
study_temp$conf_class <- res_conf_sup$estimates$c_all
# Framingham score, for each class
round(sapply(1:res_conf_sup$estimates$K_red, function(x) 
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

# study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_carrib.csv"))

### Supervised by CVD 
# Cardiometabolic indicators risk factors
x_mat <- as.matrix(study_data_dropna %>% 
                     select(bmi_c3, hyp_lab_med_sr, diab_lab_med_sr, 
                            dyslip_lab_med_sr))
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
                             save_path = paste0(wd, res_dir, "cmd_sup_carrib_full"))
res_cmd_sup$runtime  # 40 mins, 4 classes
res_cmd_sup_varadj <- baysc::swolca_var_adjust(res = res_cmd_sup, save_res = TRUE, 
                                                save_path = paste0(wd, res_dir, "cmd_sup_carrib_full_varadj"))
res_cmd_sup <- res_cmd_sup_varadj


### Try clustering all CVD risk factors at the same time
# Supervised by CVD 
# All risk factors
x_mat <- as.matrix(study_data_dropna %>% 
                     select(bmi_c3, hyp_lab_med_sr, diab_lab_med_sr, 
                            dyslip_lab_med_sr, AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL,
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
res_both_sup <- baysc::swolca(x_mat = x_mat, y_all = y_all, glm_form = glm_form, 
                             sampling_wt = sampling_wt, cluster_id = cluster_id, 
                             stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                             adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                             burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                             save_path = paste0(wd, res_dir, "both_sup_carrib_full"))
res_both_sup$runtime  # 1.5 hours, 4 classes

item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia", 
                 "Age", "Gender", "Smoking", "Alcohol", "PA", "FH_Diab", 
                 "FH_CHD", "FH_Stroke", "Rel_Angina", "Rel_MI", "Rel_CABG")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes (USE NEW REORDER FUNCTION!)
res_both_sup <- reorder_classes(res = res_both_sup, new_order = c(3, 1, 2))
# Use some of the changed functions
plot_pattern_probs(res_both_sup, item_labels = item_labels, 
                   categ_title = categ_title, y_title = y_title, num_rows = 2)
baysc::plot_pattern_profiles(res_both_sup, item_labels = item_labels,
                             categ_title = categ_title, y_title = y_title)
plot_class_dist(res_both_sup)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_both_sup$post_MCMC_out$dendrogram)
regr_coefs <- get_regr_coefs(res = res_both_sup)
plot_regr_coefs(regr_coefs = regr_coefs, res = res_both_sup)
plot_outcome_probs(res_both_sup)
# Average posterior probability: 0.95
mean(apply(res_both_sup$estimates$pred_class_probs, 1, max))


# # Try with additional confounders
# # Additional confounders
# V_data <- as.data.frame(study_data_dropna %>% 
#                           # Create composite family history score
#                           mutate(FH_CVD = case_when(
#                             DIAB_FAMHIST == 3 | FH_CHD == 3 | FH_STROKE == 3 | 
#                               EVER_ANGINA_RELATIVE == 3 | EVER_MI_RELATIVE == 3 | 
#                               EVER_CABG_RELATIVE == 3 ~ 1,
#                             .default = 0
#                           )) %>% 
#                           select(AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, 
#                                  GPAQ_LEVEL, FH_CVD) %>%
#                           mutate_all(as.factor))
# glm_form <- "~ AGE_C3 + GENDER + CIGARETTE_USE + ALCOHOL_USE + GPAQ_LEVEL + FH_CVD"
# set.seed(1)
# res_cmd_sup_cov <- baysc::swolca(x_mat = x_mat, y_all = y_all, V_data = V_data,
#                                   glm_form = glm_form, 
#                                   sampling_wt = sampling_wt, cluster_id = cluster_id, 
#                                   stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
#                                   adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
#                                   burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
#                                   save_path = paste0(wd, res_dir, "cmd_sup_carrib_cov_full"))
# res_cmd_sup_cov_varadj <- baysc::swolca_var_adjust(res = res_cmd_sup_cov, save_res = FALSE)
# 
# # Try with additional confounders VERSION 2
# # Add age_centered
# study_data_temp <- study_data_dropna %>%
#   left_join(raw_derv_vars %>% select(ID, AGE) %>% mutate(ID = as.numeric(ID)), 
#             by = join_by(ID)) %>%
#   mutate(age_cent = AGE - mean(AGE, na.rm = TRUE))
# # Additional confounders
# V_data <- as.data.frame(study_data_temp %>% 
#                           # Create composite family history score
#                           mutate(fh_bin = case_when(
#                             DIAB_FAMHIST == 3 | FH_CHD == 3 | FH_STROKE == 3 | 
#                               EVER_ANGINA_RELATIVE == 3 | EVER_MI_RELATIVE == 3 | 
#                               EVER_CABG_RELATIVE == 3 ~ 1,
#                             .default = 0
#                           ),
#                           gender_bin = ifelse(GENDER == 3, 1, 0), # 0=M, 1=F
#                           cig_bin = ifelse(CIGARETTE_USE == 3, 1, 0), # 1 = current
#                           alc_bin = ifelse(ALCOHOL_USE == 3, 1, 0), # 1 = current
#                           pa_bin = ifelse(GPAQ_LEVEL == 3, 1, 0), # 1 = sedentary
#                           ) %>%  
#                           select(age_cent, gender_bin, cig_bin, alc_bin, pa_bin, 
#                                  fh_bin) %>%
#                           mutate_at(c("gender_bin", "cig_bin", "alc_bin", 
#                                       "pa_bin", "fh_bin"), as.factor))
# glm_form <- "~ age_cent + gender_bin + cig_bin + alc_bin + pa_bin + fh_bin"
# set.seed(1)
# res_cmd_sup_cov_v2 <- baysc::swolca(x_mat = x_mat, y_all = y_all, V_data = V_data,
#                                  glm_form = glm_form, 
#                                  sampling_wt = sampling_wt, cluster_id = cluster_id, 
#                                  stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
#                                  adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
#                                  burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
#                                  save_path = paste0(wd, res_dir, "cmd_sup_carrib_cov_full_v2"))
# res_cmd_sup_cov_varadj_v2 <- baysc::swolca_var_adjust(res = res_cmd_sup_cov_v2, 
#                                                       save_res = FALSE)
# # Posterior prob: 75, K=3, order: 3,2,1
# 
# # Try with additional confounders VERSION 3
# # Additional confounders
# V_data <- as.data.frame(study_data_dropna %>% 
#                           # Create composite family history score
#                           mutate(fh_bin = case_when(
#                             DIAB_FAMHIST == 3 | FH_CHD == 3 | FH_STROKE == 3 | 
#                               EVER_ANGINA_RELATIVE == 3 | EVER_MI_RELATIVE == 3 | 
#                               EVER_CABG_RELATIVE == 3 ~ 1,
#                             .default = 0
#                           ),
#                           gender_bin = ifelse(GENDER == 3, 1, 0), # 0=M, 1=F
#                           cig_bin = ifelse(CIGARETTE_USE == 3, 1, 0), # 1 = current
#                           alc_bin = ifelse(ALCOHOL_USE == 3, 1, 0), # 1 = current
#                           pa_bin = ifelse(GPAQ_LEVEL == 3, 1, 0), # 1 = sedentary
#                           ) %>%  
#                           select(AGE_C3, gender_bin, cig_bin, alc_bin, pa_bin, 
#                                  fh_bin) %>%
#                           mutate_at(c("AGE_C3", "gender_bin", "cig_bin", "alc_bin", 
#                                       "pa_bin", "fh_bin"), as.factor))
# glm_form <- "~ AGE_C3 + gender_bin + cig_bin + alc_bin + pa_bin + fh_bin"
# set.seed(1)
# res_cmd_sup_cov_v3 <- baysc::swolca(x_mat = x_mat, y_all = y_all, V_data = V_data,
#                                     glm_form = glm_form, 
#                                     sampling_wt = sampling_wt, cluster_id = cluster_id, 
#                                     stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
#                                     adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
#                                     burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
#                                     save_path = paste0(wd, res_dir, "cmd_sup_carrib_cov_full_v3"))
# res_cmd_sup_cov_varadj_v3 <- baysc::swolca_var_adjust(res = res_cmd_sup_cov_v3, 
#                                                       save_res = FALSE)
# # Posterior prob: 80, 4 classes
# 
# 
# # Try unsupervised
# set.seed(1)
# res_cmd <- baysc::wolca(x_mat = x_mat, sampling_wt = sampling_wt, cluster_id = cluster_id, 
#                                  stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
#                                  adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
#                                  burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
#                                  save_path = paste0(wd, res_dir, "cmd_carrib_full"))
# # Instability in variance adjustment. H_inv distance 14
# res_cmd_varadj <- baysc::wolca_var_adjust(res = res_cmd, save_res = FALSE)
# res_cmd <- reorder_classes(res = res_cmd, new_order = c(3, 4, 2, 5, 1))
# # Posterior probability: 62



### Check results
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
# res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(4, 1, 3, 2))
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
# Average posterior probability: 69
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
round(sapply(1:res_cmd_sup$estimates$K_red, function(x) 
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
# load(paste0(wd, res_dir, "conf_sup_carrib_full_varadj_swolca_results.RData"))
# res_conf_sup <- res
# res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(4, 3, 2, 1))
# load(paste0(wd, res_dir, "cmd_sup_carrib_full_varadj_swolca_results.RData"))
# res_cmd_sup <- res
# res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(4, 1, 3, 2))
# study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_carrib.csv"))

# Look at grid of confounder and cardiometabolic risk factors
table(res_conf_sup$estimates$c_all, res_cmd_sup$estimates$c_all, useNA = "always")
round(prop.table(table(res_conf_sup$estimates_adjust$c_all)), 3)
round(prop.table(table(res_cmd_sup$estimates_adjust$c_all)), 3)
# By levels of the confounder risk group
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 1])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 2])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 3])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[res_conf_sup$estimates$c_all == 4])), 2)



### Restrict to those without CVD AND those without self-report
condition <- (study_data_dropna$CVD_FRAME == 0) & 
  (study_data_dropna$hyp_sr == 0) & (study_data_dropna$diab_sr == 0) & 
  (study_data_dropna$dyslip_sr == 0)
# Number of individuals remaining: 1670 / 5086
length(res_conf_sup$estimates$c_all[condition])
round(prop.table(table(res_conf_sup$estimates_adjust$c_all[condition])), 3)
round(prop.table(table(res_cmd_sup$estimates_adjust$c_all[condition])), 3)

# Look at grid of confounder and cardiometabolic risk factors
table(res_conf_sup$estimates$c_all[condition], 
      res_cmd_sup$estimates$c_all[condition], useNA = "always")
# By levels of the confounder risk group
round(prop.table(
  table(res_cmd_sup$estimates$c_all[condition & 
                                      res_conf_sup$estimates$c_all == 1])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[condition & 
                                      res_conf_sup$estimates$c_all == 2])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[condition & 
                                      res_conf_sup$estimates$c_all == 3])), 2)
round(prop.table(
  table(res_cmd_sup$estimates$c_all[condition & 
                                      res_conf_sup$estimates$c_all == 4])), 2)
# By levels of the cardiometabolic risk group
round(prop.table(
  table(res_conf_sup$estimates$c_all[condition & 
                                      res_cmd_sup$estimates$c_all == 1])), 2)
round(prop.table(
  table(res_conf_sup$estimates$c_all[condition & 
                                      res_cmd_sup$estimates$c_all == 2])), 2)
round(prop.table(
  table(res_conf_sup$estimates$c_all[condition & 
                                      res_cmd_sup$estimates$c_all == 3])), 2)
round(prop.table(
  table(res_conf_sup$estimates$c_all[condition & 
                                      res_cmd_sup$estimates$c_all == 4])), 2)


# Create subgroup
study_data_dropna <- study_data_dropna %>%
  mutate(conf_risk = res_conf_sup$estimates_adjust$c_all, 
         cmd_risk = res_cmd_sup$estimates_adjust$c_all) %>%
  mutate(subgroup = case_when(
    cmd_risk == 1 & conf_risk == 1 ~ 1,  # low : low
    cmd_risk == 1 & conf_risk == 2 ~ 2,  # low : low-med
    cmd_risk == 1 & conf_risk == 3 ~ 3,  # low : med
    cmd_risk == 1 & conf_risk == 4 ~ 4,  # low : high
    cmd_risk == 2 & conf_risk == 1 ~ 5,  # med-low : low
    cmd_risk == 2 & conf_risk == 2 ~ 6,  # med-low : low-med
    cmd_risk == 2 & conf_risk == 3 ~ 7,  # med-low : med
    cmd_risk == 2 & conf_risk == 4 ~ 8,  # med-low : high
    cmd_risk == 3 & conf_risk == 1 ~ 9,  # med : low
    cmd_risk == 3 & conf_risk == 2 ~ 10,  # med : low-med
    cmd_risk == 3 & conf_risk == 3 ~ 11,  # med : med
    cmd_risk == 3 & conf_risk == 4 ~ 12,  # med : high
    cmd_risk == 4 & conf_risk == 1 ~ 13,  # high : low
    cmd_risk == 4 & conf_risk == 2 ~ 14,  # high : low-med
    cmd_risk == 4 & conf_risk == 3 ~ 15,  # high : med
    cmd_risk == 4 & conf_risk == 4 ~ 16,  # high : high
    .default = NA
  ), 
  # Subgroup version 2, combining the two highest risk groups
  subgroup2 = case_when(
    cmd_risk == 1 & conf_risk == 1 ~ 1,  # low-low
    cmd_risk == 1 & conf_risk == 2 ~ 2,  # low-med
    cmd_risk == 1 & conf_risk %in% c(3, 4) ~ 3,  # low-high
    cmd_risk == 2 & conf_risk == 1 ~ 4,  # med-low
    cmd_risk == 2 & conf_risk == 2 ~ 5,  # med-med
    cmd_risk == 2 & conf_risk %in% c(3, 4) ~ 6,  # med-high
    cmd_risk %in% c(3, 4) & conf_risk == 1 ~ 7,  # high-low
    cmd_risk %in% c(3, 4) & conf_risk == 2 ~ 8,  # high-med
    cmd_risk %in% c(3, 4) & conf_risk %in% c(3, 4) ~ 9,  # high-high
    .default = NA
  ))
table(study_data_dropna$subgroup)
table(study_data_dropna$subgroup2)

# Save full sample data
study_data_dropna_full <- study_data_dropna

### Restrict to those without CVD AND those without self-report: n=1670 / 5086
study_data_dropna_subset <- study_data_dropna %>%
  filter(CVD_FRAME == 0 & hyp_sr == 0 & diab_sr == 0 & dyslip_sr == 0)
table(study_data_dropna_subset$subgroup)
table(study_data_dropna_subset$subgroup2)

# write.csv(study_data_dropna_full,
#           file = paste0(wd, data_dir, "study_data_dropna_full.csv"), 
#           row.names = FALSE)
# write.csv(study_data_dropna_subset,
#           file = paste0(wd, data_dir, "study_data_dropna_subset.csv"),
#           row.names = FALSE)


# Create plot of distribution of confounder groups by levels of cmd group
plot_conf_cmd_dist(data = study_data_dropna_subset, normalize = TRUE)
plot_conf_cmd_dist(data = study_data_dropna_subset, normalize = FALSE)

# Create plot of distribution of confounder groups by levels of cmd group
plot_conf_cmd_dist(data = study_data_dropna_full, normalize = TRUE)
plot_conf_cmd_dist(data = study_data_dropna_full, normalize = FALSE)


#================= WRPC model with CVD risk subgroups: subsetting ==============
# study_data_dropna_subset <- read.csv(paste0(wd, data_dir, "study_data_dropna_subset.csv"))

# Restrict to those without CVD AND those without self-report: n=1670 / 5086
study_data_dropna <- study_data_dropna_subset

### Run WRPC with subgroup version 1

# Categorical exposure matrix, nxJ, 1670x49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
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
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 20000, 
                 burn = 10000, thin = 5, update = 1000, switch = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_subset"))
res_wrpc$runtime # 3 hours, K=4


### Run UNWEIGHTED WRPC with subgroup version 1

# Categorical exposure matrix, nxJ, 1670x49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
# Stratum indicators, nx1
stratum_id <- NULL    
# Cluster indicators, nx1
cluster_id <- NULL                  
# Survey sampling weights, nx1
sampling_wt <- NULL  
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup))
seed <- 1
res_wrpc_unwt <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 20000, 
                 burn = 10000, thin = 5, update = 1000, switch = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_subset_unwt"))
res_wrpc_unwt$runtime  # 2.7 hours, K = 5


res_wrpc$runtime # 3 hours, K=4
res_wrpc$estimates$pi_med  # 4 classes
res_wrpc$post_MCMC_out$K_med
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_wrpc$post_MCMC_out$dendrogram_global)
mean(apply(res_wrpc$estimates$pred_global_class_probs, 1, max))
plot(res_wrpc$MCMC_out$pi_MCMC[, 1], type = "l")
plot(res_wrpc$post_MCMC_out$pi[, 1], type = "l")
boxplot(res_wrpc$post_MCMC_out$pi[, 1], ylim=c(0,1))
plot(res_wrpc$MCMC_out$theta_local_MCMC[, 1, 1, 4], type = "l")
plot(res_wrpc$MCMC_out$nu_MCMC[, 1, 1], type = "l")
hist(round(res_wrpc$estimates$nu_med, 2), breaks=30)
# Which foods are local
which(res_wrpc$estimates$nu_med < 0.5, arr.ind = TRUE)
plot(res_wrpc$MCMC_out$nu_MCMC[, 14, 10], type = "l")
plot(res_wrpc$MCMC_out$nu_MCMC[, 10, 16], type = "l")
res_wrpc$estimates$theta_local_modes[13, ]
res_wrpc$estimates$theta_global_modes[13, ]
res_wrpc$estimates$theta_local_modes 

## plot trace of pi
K <- res_wrpc$estimates$K_red
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for pi"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  plot(res_wrpc$estimates$pi_red[, k], type = "l", ylim = c(0,1), xlab = "", ylab = "")
}


## Plot results
# Source plotting functions
source(paste0(wd, code_dir, "wrpc_plotting_fns.R"))
# Define labels
item_labels <- colnames(fpq_49categ[, -c(1:3)])
categ_labels <- c("None", "Monthly", "Weekly", "Daily", "Daily+")
plot_wrpc_global_pattern_profiles(res = res_wrpc, 
                                  item_labels = item_labels, item_title = "Item",
                                  categ_labels = categ_labels) + 
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
subgroup_labels <- c("L:L", "L:LM", "L:M", "L:H", 
                        "LM:L", "LM:LM", "LM:M", "LM:H",
                        "M:L", "M:LM", "M:M", "M:H",
                        "H:L", "H:LM", "H:M", "H:H")
# subgroup_labels <- c("L:L", "L:M", "L:H", 
#                      "M:L", "M:M", "M:H",
#                      "H:L", "H:M", "H:H")
# subgroup_labels <- c("L", "LM", "M", "H")
plot_wrpc_local_pattern_profiles(res = res_wrpc, 
                                 item_labels = item_labels, item_title = "Item",
                                 categ_labels = categ_labels,
                                 subgroup_labels = subgroup_labels,
                                 subgroup_title = "CVD Risk Group")
plot_wrpc_allocation(res = res_wrpc, item_labels = item_labels, 
                     item_title = "Item",
                     subgroup_labels = subgroup_labels, 
                     subgroup_title = "CVD Risk Group")
plot_wrpc_local_profiles_allocation(res = res_wrpc, 
                                 item_labels = item_labels, item_title = "Item",
                                 categ_labels = categ_labels,
                                 subgroup_labels = subgroup_labels,
                                 subgroup_title = "CVD Risk Group")
# Proportion of each subpopulation following each pattern
table(res_wrpc$estimates$c_all, res_wrpc$data_vars$h_all)
plot_class_subgroup_dist(res = res_wrpc, subgroup_labels = subgroup_labels, 
                         subgroup_title = "Subgroup", normalize = TRUE)


### Run WRPC with subgroup version 2

# Categorical exposure matrix, nxJ, 2879x129
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
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
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 20000, 
                 burn = 10000, thin = 5, update = 1000, switch = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_subset_v2"))
res_wrpc$runtime # 2.2 hours



### Run WRPC with subgroup version 3: ONLY cardiometabolic risk groups

# Categorical exposure matrix, nxJ, 1670x49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$cmd_risk))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 20000, 
                 burn = 10000, thin = 5, update = 1000, switch = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_subset_cmdonly"))
res_wrpc$runtime # 1.7 hours, K=5


### Run WRPC with subgroup version 4: combined risk groups
study_data_dropna <- study_data_dropna_full %>%
  mutate(both_risk = res_both_sup$estimates$c_all) %>%
  filter(CVD_FRAME == 0 & hyp_sr == 0 & diab_sr == 0 & dyslip_sr == 0)

# Categorical exposure matrix, nxJ, 1670x49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$both_risk))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 20000, 
                 burn = 10000, thin = 5, update = 1000, switch = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_subset_both"))
res_wrpc$runtime # 1.7 hours, K=5


### Run WOLCA
# Categorical exposure matrix, nxJ, 1670x49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   

seed <- 1
res_wolca <- baysc::wolca(x_mat = x_mat, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                 burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_subset_wolca"))
res_wolca$runtime # 1.1 hours, K=5
baysc::plot_pattern_profiles(res_wolca, item_labels = item_labels)
plot_pattern_probs(res_wolca, item_labels = item_labels, num_rows = 7)
plot_class_dist(res_wolca)



#================= WRPC model with CVD risk subgroups: FULL SAMPLE =============
# study_data_dropna_full <- read.csv(paste0(wd, data_dir, "study_data_dropna_full.csv"))

# Include full sample including those with CVD and self-report
study_data_dropna <- study_data_dropna_full

### Run WRPC with subgroup version 1

# Categorical exposure matrix, nxJ, 5086x49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
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
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 20000, 
                 burn = 10000, thin = 5, update = 1000, switch = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_full"))
res_wrpc$runtime # 3 hours, K=6



### Run WRPC with subgroup version 2

# Categorical exposure matrix, nxJ, 5086x49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
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
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                 burn = 5000, thin = 5, update = 1000, switch = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_full_v2"))
res_wrpc$runtime



# study_data_dropna_full <- read.csv(paste0(wd, data_dir, "study_data_dropna_full.csv"))

### Run WRPC with subgroup version 1, excluding only those with CVD
### n = 3242
study_data_dropna <- study_data_dropna_full %>%
  filter(CVD_FRAME == 0)

# Categorical exposure matrix, nxJ, 3242x49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
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
                 burn = 5000, thin = 5, update = 1000, switch = 100, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_full_noCVD"))
res_wrpc$runtime # 2.5 hours, K=5
