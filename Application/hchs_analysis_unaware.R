#====================================================
# Analyzing HCHS data considering condition awareness
# Author: Stephanie Wu
# Date Updated: 2024/11/13
#====================================================

#==================== Preparing HCHS data w/ PR, DR, and Cuban =================

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

#================= FPQ cleaning ================================================

# Select diet variables
fpq_data <- sol_data %>% 
  select(ID, STRAT, PSU_ID, CENTER, BKGRD1_C6, fp2:fp1)

### Check site x background numbers
# Site: B=Bronx, C=Chicago, M=Miami, S=San Diego
# Background: 0=Dominican, 1=Central or South American, 2=Cuban, 3=Mexican, 
# 4=Puerto Rican, 5=More than one/Other heritage 
table(fpq_data$CENTER, fpq_data$BKGRD1_C6)

# Drop NAs: n=12761
fpq_dropna <- fpq_data %>% drop_na()

# ### Subset to Puerto Rican and Dominican background and complete data: total 3218
# fpq_PR_DR <- fpq_data %>% 
#   filter(BKGRD1_C6 %in% c(0, 4)) %>%
#   drop_na()

# Merge in derived variables, which include demographic and outcome variables.
# Merge in self-report variables
# Combined study data: 12761 x 521
study_data <- raw_derv_vars %>% 
  mutate(ID = as.integer(ID)) %>%
  left_join(raw_mhea_vars %>% mutate(ID = as.integer(ID)), by = join_by(ID)) %>%
  right_join(fpq_dropna %>% 
               select(ID, fp2:fp1) %>% 
               mutate(across(fp2:fp1, as.numeric)), by = join_by(ID)) 
  
# Find 5-level foods
which(apply(fpq_dropna[, -c(1:5)], 2, 
            function(x) max(x, na.rm = TRUE) == 5))
# Convert 5-level foods to 4 levels
study_data <- study_data %>%
  mutate_at(c("fp2", "fp116", "fp7", "fp10", "fp11", "fp12", "fp13", "fp17", 
              "fp53", "fp101"), ~ifelse(.x == 5, 4, .x))
# Check no 5-level foods
apply(study_data %>% select(fp2:fp1), 2, function(x) max(x, na.rm = TRUE))


#==================== Confounders cleaning ========================================

# # Look at distribution of various variables
# vars_summ <- study_data %>% 
#   select(AGE, GENDER, WAIST_HIP,
#          CIGARETTE_PACK_YEARS_C3, # 1: 0, 2: 0-10, 3: >10. 37 NAs
#          CIGARETTE_USE, # 1: never, 2: former, 3: current. 5 NAs
#          ALCOHOL_USE_LEVEL, # 1: no current, 2: low level, 3: high level. 4 NAs
#          ALCOHOL_USE, # 1: never, 2: former, 3: current. 3 NAs
#          GPAQ_LEVEL, # Physical activity. 1:high, 2: moderate, 3: low. 11 NAs
#          DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE, # 0: No, 1: Yes. ~44 NAs
#          EVER_MI_RELATIVE, EVER_CABG_RELATIVE,  # 0: No, 1: Yes. ~24 NAs
#          CKD, # 1:normal, 2:mild, 3:moderate, 4:severe, 5:end-stage. 15 NAs
#          AHI_GE15, # Sleep apnea. 1=yes. # 412 NAs
#          COPD_BY_BD,  # 1=yes. 2944 NAs
#          MED_LLD, # lipid lowering drugs/antihyperlipidemics. 1=yes. 64 NAs
#          MED_ANTIDIAB, # Antidiabetics. 1=yes. 64 NAs
#          MED_ANTIHYPERT # Antihypertensives. 1=yes. 64 NAs
#   )
# dim(vars_summ)  # n=3218
# summary(vars_summ$AGE)
# summary(vars_summ$WAIST_HIP)
# apply(vars_summ[, -1], 2, function(x) 
#   round(prop.table(table(x, useNA = "always")), 2))
# apply(vars_summ[, -1], 2, function(x) table(x, useNA = "always"))

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
    bmi_c4 = case_when(
      BMIGRP_C6 %in% c(1, 2) ~ 1, # normal or underweight BMI<25
      BMIGRP_C6 == 3 ~ 2, # overweight 25<=BMI<30
      BMIGRP_C6 == 4 ~ 3, # obese 30<=BMI<35
      BMIGRP_C6 %in% c(5, 6) ~ 4, # obese 2 or 3 (BMI>=35)
      .default = NA),   
    hyp_c3 = case_when(
      HYPERTENSION == 0 ~ 1, # none
      HYPERTENSION == 1 & MHEA1 == 1 ~ 3, # hypertension aware (lab-based or scanned meds, & self-report)
      HYPERTENSION == 1 & MHEA1 == 0 ~ 4, # hypertension unaware
      .default = NA),
    diab_c4 = case_when(
      DIABETES2 == 1 ~ 1, # none 
      DIABETES2 == 2 ~ 2, # pre-diabetic (fasting time > 8 AND fasting glucose in range 100-125 mg/dL) or (post-OGTT glucose in range 140-199 mg/dL) or (5.7%≤A1C<6.5%)
      DIABETES2 == 3 & MHEA16 == 1 ~ 3, # diabetes aware (FBS>=>200 or post-OGTT>=200 or A1C>=6.5% or scanned meds, & self-report)
      DIABETES2 == 3 & MHEA16 == 0 ~ 4, # diabetes unaware
      .default = NA),
    dyslip_c3 = case_when(
      DYSLIPIDEMIA == 0 ~ 1, # none
      (DYSLIPIDEMIA == 1 | MED_LLD == 1) & (MHEA2 == 1) ~ 3, # high chol aware (LDL>=160 or HDL<=40 or trigl>=200 or scanned meds, & self-report)
      (DYSLIPIDEMIA == 1 | MED_LLD == 1) & (MHEA2 == 0) ~ 4, # high chol unaware 
      .default = NA),
    
    # Additional old variables
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
      CKD %in% c(3, 4, 5) ~ 3), # moderate, severe, end-stage, GFR<60
    hyp_c5 = case_when(
      HYPERTENSION_C4 == 1 ~ 1, # none 
      HYPERTENSION_C4 == 2 ~ 2, # pre-hypertension 120<=SBP<140
      HYPERTENSION_C4 == 3 ~ 3, # hypertension treated 
      HYPERTENSION_C4 == 4 & MHEA1 == 1 ~ 4, # hypertension untreated but aware
      HYPERTENSION_C4 == 4 & MHEA1 == 0 ~ 5, # hypertension untreated and unaware
      .default = NA),
    diab_c5 = case_when(
      DIABETES_C4 == 1 ~ 1, # none 
      DIABETES_C4 == 2 ~ 2, # pre-diabetic (fasting time > 8 AND fasting glucose in range 100-125 mg/dL) or (post-OGTT glucose in range 140-199 mg/dL) or (5.7%≤A1C<6.5%)
      DIABETES_C4 == 3 ~ 3, # diabetes treated (FBS>=>200 or post-OGTT>=200 or A1C>=6.5% or scanned meds)
      DIABETES_C4 == 4 & MHEA16 == 1 ~ 4, # diabetes untreated but aware
      DIABETES_C4 == 4 & MHEA16 == 0 ~ 5, # diabetes untreated and unaware
      .default = NA)
    )  

# Check missingness
naniar::gg_miss_var(study_data_conf_cmd %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL, 
                             DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                             EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
                             bmi_c4, hyp_c3, diab_c4, dyslip_c3))
naniar::vis_miss(study_data_conf_cmd %>% 
                   select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                          AGE_C3, GENDER, CIGARETTE_USE, ALCOHOL_USE, GPAQ_LEVEL, 
                          DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE,
                          EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
                          bmi_c4, hyp_c3, diab_c4, dyslip_c3))

# Sample size after dropping all missingness: n = 12480
study_data_dropna <- study_data_conf_cmd %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, CENTER, BKGRD1_C6,
         AGE_C3, GENDER, CIGARETTE_USE, 
         ALCOHOL_USE, GPAQ_LEVEL, DIAB_FAMHIST, FH_CHD, FH_STROKE, 
         EVER_ANGINA_RELATIVE, EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
         bmi_c4, hyp_c3, diab_c4, dyslip_c3, fp2:fp1) %>% 
  drop_na()

# CVD levels: 3609 w/ CVD (will exclude from RPC)
table(study_data_dropna$CVD_FRAME, useNA = "always")

# # Save cleaned data
# write.csv(study_data_dropna, file = paste0(wd, data_dir, "cleaned_data_unaware_full.csv"),
#           row.names = FALSE)


### Sample size restricted to exclude those with any self-reported conditions 
# n = 7816
study_data_unaware <- study_data_dropna %>%
  filter((hyp_c3 != 3) & (diab_c4 != 3) & (dyslip_c3 != 3))

table(study_data_dropna$hyp_c3)
table(study_data_dropna$diab_c4)
table(study_data_dropna$dyslip_c3)

table(study_data_unaware$hyp_c3)
table(study_data_unaware$diab_c4)
table(study_data_unaware$dyslip_c3)

### Check site x background numbers
# Site: B=Bronx, C=Chicago, M=Miami, S=San Diego
# Background: 0=Dominican, 1=Central or South American, 2=Cuban, 3=Mexican, 
# 4=Puerto Rican, 5=More than one/Other heritage 
table(study_data_unaware$CENTER, study_data_unaware$BKGRD1_C6)
table(study_data_unaware$BKGRD1_C6)

# CVD levels: 1578 w/ CVD (will exclude from RPC)
table(study_data_unaware$CVD_FRAME, useNA = "always")

# # Save cleaned data
# write.csv(study_data_unaware, file = paste0(wd, data_dir, "cleaned_data_unaware.csv"),
#           row.names = FALSE)


### Subset to Carribean backgrounds (Dominican, Cuban, and Puerto Rican)
# n = 2879
study_data_unaware_carrib <- study_data_unaware %>%
  filter(BKGRD1_C6 %in% c(0, 2, 4))

# CVD levels: 710 w/ CVD (will exclude from RPC)
table(study_data_unaware_carrib$CVD_FRAME, useNA = "always")

# # Save cleaned data
# write.csv(study_data_unaware_carrib, 
#           file = paste0(wd, data_dir, "cleaned_data_unaware_carrib.csv"),
#           row.names = FALSE)

study_data_dropna <- study_data_unaware_carrib

#============== Derive confounder risk groups supervised by CVD ================
# study_data_dropna <- read.csv(paste0(wd, data_dir, "cleaned_data_unaware_carrib.csv"))

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
seed <- 2
res_conf_sup <- baysc::swolca(x_mat = x_mat, y_all = y_all, glm_form = glm_form, 
                              sampling_wt = sampling_wt, cluster_id = cluster_id, 
                              stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                              adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                              burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                              save_path = paste0(wd, res_dir, "conf_sup_unaware_carrib_seed2"))
res_conf_sup$runtime  # 40 mins, 3 classes
item_labels <- c("Age", "Gender", "Smoking", "Alcohol", "PA", "FH_Diab", 
                 "FH_CHD", "FH_Stroke", "Rel_Angina", "Rel_MI", "Rel_CABG")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes (USE NEW REORDER FUNCTION!)
res_conf_sup <- reorder_classes(res = res_conf_sup, new_order = c(4, 3, 2, 1))
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
# Average posterior probability: 0.96
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
seed <- 2
res_cmd_sup <- baysc::swolca(x_mat = x_mat, y_all = y_all, glm_form = glm_form, 
                             sampling_wt = sampling_wt, cluster_id = cluster_id, 
                             stratum_id = stratum_id, run_sampler = "both", K_max = 30, 
                             adapt_seed = seed, class_cutoff = 0.05, n_runs = 20000, 
                             burn = 10000, thin = 5, update = 1000, save_res = TRUE, 
                             save_path = paste0(wd, res_dir, "cmd_sup_unaware_carrib_seed2"))
res_cmd_sup$runtime  # 23 mins, 4 classes
item_labels <- c("BMI", "Hypertension", "Diabetes", "Dyslipidemia")
categ_title <- "Risk Level"
y_title <- "Risk Level Probability"
# Reorder classes
res_cmd_sup <- reorder_classes(res = res_cmd_sup, new_order = c(2, 3, 1))
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
# Average posterior probability: 73
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





