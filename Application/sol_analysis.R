#====================
# Analyzing HCHS data
#====================

rm(list = ls())

library(tidyverse)  # data wrangling
library(baysc)      # bayesian survey clustering
library(haven)      # read sas file
library(naniar)     # missingness

wd <- "~/Documents/Github/WRPC/"
code_dir <- "Code/"
data_dir <- "Application/HCHCS_Data/"
res_dir <- "Results/"
sol_data <- read.csv(paste0(wd, data_dir, "FPQ129_CVD07FEB2023.csv"))
# Select diet variables
fpq_data <- sol_data %>% 
  select(ID, STRAT, PSU_ID, CENTER, BKGRD1_C6, fp2:fp1)

### Check site x background numbers
# Site: B=Bronx, C=Chicago, M=Miami, S=San Diego
# Background: 0=Dominican, 1=Central or South American, 2=Cuban, 3=Mexican, 4=Puerto Rican
# 5=More than one/Other heritage 
# For Mexicans, Bronx: 163, Chicago: 2006, Miami: 29, San Diego: 3212
table(fpq_data$CENTER, fpq_data$BKGRD1_C6)
# Subset to Mexican background for Chicago and San Diego: total 5218
fpq_mex <- fpq_data %>%
  filter(CENTER %in% c("C", "S") & BKGRD1_C6 == 3)


# Large population centers: 
# Bronx-DR (1148)
# Bronx-PR (1449)
# [Chicago-C/SA (647)]
# [Chicago-PR (674)]
# Chicago-Mex (2006) 
# Miami-C/SA (1357)
# Miami-Cuban (2089)
# San Diego-Mex (3212)
# Total: 6 groups = 11261; 8 groups = 12582

raw_derv_vars <- read_sas(paste0(wd, data_dir, "part_derv_inv4.sas7bdat"))
raw_fpq_vars <- read_sas(paste0(wd, data_dir, "fpe_inv4.sas7bdat"))

table(raw_derv_vars$HYPERTENSION2, useNA = "always")
table(raw_derv_vars$HYPERT_AWARENESS, useNA = "always")
table(raw_derv_vars$HYPERTENSION2, raw_derv_vars$HYPERT_AWARENESS, 
      useNA = "always")

table(raw_derv_vars$DYSLIPIDEMIA, useNA = "always")

table(raw_derv_vars$DIABETES2, useNA = "always")

par(mfrow = c(1,1))
par(mar=c(5,4,4,2)+0.1)
hist(raw_derv_vars$YRSUS)
hist(raw_derv_vars$AGE_IMMI)
table(raw_derv_vars$US_BORN)
yrsus_3cat <- raw_derv_vars %>% mutate(
  YRSUS_C3 = case_when(
    US_BORN == 1 ~ 1, 
    US_BORN == 0 & YRSUS >= 10 ~ 2,
    US_BORN == 0 & YRSUS < 10 ~ 3,
    .default = NA
))
table(yrsus_3cat$YRSUS_C3, useNA = "always")


# Subset to those with FPQ data
raw_derv_sub <- raw_derv_vars %>% filter(ID %in% fpq_data$ID)
# Create variable combining site and ethnicity for the largest sample sizes
raw_derv_sub <- raw_derv_sub %>% mutate(
  site_eth = case_when(
    CENTER == "B" & BKGRD1_C6 == 0 ~ "Bronx_DR",
    CENTER == "B" & BKGRD1_C6 == 4 ~ "Bronx_PR",
    CENTER == "C" & BKGRD1_C6 == 1 ~ "Chicago_CA_SA",
    CENTER == "C" & BKGRD1_C6 == 4 ~ "Chicago_PR",
    CENTER == "C" & BKGRD1_C6 == 3 ~ "Chicago_Mex",
    CENTER == "M" & BKGRD1_C6 == 1 ~ "Miami_CA_SA",
    CENTER == "M" & BKGRD1_C6 == 2 ~ "Miami_Cuban",
    CENTER == "S" & BKGRD1_C6 == 3 ~ "SanDiego_Mex",
  )
)
# Distribution of various subgroups by site/ethnicity
table(raw_derv_sub$site_eth, raw_derv_sub$HYPERTENSION2)
table(raw_derv_sub$site_eth, raw_derv_sub$DIABETES2)
table(raw_derv_sub$site_eth, raw_derv_sub$DYSLIPIDEMIA)


# Look at distribution of various variables
vars_summ <- raw_derv_vars %>% 
  select(AGE, GENDER, 
         CIGARETTE_PACK_YEARS_C3, # 1: 0, 2: 0-10, 3: >10. 431 NAs
         CIGARETTE_USE, # 1: never, 2: former, 3: current. 93 NAs
         ALCOHOL_USE_LEVEL, # 1: no current, 2: low level, 3: high level. 140 NAs
         ALCOHOL_USE, # 1: never, 2: former, 3: current. 65 NAs
         GPAQ_LEVEL, # Physical activity. 1:high, 2: moderate, 3: low. 140 NAs
         DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE, # 0: No, 1: Yes. ~200 NAs
         EVER_MI_RELATIVE, EVER_CABG_RELATIVE,  # 0: No, 1: Yes. ~200 NAs
         CKD, # 1:normal, 2:mild, 3:moderate, 4:severe, 5:end-stage. 168 NAs
         AHI_GE15, # Sleep apnea. 1=yes. # 1946 NAs
         COPD_BY_BD,  # 1=yes. 15244 NAs
         MED_LLD, # lipid lowering drugs/antihyperlipidemics. 1=yes. 394 NAs
         MED_ANTIDIAB, # Antidiabetics. 1=yes. 394 NAs
         MED_ANTIHYPERT, # Antihypertensives. 1=yes. 394 NAs
  )
dim(vars_summ)
summary(vars_summ$AGE)
apply(vars_summ[, -1], 2, function(x) 
  round(prop.table(table(x, useNA = "always")), 2))
apply(vars_summ[, -1], 2, function(x) table(x, useNA = "always"))

# Repeat for ages 30 to 74
raw_derv_vars_agerestr <- raw_derv_vars %>%
  filter(AGE >= 30 & AGE <= 74)
vars_summ <- raw_derv_vars_agerestr %>% 
  select(AGE, GENDER, 
         CIGARETTE_PACK_YEARS_C3, # 1: 0, 2: 0-10, 3: >10. 360 NAs
         CIGARETTE_USE, # 1: never, 2: former, 3: current. 80 NAs
         ALCOHOL_USE_LEVEL, # 1: no current, 2: low level, 3: high level. 66 NAs
         ALCOHOL_USE, # 1: never, 2: former, 3: current. 61 NAs
         GPAQ_LEVEL, # Physical activity. 1:high, 2: moderate, 3: low. 106 NAs
         DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE, # 0: No, 1: Yes. ~150 NAs
         EVER_MI_RELATIVE, EVER_CABG_RELATIVE,  # 0: No, 1: Yes. ~150 NAs
         CKD, # 1:normal, 2:mild, 3:moderate, 4:severe, 5:end-stage. 125 NAs
         AHI_GE15, # Sleep apnea. 1=yes. # 1577 NAs
         COPD_BY_BD,  # 1=yes. 12736 NAs
         MED_LLD, # lipid lowering drugs/antihyperlipidemics. 1=yes. 338 NAs
         MED_ANTIDIAB, # Antidiabetics. 1=yes. 338 NAs
         MED_ANTIHYPERT, # Antihypertensives. 1=yes. 338 NAs
  )
dim(vars_summ)  # n=13730
summary(vars_summ$AGE)
apply(vars_summ[, -1], 2, function(x) 
  round(prop.table(table(x, useNA = "always")), 2))
apply(vars_summ[, -1], 2, function(x) table(x, useNA = "always"))


# Repeat for ages 30 to 74 and Puerto Rican ethnicity
raw_derv_vars_age_PR <- raw_derv_vars %>%
  filter(AGE >= 30 & AGE <= 74 & BKGRD1_C6 == 4)
vars_summ <- raw_derv_vars_age_PR %>% 
  select(AGE, GENDER, 
         CIGARETTE_PACK_YEARS_C3, # 1: 0, 2: 0-10, 3: >10. 43 NAs
         CIGARETTE_USE, # 1: never, 2: former, 3: current. 6 NAs
         ALCOHOL_USE_LEVEL, # 1: no current, 2: low level, 3: high level. 4 NAs
         ALCOHOL_USE, # 1: never, 2: former, 3: current. 3 NAs
         GPAQ_LEVEL, # Physical activity. 1:high, 2: moderate, 3: low. 12 NAs
         DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE, # 0: No, 1: Yes. ~40 NAs
         EVER_MI_RELATIVE, EVER_CABG_RELATIVE,  # 0: No, 1: Yes. ~30 NAs
         CKD, # 1:normal, 2:mild, 3:moderate, 4:severe, 5:end-stage. 40 NAs
         AHI_GE15, # Sleep apnea. 1=yes. # 356 NAs
         COPD_BY_BD,  # 1=yes. 2131 NAs
         MED_LLD, # lipid lowering drugs/antihyperlipidemics. 1=yes. 64 NAs
         MED_ANTIDIAB, # Antidiabetics. 1=yes. 64 NAs
         MED_ANTIHYPERT, # Antihypertensives. 1=yes. 64 NAs
  )
dim(vars_summ)  # n=2337
summary(vars_summ$AGE)
apply(vars_summ[, -1], 2, function(x) 
  round(prop.table(table(x, useNA = "always")), 2))
apply(vars_summ[, -1], 2, function(x) table(x, useNA = "always"))


# Repeat for Puerto Rican ethnicity
raw_derv_vars_PR <- raw_derv_vars %>%
  filter(BKGRD1_C6 == 4)
vars_summ <- raw_derv_vars_PR %>% 
  select(AGE, GENDER, 
         CIGARETTE_PACK_YEARS_C3, # 1: 0, 2: 0-10, 3: >10. 47 NAs
         CIGARETTE_USE, # 1: never, 2: former, 3: current. 7 NAs
         ALCOHOL_USE_LEVEL, # 1: no current, 2: low level, 3: high level. 4 NAs
         ALCOHOL_USE, # 1: never, 2: former, 3: current. 3 NAs
         GPAQ_LEVEL, # Physical activity. 1:high, 2: moderate, 3: low. 16 NAs
         DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE, # 0: No, 1: Yes. ~50 NAs
         EVER_MI_RELATIVE, EVER_CABG_RELATIVE,  # 0: No, 1: Yes. ~35 NAs
         CKD, # 1:normal, 2:mild, 3:moderate, 4:severe, 5:end-stage. 52 NAs
         AHI_GE15, # Sleep apnea. 1=yes. # 437 NAs
         COPD_BY_BD,  # 1=yes. 2489 NAs
         MED_LLD, # lipid lowering drugs/antihyperlipidemics. 1=yes. 70 NAs
         MED_ANTIDIAB, # Antidiabetics. 1=yes. 70 NAs
         MED_ANTIHYPERT, # Antihypertensives. 1=yes. 70 NAs
  )
dim(vars_summ)  # n=2728
summary(vars_summ$AGE)
apply(vars_summ[, -1], 2, function(x) 
  round(prop.table(table(x, useNA = "always")), 2))
apply(vars_summ[, -1], 2, function(x) table(x, useNA = "always"))

#================= FPQ cleaning ================================================
# Combined study data: 12519 x 394
study_data <- raw_derv_sub %>% 
  mutate(ID = as.integer(ID)) %>%
  left_join(fpq_data %>% 
              select(ID, fp2:fp1) %>% 
              mutate(across(fp2:fp1, as.numeric)), by = join_by(ID))

# Find 5-level foods
which(apply(fpq_data[, -c(1:5)], 2, 
            function(x) max(x, na.rm = TRUE) == 5))
# Convert 5-level foods to 4 levels
study_data <- study_data %>%
  mutate_at(c("fp2", "fp116", "fp7", "fp10", "fp11", "fp12", "fp13", "fp17", 
              "fp53", "fp101"), ~ifelse(.x == 5, 4, .x))
# Check no 5-level foods
apply(study_data %>% select(fp2:fp1), 2, function(x) max(x, na.rm = TRUE))


### Define subgroups
study_data <- study_data %>%
  mutate(subgroup = case_when(
    AGE < 45 & GENDER == "F" ~ 1, #"f_young",
    AGE >= 45 & AGE < 65 & GENDER == "F" ~ 2, #"f_mid",
    AGE >= 54 & GENDER == "F" ~ 3, #"f_old",
    AGE < 45 & GENDER == "M" ~ 4, #"m_young",
    AGE >= 45 & AGE < 65 & GENDER == "M" ~ 5, #"m_mid",
    AGE >= 54 & GENDER == "M" ~ 6, #"m_old",
    .default = NA
  ))


# Check missingness
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             subgroup, fp2:fp1))
# naniar::vis_miss(study_data %>% 
#                    select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
#                           subgroup, fp2:fp1))

# Drop missing: sample size is 11519 x 394
# Note: no missingness for ID, STRAT, PSU_ID, WEIGHT
# DIABETES2 has 3 missing
# fp2:fp1 has ~1000 missing
study_data <- study_data %>%
  drop_na(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, subgroup, fp2:fp1)

# Distribution of various subgroups
table(study_data$subgroup)
# Distribution of various subgroups by site/ethnicity
table(study_data$site_eth, study_data$subgroup)


#================= Process CVD risk factor variables ===========================
cvd_vars <- raw_derv_vars %>%
  mutate(bmi_c4 = BMIGRP_C4,   # 1=underweight, 2=normal, 3=overweight, 4=obese
         dyslipidemia_c3 = case_when(
           DYSLIPIDEMIA == 0 ~ 1,  # none
           DYSLIPIDEMIA == 1 & MED_LLD == 1 ~ 2,  # treated
           DYSLIPIDEMIA == 1 & MED_LLD == 0 ~ 3,  # untreated
           .default = NA
         ),
         diabetes_c4 = DIABETES_C4,  # 1=none, 2=pre, 3=treated, 4=untreated
         hypertension_c4 = HYPERTENSION_C4)  


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
