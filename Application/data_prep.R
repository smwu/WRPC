#====================
# Preparing HCHS data
#====================

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

#================= FPQ cleaning ================================================

# Select diet variables
fpq_data <- sol_data %>% 
  select(ID, STRAT, PSU_ID, CENTER, BKGRD1_C6, fp2:fp1)

### Check site x background numbers
# Site: B=Bronx, C=Chicago, M=Miami, S=San Diego
# Background: 0=Dominican, 1=Central or South American, 2=Cuban, 3=Mexican, 
# 4=Puerto Rican, 5=More than one/Other heritage 
table(fpq_data$CENTER, fpq_data$BKGRD1_C6)


### Subset to Puerto Rican background and complete data: total 2066
fpq_PR <- fpq_data %>% 
  filter(BKGRD1_C6 == 4) %>%
  drop_na()

# Merge in derived variables, which include demographic and outcome variables
# Combined study data: 2219 x 393
study_data <- raw_derv_vars %>% 
  mutate(ID = as.integer(ID)) %>%
  right_join(fpq_PR %>% 
              select(ID, fp2:fp1) %>% 
              mutate(across(fp2:fp1, as.numeric)), by = join_by(ID))

# Find 5-level foods
which(apply(fpq_PR[, -c(1:5)], 2, 
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
         CIGARETTE_PACK_YEARS_C3, # 1: 0, 2: 0-10, 3: >10. 31 NAs
         CIGARETTE_USE, # 1: never, 2: former, 3: current. 4 NAs
         ALCOHOL_USE_LEVEL, # 1: no current, 2: low level, 3: high level. 4 NAs
         ALCOHOL_USE, # 1: never, 2: former, 3: current. 3 NAs
         GPAQ_LEVEL, # Physical activity. 1:high, 2: moderate, 3: low. 9 NAs
         DIAB_FAMHIST, FH_CHD, FH_STROKE, EVER_ANGINA_RELATIVE, # 0: No, 1: Yes. ~35 NAs
         EVER_MI_RELATIVE, EVER_CABG_RELATIVE,  # 0: No, 1: Yes. ~20 NAs
         CKD, # 1:normal, 2:mild, 3:moderate, 4:severe, 5:end-stage. 12 NAs
         AHI_GE15, # Sleep apnea. 1=yes. # 307 NAs
         COPD_BY_BD,  # 1=yes. 2026 NAs
         MED_LLD, # lipid lowering drugs/antihyperlipidemics. 1=yes. 48 NAs
         MED_ANTIDIAB, # Antidiabetics. 1=yes. 48 NAs
         MED_ANTIHYPERT # Antihypertensives. 1=yes. 48 NAs
  )
dim(vars_summ)  # n=2066
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

# # Read in self-report variables
# self_report <- read_sas("~/Downloads/mhea_inv4.sas7bdat")
# # Combine into derived variables
# all_vars <- study_data %>%
#   left_join(self_report %>% mutate(ID = as.numeric(ID)), by = join_by(ID))
# # Subset to PR
# all_vars_PR <- all_vars %>%
#   filter(ID %in% study_data$ID)
# # Hypertension
# # self-report
# round(prop.table(table(all_vars_PR$MHEA1, useNA = "always")), 3)  
# # lab or meds
# round(prop.table(table(all_vars_PR$HYPERTENSION, useNA = "always")), 3)  
# round(prop.table(table(all_vars_PR$HYPERT_AWARENESS, useNA = "always")), 3)  
# round(prop.table(table(all_vars_PR$DM_AWARE, useNA = "always")), 3)  


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

# Sample size after dropping all missingness: n = 1988
study_data_dropna <- study_data_conf_cmd %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, AGE_C3, GENDER, CIGARETTE_USE, 
         ALCOHOL_USE, GPAQ_LEVEL, DIAB_FAMHIST, FH_CHD, FH_STROKE, 
         EVER_ANGINA_RELATIVE, EVER_MI_RELATIVE, EVER_CABG_RELATIVE, CVD_FRAME,
         bmi_c3, hyp_lab_med, diab_lab_med, dyslip_lab_med, ckd_lab, fp2:fp1) %>% 
  drop_na()

# Save cleaned data
write.csv(study_data_dropna, file = paste0(wd, data_dir, "cleaned_data.csv"),
          row.names = FALSE)


#============ Exploratory data analysis ========================================
# Ignore strata and clustering, or else get error "stratum 23 has only one PSU 
# at stage 1
svydes <- survey::svydesign(ids=~1, weights=~WEIGHT_FINAL_EXPANDED, 
                            data = study_data_dropna)
# Get hypertension overall and by gender
# Higher prevalence among women
survey::svyby(~hyp_lab_med, ~GENDER, svydes, survey::svymean)

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



#============== Updated plotting functions =====================================

get_regr_coefs <- function (res, ci_level = 0.95, digits = 2) {
  if (!(inherits(res, c("swolca", "wolca")))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting \n         from a call to one of these functions")
  }
  if (!(ci_level > 0 & ci_level < 1)) {
    stop("ci_level must be between 0 and 1")
  }
  quant_lb <- (1 - ci_level)/2
  quant_ub <- 1 - quant_lb
  if (!is.null(res$estimates_adjust)) {
    estimates <- res$estimates_adjust
  } else {
    estimates <- res$estimates
  }
  terms <- labels(stats::terms(stats::as.formula(res$data_vars$glm_form)))
  if (length(terms) > 0) {
    full_glm_form <- paste0("y_all ~ ", paste0("c_all * ", 
                                               terms, collapse = " + "))
  } else {
    full_glm_form <- paste0("y_all ~ c_all")
  }
  full_data <- data.frame(c_all = as.factor(res$estimates$c_all), 
                          y_all = res$data_vars$y_all, res$data_vars$V_data)
  model_matrix <- model.matrix(as.formula(full_glm_form), 
                               data = full_data)
  beta <- as.data.frame(matrix(NA, nrow = ncol(model_matrix), 
                               ncol = 5))
  beta[, 1] <- colnames(model_matrix)
  if (is(res, "wolca")) {
    if (ci_level != res$data_vars$ci_level) {
      stop("ci_level must match the specified ci_level in the wolca() function")
    }
    colnames(beta) <- c("Covariate", "Estimate", "LB", "UB", 
                        "p-value")
    beta[, c(2, 5)] <- res$estimates_svyglm$fit_summary$coefficients[, 
                                                                     c(1, 4)]
    beta[, 2] <- format(round(beta[, 2], digits), digits)
    beta[, 3] <- format(round(convert_mix_to_ref(res$estimates_svyglm$xi_est_lb), 
                              digits), digits)
    beta[, 4] <- format(round(convert_mix_to_ref(res$estimates_svyglm$xi_est_ub), 
                              digits), digits)
    beta[, 5] <- ifelse(beta[, 5] < 10^(-digits), 
                        paste0("<", 10^(-digits)), 
                        format(round(beta[, 5], digits), digits))
  } else {
    est_xi <- estimates$xi_med
    est_lb <- apply(estimates$xi_red, c(2, 3), 
                    function(x) stats::quantile(x, quant_lb))
    est_ub <- apply(estimates$xi_red, c(2, 3), 
                    function(x) stats::quantile(x, quant_ub))
    est_red <- estimates$xi_red
    K <- nrow(est_xi)
    Q <- ncol(est_xi)
    colnames(beta) <- c("Covariate", "Estimate", "LB", "UB", 
                        "P(xi > 0)")
    beta[1, -1] <- c(est_xi[1, 1], 
                     get_ci(post_samp = est_red[, 1, 1]), 
                     get_prob_pos(est_red[, 1, 1]))
    ### CHANGED
    if (K > 1) {
      for (i in 2:K) {
        beta[i, -1] <- c(stats::median(est_red[, i, 1] - est_red[, 1, 1]), 
                         get_ci(est_red[, i, 1] - est_red[, 1, 1], digits = digits), 
                         get_prob_pos(est_red[, i, 1] - est_red[, 1, 1], digits = digits))
      }
    }
    if (Q > 1) {
      for (i in 2:Q) {
        beta[K + (i - 1), -1] <- c(est_xi[1, i], 
                                   get_ci(est_red[, 1, i]), 
                                   get_prob_pos(est_red[, 1, i]))
      }
      for (i in 2:Q) {
        for (j in 2:K) {
          beta[Q + (i - 1) * (K - 1) + (j - 1), -1] <- 
            c(stats::median(est_red[, j, i] - est_red[, 1, i]), 
              get_ci(est_red[, j, i] - est_red[, 1, i], digits = digits), 
              get_prob_pos(est_red[, j, i] - est_red[, 1, i], digits = digits))
        }
      }
    }
    ### END CHANGED
    beta$Estimate <- as.numeric(beta$Estimate)
    beta$LB <- as.numeric(beta$LB)
    beta$UB <- as.numeric(beta$UB)
  }
  beta <- dplyr::mutate_if(beta, is.numeric, round, digits = digits)
  beta
}


convert_to_probs <- function (est_xi, glm_form, V_data, cov_name = NULL) {
  ### CHANGED 
  if (!is.null(cov_name)) {
    if (!all(sapply(cov_name, function(x) grepl(x, glm_form)))) {
      stop("all variables in cov_name must be specified in glm_form")
    }
    else if (!all(sapply(cov_name, function(x) x %in% colnames(V_data)))) {
      stop("all variables in cov_name must be found in V_data")
    }
  }
  ### END CHANGED
  if (grepl("c_all", glm_form)) {
    stop("glm_form must not contain the latent class variable c_all")
  }
  K <- nrow(est_xi)
  ### CHANGED
  # Case with no covariates
  if (is.null(cov_name)) {
    Phi_df <- as.data.frame(t(stats::pnorm(est_xi)))
    # all_Phi_df <- as.data.frame(sapply(1:K, function(k) 
    #   stats::pnorm(est_xi[k, ])))
    colnames(Phi_df) <- paste0("Class", 1:K)
    # Phi_df <- Phi_df[1:num_cov_levels, ]
    
  } else { # Case with additional covariates
    cov_names <- labels(stats::terms(stats::as.formula(glm_form)))
    oth_names <- cov_names[cov_names != cov_name]
    if (!is.factor(V_data[[cov_name]])) { # factor
      stop("Please convert V_data cov_name data type to factor")
    } else {
    cov_levels <- lapply(cov_name, function(x) levels(V_data[[x]]))
    num_cov_levels <- nrow(expand.grid(cov_levels))
    oth_levels <- lapply(oth_names, function(x) levels(V_data[[x]]))
    all_levels <- append(cov_levels, oth_levels)
    all_level_comb <- expand.grid(all_levels)
    colnames(all_level_comb) <- c(cov_name, oth_names)
    all_model_mat <- stats::model.matrix(stats::as.formula(glm_form), 
                                         all_level_comb)
    all_Phi_df <- as.data.frame(sapply(1:K, function(k) stats::pnorm(all_model_mat %*% 
                                                                       est_xi[k, ])))
    colnames(all_Phi_df) <- paste0("Class", 1:K)
    colnames(all_level_comb)[1:length(cov_name)] <- paste0("Cov", 
                                                           1:length(cov_name))
    Phi_df <- cbind(all_Phi_df, all_level_comb)
    Phi_df <- Phi_df[1:num_cov_levels, ]
    }
  }
  ### END CHANGED
  
  
  return(Phi_df)
}


plot_outcome_probs <- function (res, cov_name = NULL, ci_level = 0.95, add_lines = FALSE, 
                                cov_labels = NULL, class_labels = NULL, class_title = "Dietary Pattern", 
                                x_title = NULL, y_title = "Probability of Outcome", ...) {
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting \n         from a call to one of these functions")
  }
  if (!is.null(ci_level)) {
    if (!(ci_level > 0 & ci_level < 1)) {
      stop("ci_level must be between 0 and 1")
    }
    quant_lb <- (1 - ci_level)/2
    quant_ub <- 1 - quant_lb
  }
  ### CHANGED
  if (!is.null(cov_name)) {
    if (!(length(cov_name) %in% c(1, 2)) | !(is.character(cov_name[[1]]))) {
      stop("cov_name must be `NULL` or a string vector of length 1 or 2")
    }
  }
  if (!is.null(cov_labels)) {
    if (!is.list(cov_labels)) {
      cov_labels <- list(cov_labels)
    }
    if (!(length(cov_labels) %in% c(1, 2)) | !(is.character(cov_labels[[1]]))) {
      print(cov_labels)
      stop("cov_labels must be a list of length 1 or 2 composed of string vectors")
    }
    else if (length(cov_name) != length(cov_labels)) {
      stop(paste0("cov_name is a vector of length ", length(cov_name), 
                  ", while cov_labels is a list of length ", length(cov_labels), 
                  ". The two must be of the same length."))
    }
  }
  ### END CHANGED
  if (is(res, "wolca")) {
    if (is.null(res$estimates_svyglm)) {
      stop("wolca object does not have regression estimates. Please run wolca_svyglm().")
    }
    estimates <- res$estimates_svyglm
    est_xi <- estimates$xi_est
    if (!is.null(ci_level)) {
      if (ci_level != res$data_vars$ci_level) {
        stop("ci_level must match the specified ci_level in the wolca() function")
      }
      est_lb <- estimates$xi_est_lb
      est_ub <- estimates$xi_est_ub
    }
  } else {
    if (!is.null(res$estimates_adjust)) {
      estimates <- res$estimates_adjust
    }
    else {
      estimates <- res$estimates
    }
    est_xi <- estimates$xi_med
    if (!is.null(ci_level)) {
      est_lb <- apply(estimates$xi_red, c(2, 3), function(x) stats::quantile(x, 
                                                                             quant_lb))
      est_ub <- apply(estimates$xi_red, c(2, 3), function(x) stats::quantile(x, 
                                                                             quant_ub))
    }
  }
  K <- nrow(est_xi)
  Phi_df <- convert_to_probs(est_xi = est_xi, glm_form = res$data_vars$glm_form, 
                             V_data = res$data_vars$V_data, cov_name = cov_name)
  if (!is.null(ci_level)) {
    Phi_lb <- convert_to_probs(est_xi = est_lb, glm_form = res$data_vars$glm_form, 
                               V_data = res$data_vars$V_data, cov_name = cov_name)
    Phi_ub <- convert_to_probs(est_xi = est_ub, glm_form = res$data_vars$glm_form, 
                               V_data = res$data_vars$V_data, cov_name = cov_name)
  }
  Phi_df_long <- Phi_df %>% tidyr::pivot_longer(cols = 1:K, 
                                                names_to = "Class", values_to = "Phi")
  if (!is.null(ci_level)) {
    Phi_lb_long <- Phi_lb %>% tidyr::pivot_longer(cols = 1:K, 
                                                  names_to = "Class", values_to = "Phi_lb")
    Phi_ub_long <- Phi_ub %>% tidyr::pivot_longer(cols = 1:K, 
                                                  names_to = "Class", values_to = "Phi_ub")
    col_names <- colnames(Phi_df_long)[-ncol(Phi_df_long)]
    Phi_df_long <- Phi_df_long %>% dplyr::left_join(Phi_lb_long, 
                                                    by = col_names) %>% dplyr::left_join(Phi_ub_long, 
                                                                                         by = col_names)
  }
  ### CHANGED
  if (is.null(x_title)) {
    if (!is.null(cov_name)) {
      x_title <- cov_name[1]
    } else {
      x_title <- "Classes"
    }
    
  }
  if (!is.null(cov_name)) {
    if (is.null(cov_labels)) {
      cov_labels <- lapply(cov_name, function(x) levels(res$data_vars$V_data[[x]]))
    }
    else {
      for (i in 1:length(cov_labels)) {
        num_categs <- length(levels(res$data_vars$V_data[[cov_name[i]]]))
        if (length(cov_labels[[i]]) != num_categs) {
          stop(paste0("length of cov_labels for covariate ", 
                      cov_name[i], " must equal the number of categories: ", 
                      num_categs))
        }
      }
    }
  }
  ### END CHANGED
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", 
                K))
  }
  Class <- Cov1 <- Cov2 <- Phi <- NULL
  if (is.null(cov_name)) {
    g <- Phi_df_long %>% 
      ggplot2::ggplot(ggplot2::aes(x = Class, y = Phi, group = Class, col = Class)) + 
      ggplot2::theme_bw() + 
      ggplot2::scale_color_brewer(palette = "Set2", labels = class_labels) + 
      ggplot2::labs(col = class_title, x = x_title, y = y_title) + 
      ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.5)) + 
      ggplot2::theme(text = ggplot2::element_text(size = 15), 
                     axis.text.x = ggplot2::element_text(size = 10, color = "black"), 
                     axis.text.y = ggplot2::element_text(size = 10, color = "black"), 
                     axis.title.x = ggplot2::element_text(size = 12, color = "black", face = "bold"), 
                     axis.title.y = ggplot2::element_text(size = 12, color = "black", face = "bold"), 
                     legend.title = ggplot2::element_text(size = 12, color = "black"), 
                     legend.text = ggplot2::element_text(size = 11, color = "black"), 
                     legend.position = "top")   
  } else {
    Phi_df_long$Cov1 <- factor(Phi_df_long$Cov1, labels = cov_labels[[1]])
    g <- Phi_df_long %>% 
      ggplot2::ggplot(ggplot2::aes(x = Cov1, y = Phi, group = Class, col = Class)) + 
      ggplot2::theme_bw() + 
      ggplot2::scale_color_brewer(palette = "Set2", labels = class_labels) + 
      ggplot2::labs(col = class_title, x = x_title, y = y_title) + 
      ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.5)) + 
      ggplot2::theme(text = ggplot2::element_text(size = 15), 
                     axis.text.x = ggplot2::element_text(size = 10, color = "black"), 
                     axis.text.y = ggplot2::element_text(size = 10, color = "black"), 
                     axis.title.x = ggplot2::element_text(size = 12, color = "black", face = "bold"), 
                     axis.title.y = ggplot2::element_text(size = 12, color = "black", face = "bold"), 
                     legend.title = ggplot2::element_text(size = 12, color = "black"), 
                     legend.text = ggplot2::element_text(size = 11, color = "black"), 
                     legend.position = "top")
    
    if (length(cov_name) == 2) {
      Phi_df_long$Cov2 <- factor(Phi_df_long$Cov2, labels = cov_labels[[2]])
      g <- g + ggplot2::facet_grid(~Cov2)
    }
  }
  
  if (!is.null(ci_level)) {
    g <- g + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Phi_lb, ymax = Phi_ub, col = Class), 
      width = 0.5, alpha = 0.7, position = ggplot2::position_dodge(width = 0.5))
  }
  if (add_lines) {
    g <- g + ggplot2::geom_line(linewidth = 0.7, alpha = 0.3, 
                                position = ggplot2::position_dodge(width = 0.5))
  }
  return(g)
}


plot_pattern_probs <- function (res, item_labels = NULL, categ_labels = NULL, categ_title = "Consumption Level", 
                                class_labels = NULL, class_title = "Dietary Pattern", y_title = "Consumption Level Probability", 
                                num_rows = 4,
                                ...) {
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting \n         from a call to one of these functions")
  }
  if (!is.null(res$estimates_adjust)) {
    est_item_probs <- res$estimates_adjust$theta_med
  }
  else {
    est_item_probs <- res$estimates$theta_med
  }
  K <- dim(est_item_probs)[2]
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  }
  else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ", 
                res$data_vars$J))
  }
  if (is.null(class_labels)) {
    class_labels <- 1:K
  }
  else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", 
                K))
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  }
  else if (length(categ_labels) != res$data_vars$R) {
    stop(paste0("length of categ_labels must equal the number of exposure categories, R = ", 
                res$data_vars$R))
  }
  dimnames(est_item_probs)[[1]] <- item_labels
  dimnames(est_item_probs)[[2]] <- class_labels
  theta_plot <- data.frame(expand.grid(lapply(dim(est_item_probs), seq_len)), 
                           value = as.vector(est_item_probs))
  Item <- Class <- Probability <- Level <- NULL
  colnames(theta_plot) <- c("Item", "Class", "Level", "Probability")
  theta_plot %>% 
    ggplot2::ggplot(ggplot2::aes(x = factor(Class, labels = class_labels), 
                                 y = Probability, fill = factor(Level))) + 
    ggplot2::geom_bar(stat = "identity", position = "stack") + 
    ggplot2::facet_wrap(factor(Item, labels = item_labels) ~ ., nrow = num_rows) + 
    ggplot2::scale_fill_brewer(type = "seq", palette = "RdYlBu", 
                               direction = -1, name = categ_title, labels = categ_labels) + 
    ggplot2::theme_bw() + ggplot2::labs(x = class_title, y = y_title) + 
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_text(size = 11, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"), 
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"), 
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"), 
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"), 
                   legend.text = ggplot2::element_text(size = 11, color = "black"), 
                   legend.position = "top", strip.text = ggplot2::element_text(size = 9), 
                   strip.background = ggplot2::element_rect(fill = "gray90"))
}


reorder_classes <- function (res, new_order) {
  if (!inherits(res, c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting \n         from a call to one of these functions")
  }
  else if ((inherits(res, "wolca")) & !is.null(res$estimates_svyglm)) {
    warning(paste0("For WOLCA, reordering of classes should be done before ", 
                   "calling wolca_svyglm(). res$estimates_svyglm should be NULL ", 
                   "prior to running this function."))
  }
  res_new <- res
  if (!is.null(res$estimates_adjust)) {
    ### CHANGED
    res_new$estimates_adjust$pi_red <- res$estimates_adjust$pi_red[, 
                                                                   new_order, drop = FALSE]
    res_new$estimates_adjust$theta_red <- res$estimates_adjust$theta_red[, 
                                                                         , new_order, , drop = FALSE]
    res_new$estimates_adjust$pi_med <- res$estimates_adjust$pi_med[new_order, drop = FALSE]
    res_new$estimates_adjust$theta_med <- res$estimates_adjust$theta_med[, 
                                                                         new_order, , drop = FALSE]
    ### END CHANGED
    for (i in 1:5) {
      res_new$estimates_adjust$c_all[res$estimates_adjust$c_all == 
                                       new_order[i]] <- i
    }
    if (is(res, "swolca")) {
      ### CHANGED
      res_new$estimates_adjust$xi_red <- res$estimates_adjust$xi_red[, 
                                                                     new_order, , drop = FALSE]
      res_new$estimates_adjust$xi_med <- res$estimates_adjust$xi_med[new_order, , drop = FALSE]
    }
  }
  ### CHANGED
  res_new$estimates$pi_red <- res$estimates$pi_red[, new_order, drop = FALSE]
  res_new$estimates$theta_red <- res$estimates$theta_red[, 
                                                         , new_order, , drop = FALSE]
  res_new$estimates$pi_med <- res$estimates$pi_med[new_order, drop = FALSE]
  res_new$estimates$theta_med <- res$estimates$theta_med[, 
                                                         new_order, , drop = FALSE]
  for (i in 1:5) {
    res_new$estimates$c_all[res$estimates$c_all == new_order[i]] <- i
  }
  if (is(res, "swolca")) {
    res_new$estimates$xi_red <- res$estimates$xi_red[, 
                                                     new_order, , drop = FALSE]
    res_new$estimates$xi_med <- res$estimates$xi_med[new_order, , drop = FALSE]
    ### END CHANGED
  }
  return(res_new)
}


plot_class_dist <- function (res, class_labels = NULL, class_title = "Dietary Pattern", 
          y_title = "Class Membership Probability", ...) 
{
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting \n         from a call to one of these functions")
  }
  ### CHANGED
  if (!is.null(res$estimates_adjust)) {
    pi_red <- as.data.frame(res$estimates_adjust$pi_red)
  }
  ### END CHANGED
  else {
    pi_red <- as.data.frame(res$estimates$pi_red)
  }
  K <- dim(pi_red)[2]
  if (is.null(class_labels)) {
    class_labels <- 1:K
  }
  else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", 
                K))
  }
  colnames(pi_red) <- class_labels
  pi_red_plot <- pi_red %>% tidyr::pivot_longer(cols = tidyselect::everything(), 
                                                names_to = "pi_comp", values_to = "value")
  pi_comp <- value <- NULL
  pi_red_plot %>% ggplot2::ggplot(ggplot2::aes(x = pi_comp, 
                                               y = value)) + ggplot2::theme_bw() + ggplot2::scale_fill_brewer(palette = "Set2") + 
    ggplot2::geom_boxplot() + ggplot2::labs(x = class_title, 
                                            y = y_title) + ggplot2::theme(text = ggplot2::element_text(size = 15), 
                                                                          axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                                                                          axis.text.y = ggplot2::element_text(size = 11, color = "black"), 
                                                                          axis.title.x = ggplot2::element_text(size = 13, color = "black", 
                                                                                                               face = "bold"), axis.title.y = ggplot2::element_text(size = 13, 
                                                                                                                                                                    color = "black", face = "bold"))
}

plot_regr_coefs <- function (regr_coefs, res, cov_labels = NULL, ...) {
  ### CHANGED
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop(paste0("res must be an object of class `swolca` or `wolca`, resulting ",
    "from a call to one of these functions"))
  }
  if (regr_coefs$Covariate[2] != "c_all2") {
    stop(paste0("regr_coefs must have latent class covariates with names of the ",
    "form `c_allk`, where k ranges from 2 to the total number of classes"))
  }
  if (is.null(cov_labels)) {
    cov_labels <- levels(factor(regr_coefs$Covariate, 
                                levels = regr_coefs$Covariate))
  }
  else {
    if ((length(cov_labels) != length(regr_coefs$Covariate)) | 
        !is.character(cov_labels)) {
      stop(paste0("cov_labels must be a string vector that has the same values as ",
                  "`regr_coefs$Covariate` and is of the same length"))
    } else if (length(setdiff(cov_labels, regr_coefs$Covariate))) {
      stop(paste0("cov_labels must be a string vector that has the same values as  ",
                  "`regr_coefs$Covariate` and is of the same length"))
    }
  }
  ### END CHANGED
  if (!is.null(res$estimates_adjust)) {
    K <- length(res$estimates_adjust$pi_med)
  }
  else {
    K <- length(res$estimates$pi_med)
  }
  plot_df <- regr_coefs
  plot_df$Class <- 1
  for (k in 1:K) {
    class_str <- paste0("c_all", k)
    plot_df$Class[grepl(class_str, plot_df$Covariate)] <- k
  }
  plot_df$Class <- as.factor(plot_df$Class)
  ### CHANGED
  plot_df$Covariate <- factor(plot_df$Covariate, levels = cov_labels, 
                              labels = cov_labels)
  ### END CHANGED
  Class <- Covariate <- Estimate <- LB <- UB <- NULL
  plot_df %>% ggplot2::ggplot(ggplot2::aes(x = Covariate, 
                                           y = Estimate, col = Class)) + ggplot2::theme_bw() + 
    ggplot2::geom_point() + ggplot2::scale_color_brewer(palette = "Set2", 
                                                        labels = Class, aesthetics = c("color")) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 70, 
                                                                                                                                                      vjust = 1, hjust = 1)) + ggplot2::geom_hline(yintercept = 0, 
                                                                                                                                                                                                   linetype = "dashed") + ggplot2::geom_errorbar(ggplot2::aes(ymin = LB, 
                                                                                                                                                                                                                                                              ymax = UB, color = Class, width = 0.2))
}

#================ OLD CODE =====================================================
# # Exploration of dataset variables
# # Large population centers: 
# Bronx-DR (1148)
# Bronx-PR (1449)
# [Chicago-C/SA (647)]
# [Chicago-PR (674)]
# Chicago-Mex (2006) 
# Miami-C/SA (1357)
# Miami-Cuban (2089)
# San Diego-Mex (3212)
# Total: 6 groups = 11261; 8 groups = 12582

# For Mexicans, Bronx: 163, Chicago: 2006, Miami: 29, San Diego: 3212
# # Subset to Mexican background for Chicago and San Diego: total 5218
# fpq_mex <- fpq_data %>%
#   filter(CENTER %in% c("C", "S") & BKGRD1_C6 == 3)
#
#
# table(raw_derv_vars$HYPERTENSION2, useNA = "always")
# table(raw_derv_vars$HYPERT_AWARENESS, useNA = "always")
# table(raw_derv_vars$HYPERTENSION2, raw_derv_vars$HYPERT_AWARENESS, 
#       useNA = "always")
# 
# table(raw_derv_vars$DYSLIPIDEMIA, useNA = "always")
# 
# table(raw_derv_vars$DIABETES2, useNA = "always")
# 
# par(mfrow = c(1,1))
# par(mar=c(5,4,4,2)+0.1)
# hist(raw_derv_vars$YRSUS)
# hist(raw_derv_vars$AGE_IMMI)
# table(raw_derv_vars$US_BORN)
# yrsus_3cat <- raw_derv_vars %>% mutate(
#   YRSUS_C3 = case_when(
#     US_BORN == 1 ~ 1, 
#     US_BORN == 0 & YRSUS >= 10 ~ 2,
#     US_BORN == 0 & YRSUS < 10 ~ 3,
#     .default = NA
#   ))
# table(yrsus_3cat$YRSUS_C3, useNA = "always")


## Previous analyses for age x gender
# # Combined study data: 12519 x 394
# study_data <- raw_derv_sub %>% 
#   mutate(ID = as.integer(ID)) %>%
#   left_join(fpq_data %>% 
#               select(ID, fp2:fp1) %>% 
#               mutate(across(fp2:fp1, as.numeric)), by = join_by(ID))
# 
# study_data <- study_data %>%
#   mutate(subgroup = case_when(
#     AGE < 45 & GENDER == "F" ~ 1, #"f_young",
#     AGE >= 45 & AGE < 65 & GENDER == "F" ~ 2, #"f_mid",
#     AGE >= 54 & GENDER == "F" ~ 3, #"f_old",
#     AGE < 45 & GENDER == "M" ~ 4, #"m_young",
#     AGE >= 45 & AGE < 65 & GENDER == "M" ~ 5, #"m_mid",
#     AGE >= 54 & GENDER == "M" ~ 6, #"m_old",
#     .default = NA
# ))

### Define subgroups
# study_data <- study_data %>%
#   mutate(subgroup = case_when(
#     AGE < 45 & GENDER == "F" ~ 1, #"f_young",
#     AGE >= 45 & AGE < 65 & GENDER == "F" ~ 2, #"f_mid",
#     AGE >= 54 & GENDER == "F" ~ 3, #"f_old",
#     AGE < 45 & GENDER == "M" ~ 4, #"m_young",
#     AGE >= 45 & AGE < 65 & GENDER == "M" ~ 5, #"m_mid",
#     AGE >= 54 & GENDER == "M" ~ 6, #"m_old",
#     .default = NA
#   ))


# # Check missingness
# naniar::gg_miss_var(study_data %>% 
#                       select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
#                              subgroup, fp2:fp1))
# # naniar::vis_miss(study_data %>% 
# #                    select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
# #                           subgroup, fp2:fp1))
# 
# # Drop missing: sample size is 11519 x 394
# # Note: no missingness for ID, STRAT, PSU_ID, WEIGHT
# # DIABETES2 has 3 missing
# # fp2:fp1 has ~1000 missing
# study_data <- study_data %>%
#   drop_na(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, subgroup, fp2:fp1)
# 
# # Distribution of various subgroups
# table(study_data$subgroup)
# # Distribution of various subgroups by site/ethnicity
# table(study_data$site_eth, study_data$subgroup)