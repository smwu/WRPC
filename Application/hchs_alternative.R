#====================================================
# Alternative HCHS analyses
# Author: Stephanie Wu
# Date Updated: 2024/11/26
#====================================================

rm(list = ls())

library(tidyverse)  # data wrangling
library(baysc)      # bayesian survey clustering
library(haven)      # read sas file
library(naniar)     # missingness
library(sfsmisc)    # mult.fig

wd <- "~/Documents/Github/WRPC/"
wd <- "/n/netscratch/stephenson_lab/Lab/stephwu18/WRPC/"
code_dir <- "Model_Code/"
data_dir <- "Application/HCHS_Data/"
res_dir <- "Results/"

# Read in raw derived variables
raw_derv_vars <- read_sas(paste0(wd, data_dir, "part_derv_inv4.sas7bdat"))
# Read in raw self-report variables
raw_mhea_vars <- read_sas(paste0(wd, data_dir, "mhea_inv4.sas7bdat"))
# Read in raw FPQ data
raw_fpq_vars <- read_sas(paste0(wd, data_dir, "fpe_inv4.sas7bdat"))
# Read in cleaned FPQ data w/ 49 categories
fpq_49categ <- read.csv(paste0(wd, data_dir, "fpq_49categ_2024Nov19.csv"))

# Source WRPC functions
source(paste0(wd, code_dir, "wrpc.R"))
source(paste0(wd, code_dir, "wrpc_mcmc_fns.R"))
source(paste0(wd, code_dir, "wrpc_utilities.R"))
Rcpp::sourceCpp(paste0(wd, code_dir, "wrpc_mcmc.cpp"))


#================= FPQ cleaning ================================================

# Select diet variables
fpq_data <- fpq_49categ %>% 
  select(ID, citrus_juice:soup_oth)

# Drop NAs: n=12761
fpq_dropna <- fpq_data %>% drop_na()

# Merge in derived variables, which include demographic and outcome variables.
# Merge in self-report variables
# Combined study data: 12761 x 441
comb_data <- raw_derv_vars %>% 
  mutate(ID = as.integer(ID)) %>%
  left_join(raw_mhea_vars %>% mutate(ID = as.integer(ID)), by = join_by(ID)) %>%
  right_join(fpq_dropna %>% 
               mutate(across(citrus_juice:soup_oth, as.numeric)), by = join_by(ID)) 

# Check all 5-level foods
apply(comb_data %>% select(citrus_juice:soup_oth), 2, 
      function(x) max(x, na.rm = TRUE))


#============ Variable conversions =============================================

# Variable conversions
comb_data <- comb_data %>%
  mutate(# Years in US categorized,
    YRSUS_C3 = case_when(
      US_BORN == 1 ~ 1, # "US_born",
      US_BORN == 0 & YRSUS >= 10 ~ 2, # ">=10",
      US_BORN == 0 & YRSUS < 10 ~ 3, # "<10",
      .default = NA
    ),
    # Age categories
    AGE_C3 = case_when(
      AGE < 45 ~ 1, # "<45", 
      (AGE >= 45) & (AGE < 65) ~ 2, # "45-64",
      AGE >= 65 ~ 3, # ">=65",
      .default = NA),
    GENDER = case_match(GENDER, 
                        "M" ~ 1, 
                        "F" ~ 2,
                        .default = NA),
    ETHIS_C3 = case_when(
      ETHIS < 3 ~ 1, # low
      ETHIS == 3 ~ 2, # medium
      ETHIS > 3 ~ 3, # high
      .default = NA
    ),
    ETHIS_C2 = case_when(
      ETHIS <= 2.5 ~ 1, # low
      ETHIS > 2.5 ~ 2, # high
      .default = NA
    ),
    INCOME_C2 = case_when(
      INCOME %in% c(1, 2, 3, 4, 5) ~ 1, # <30k annual household
      INCOME %in% c(6, 7, 8, 9, 10) ~ 2, # >=30k
      .default = NA
    ),
    INCOME_NEW_C3 = case_when(
      INCOME %in% c(1, 2, 3) ~ 1, # <= 20k annual household
      INCOME %in% c(4, 5, 6) ~ 2, # >20k and <= 40k
      INCOME %in% c(7, 8, 9, 10) ~ 3, # >40k
      .default = NA
    ))

# Note: can't create poverty thresholds b/c don't have ECEA5 variable


#================== Carribean, MARITAL-GENDER-AGE ==============================

# Subset to Carribean background: n = 5209
study_data <- comb_data %>%
  filter(BKGRD1_C6 %in% c(0, 2, 4))

# Exploratory proportions
# Marital: 1=single, 2=married/cohabiting, 3=separated/divorced/widowed
# Gender: 1=M, 2=F
# Age: 1=18-44, 2=45+
table(study_data$MARITAL_STATUS, study_data$GENDER)  
table(study_data$MARITAL_STATUS, study_data$AGEGROUP_C2) 

### Create subgroup based on marital status, gender, age
subgroup_labels <- 
  c("single-M-young", "single-M-old", "single-F-young", "single-F-old",
    "married-M-young", "married-M-old", "married-F-young", "married-F-old",
    "divorced-M-young", "divorced-M-old", "divorced-F-young", "divorced-F-old")
study_data <- study_data %>%
  mutate(subgroup_marital_gender_age = case_when(
    MARITAL_STATUS == 1 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 1,  # single-M-young
    MARITAL_STATUS == 1 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 2,  # single-M-old
    MARITAL_STATUS == 1 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 3,  # single-F-young
    MARITAL_STATUS == 1 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 4,  # single-F-old
    MARITAL_STATUS == 2 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 5,  # married-M-young
    MARITAL_STATUS == 2 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 6,  # married-M-old
    MARITAL_STATUS == 2 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 7,  # married-F-young
    MARITAL_STATUS == 2 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 8,  # married-F-old
    MARITAL_STATUS == 3 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 9,  # divorced-M-young
    MARITAL_STATUS == 3 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 10, # divorced-M-old
    MARITAL_STATUS == 3 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 11, # divorced-F-young
    MARITAL_STATUS == 3 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 12, # divorced-F-old
    .default = NA
  ))
table(study_data$subgroup_marital_gender_age)

# Check missingness
# Note: FPQ is already subsetted to complete cases
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             MARITAL_STATUS, GENDER, AGEGROUP_C2, 
                             subgroup_marital_gender_age))

# Drop NAs
# Sample size after dropping all missingness: n = 5203
study_data_dropna <- study_data %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
         subgroup_marital_gender_age, citrus_juice:soup_oth) %>% 
  drop_na()


### Run WRPC with subgroup 

# Categorical exposure matrix, nxJ, 5203 x 49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup_marital_gender_age))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                 burn = 5000, thin = 5, update = 1000, switch = 50, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_marital_age_gender"))
res_wrpc$runtime # 3.6 hours, K_fixed = 6, K_red = 9


#================== Puerto Rican, NATIVITY-ETHIS-AGE ==============================

# Subset to Puerto Rican background: n = 2066
study_data <- comb_data %>%
  filter(BKGRD1_C6 %in% c(4))

# Exploratory proportions
# YRSUS_C3: 1=US-born, 2= >=10 years, 3= <10 years (very few in 3)
# US_BORN: 0=no, 1=yes (born in 50 US states/continental)
# ETHIS_C3: 1=low, 2=median, 3=high
# Age: 1=18-44, 2=45+
table(study_data$AGEGROUP_C2, study_data$YRSUS_C3)  
table(study_data$AGEGROUP_C2, study_data$ETHIS_C3) 
table(study_data$YRSUS_C3, study_data$ETHIS_C3)
table(study_data$AGEGROUP_C2, study_data$US_BORN) 
table(study_data$US_BORN, study_data$ETHIS_C3)

### Create subgroup based on age, nativity, ethnic identity score
subgroup_labels <- 
  c("nonUS-ethL-young", "nonUS-ethL-old", "nonUS-ethM-young", "nonUS-ethM-old", 
    "nonUS-ethH-young", "nonUS-ethH-old",
    "US-ethL-young", "US-ethL-old", "US-ethM-young", "US-ethM-old",
    "US-ethH-young", "US-ethH-old")
study_data <- study_data %>%
  mutate(subgroup_nativity_ethis_age = case_when(
    US_BORN == 0 & ETHIS_C3 == 1 & AGEGROUP_C2 == 1 ~ 1,  # nonUS-ethL-young
    US_BORN == 0 & ETHIS_C3 == 1 & AGEGROUP_C2 == 2 ~ 2,  # nonUS-ethL-old
    US_BORN == 0 & ETHIS_C3 == 2 & AGEGROUP_C2 == 1 ~ 3,  # nonUS-ethM-young
    US_BORN == 0 & ETHIS_C3 == 2 & AGEGROUP_C2 == 2 ~ 4,  # nonUS-ethM-old
    US_BORN == 0 & ETHIS_C3 == 3 & AGEGROUP_C2 == 1 ~ 5,  # nonUS-ethH-young
    US_BORN == 0 & ETHIS_C3 == 3 & AGEGROUP_C2 == 2 ~ 6,  # nonUS-ethH-old
    US_BORN == 1 & ETHIS_C3 == 1 & AGEGROUP_C2 == 1 ~ 7,  # US-ethL-young
    US_BORN == 1 & ETHIS_C3 == 1 & AGEGROUP_C2 == 2 ~ 8,  # US-ethL-old
    US_BORN == 1 & ETHIS_C3 == 2 & AGEGROUP_C2 == 1 ~ 9,  # US-ethM-young
    US_BORN == 1 & ETHIS_C3 == 2 & AGEGROUP_C2 == 2 ~ 10, # US-ethM-old
    US_BORN == 1 & ETHIS_C3 == 3 & AGEGROUP_C2 == 1 ~ 11, # US-ethH-young
    US_BORN == 1 & ETHIS_C3 == 3 & AGEGROUP_C2 == 2 ~ 12, # US-ethH-old
    .default = NA
  ))
table(study_data$subgroup_nativity_ethis_age)
round(prop.table(table(study_data$subgroup_nativity_ethis_age)), 2)


# Check missingness
# Note: FPQ is already subsetted to complete cases
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             US_BORN, ETHIS_C3, AGEGROUP_C2, 
                             subgroup_nativity_ethis_age))

# Drop NAs
# Sample size after dropping all missingness: n = 2042
study_data_dropna <- study_data %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
         subgroup_nativity_ethis_age, citrus_juice:soup_oth) %>% 
  drop_na()


### Run WRPC with subgroup 

# Categorical exposure matrix, nxJ, 2042 x 49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup_nativity_ethis_age))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                 burn = 5000, thin = 5, update = 1000, switch = 50, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_nativity_ethis_age"))
res_wrpc$runtime # 1.5 hours


#================== Puerto Rican, INCOME-NATIVITY-GENDER-AGE ===================

# Subset to Puerto Rican background: n = 2066
study_data <- comb_data %>%
  filter(BKGRD1_C6 %in% c(4))

# Exploratory proportions
# Income: 1= <30k, 2 = >=30k
# US_BORN: 0=no, 1=yes (born in 50 US states/continental)
# Gender: 1=M, 2=F
# Age: 1=18-44, 2=45+
table(study_data$INCOME_C5)
table(study_data$AGE_C3)
table(study_data$INCOME_C2, study_data$AGEGROUP_C2)
table(study_data$INCOME_C2, study_data$GENDER)  
table(study_data$INCOME_C2, study_data$US_BORN) 

### Create subgroup based on income, nativity, gender, age
subgroup_labels <- 
  c("incL-nonUS-M-young", "incL-nonUS-M-old", "incL-nonUS-F-young", "incL-nonUS-F-old",
    "incL-US-M-young", "incL-US-M-old", "incL-US-F-young", "incL-US-F-old",
    "incH-nonUS-M-young", "incH-nonUS-M-old", "incH-nonUS-F-young", "incH-nonUS-F-old", 
    "incH-US-M-young", "incH-US-M-old", "incH-US-F-young", "incH-US-F-old")
study_data <- study_data %>%
  mutate(subgroup_income_nativity_gender_age = case_when(
    INCOME_C2 == 1 & US_BORN == 0 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 1,  # incL-nonUS-M-young
    INCOME_C2 == 1 & US_BORN == 0 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 2,  # incL-nonUS-M-old
    INCOME_C2 == 1 & US_BORN == 0 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 3,  # incL-nonUS-F-young
    INCOME_C2 == 1 & US_BORN == 0 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 4,  # incL-nonUS-F-old
    INCOME_C2 == 1 & US_BORN == 1 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 5,  # incL-US-M-young
    INCOME_C2 == 1 & US_BORN == 1 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 6,  # incL-US-M-old
    INCOME_C2 == 1 & US_BORN == 1 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 7,  # incL-US-F-young
    INCOME_C2 == 1 & US_BORN == 1 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 8,  # incL-US-F-old
    INCOME_C2 == 2 & US_BORN == 0 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 9,  # incH-nonUS-M-young
    INCOME_C2 == 2 & US_BORN == 0 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 10, # incH-nonUS-M-old
    INCOME_C2 == 2 & US_BORN == 0 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 11, # incH-nonUS-F-young
    INCOME_C2 == 2 & US_BORN == 0 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 12, # incH-nonUS-F-old
    INCOME_C2 == 2 & US_BORN == 1 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 13,  # incH-US-M-young
    INCOME_C2 == 2 & US_BORN == 1 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 14, # incH-US-M-old
    INCOME_C2 == 2 & US_BORN == 1 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 15, # incH-US-F-young
    INCOME_C2 == 2 & US_BORN == 1 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 16, # incH-US-F-old
    .default = NA
  ))
table(study_data$subgroup_income_nativity_gender_age)

# Check missingness
# Note: FPQ is already subsetted to complete cases
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             INCOME_C2, US_BORN, GENDER, AGEGROUP_C2, 
                             subgroup_income_nativity_gender_age))

# Drop NAs
# Sample size after dropping all missingness: n = 1896
study_data_dropna <- study_data %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
         subgroup_income_nativity_gender_age, citrus_juice:soup_oth) %>% 
  drop_na()


### Run WRPC with subgroup 

# Categorical exposure matrix, nxJ, 1896 x 49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
dim(x_mat)
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup_income_nativity_gender_age))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                 burn = 5000, thin = 5, update = 1000, switch = 50, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_income_nativity_gender_age"))
res_wrpc$runtime # 1.5 hours, K = 4



### Run Unweighted WRPC with subgroup 

seed <- 1
res_unwt_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = NULL, cluster_id = NULL, 
                 stratum_id = NULL, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                 burn = 5000, thin = 5, update = 1000, switch = 50, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_income_nativity_gender_age_unwt"))
res_unwt_wrpc$runtime # 1.5 hours, K = 4
plot(res_unwt_wrpc$post_MCMC_out$dendrogram_global)
# Average posterior probability: 0.94
mean(apply(res_unwt_wrpc$estimates$pred_global_class_probs, 1, max))


### Run WOLCA without subgroup 

seed <- 1
res_wolca <- baysc::wolca(x_mat = x_mat, 
                          sampling_wt = sampling_wt, cluster_id = cluster_id, 
                          stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff = 0.05, n_runs = 10000, 
                 burn = 5000, thin = 5, update = 1000, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_income_nativity_gender_age"))
res_wolca$runtime # 40 mins, K = 3
plot(res_wolca$post_MCMC_out$dendrogram)
# Average posterior probability: 0.94
mean(apply(res_wolca$estimates$pred_class_probs, 1, max))
item_labels <- colnames(fpq_49categ[, -c(1:3)])
categ_labels <- c("None", "Monthly", "Weekly", "Daily", "Daily+")
baysc::plot_pattern_profiles(res = res_wolca, item_labels = item_labels, 
                             categ_labels = categ_labels)


#================== Puerto Rican, INCOME (2-cat)-GENDER-AGE ====================

# Subset to Puerto Rican background: n = 2066
study_data <- comb_data %>%
  filter(BKGRD1_C6 %in% c(4))

# Exploratory proportions
# Income: 1= <30k, 2 = >=30k
# Gender: 1=M, 2=F
# Age: 1=18-44, 2=45+
table(study_data$INCOME_C2, study_data$AGEGROUP_C2)
table(study_data$INCOME_C2, study_data$GENDER)
table(study_data$AGEGROUP_C2, study_data$GENDER)

### Create subgroup based on income, gender, age
subgroup_labels <- 
  c("incL-M-young", "incL-M-old", "incL-F-young", "incL-F-old",
    "incH-M-young", "incH-M-old", "incH-F-young", "incH-F-old")
study_data <- study_data %>%
  mutate(subgroup_incomeC2_gender_age = case_when(
    INCOME_C2 == 1 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 1,  # incL-M-young
    INCOME_C2 == 1 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 2,  # incL-M-old
    INCOME_C2 == 1 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 3,  # incL-F-young
    INCOME_C2 == 1 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 4,  # incL-F-old
    INCOME_C2 == 2 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 5,  # incH-M-young
    INCOME_C2 == 2 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 6,  # incH-M-old
    INCOME_C2 == 2 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 7,  # incH-F-young
    INCOME_C2 == 2 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 8, # incH-F-old
    .default = NA
  ))
table(study_data$subgroup_incomeC2_gender_age)

# Check missingness
# Note: FPQ is already subsetted to complete cases
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             INCOME_NEW_C3, GENDER, AGEGROUP_C2, 
                             subgroup_incomeC2_gender_age))

# Drop NAs
# Sample size after dropping all missingness: n = 1896
study_data_dropna <- study_data %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
         subgroup_incomeC2_gender_age, citrus_juice:soup_oth) %>% 
  drop_na()


### Run WRPC with subgroup 

# Categorical exposure matrix, nxJ, 1896 x 49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
dim(x_mat)
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup_incomeC2_gender_age))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                 burn = 5000, thin = 5, update = 1000, switch = 50, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_incomeC2_gender_age"))
res_wrpc$runtime # 1.3 hours, K = 5



#================== Puerto Rican, INCOME-NATIVITY-AGE AMONG WOMEN ==============

# Subset to Puerto Rican background and women: n = 1258
study_data <- comb_data %>%
  filter(BKGRD1_C6 %in% c(4) & GENDER == 2)

# Exploratory proportions
# Income: 1= <30k, 2 = >=30k
# US_BORN: 0=no, 1=yes (born in 50 US states/continental)
# Age: 1=18-44, 2=45+
table(study_data$INCOME_C2, study_data$AGEGROUP_C2)
table(study_data$INCOME_C2, study_data$US_BORN) 

### Create subgroup based on income, nativity, gender, age
subgroup_labels <- 
  c("incL-nonUS-F-young", "incL-nonUS-F-old", "incL-US-F-young", "incL-US-F-old",
    "incH-nonUS-F-young", "incH-nonUS-F-old", "incH-US-F-young", "incH-US-F-old")
study_data <- study_data %>%
  mutate(subgroup_income_nativity_age = case_when(
    INCOME_C2 == 1 & US_BORN == 0 & AGEGROUP_C2 == 1 ~ 1,  # incL-nonUS-F-young
    INCOME_C2 == 1 & US_BORN == 0 & AGEGROUP_C2 == 2 ~ 2,  # incL-nonUS-F-old
    INCOME_C2 == 1 & US_BORN == 1 & AGEGROUP_C2 == 1 ~ 3,  # incL-US-F-young
    INCOME_C2 == 1 & US_BORN == 1 & AGEGROUP_C2 == 2 ~ 4,  # incL-US-F-old
    INCOME_C2 == 2 & US_BORN == 0 & AGEGROUP_C2 == 1 ~ 5, # incH-nonUS-F-young
    INCOME_C2 == 2 & US_BORN == 0 & AGEGROUP_C2 == 2 ~ 6, # incH-nonUS-F-old
    INCOME_C2 == 2 & US_BORN == 1 & AGEGROUP_C2 == 1 ~ 7, # incH-US-F-young
    INCOME_C2 == 2 & US_BORN == 1 & AGEGROUP_C2 == 2 ~ 8, # incH-US-F-old
    .default = NA
  ))
table(study_data$subgroup_income_nativity_age)

# Check missingness
# Note: FPQ is already subsetted to complete cases
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             INCOME_C2, US_BORN, GENDER, AGEGROUP_C2, 
                             subgroup_income_nativity_age))

# Drop NAs
# Sample size after dropping all missingness: n = 1150
study_data_dropna <- study_data %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
         subgroup_income_nativity_age, citrus_juice:soup_oth) %>% 
  drop_na()


### Run WRPC with subgroup 

# Categorical exposure matrix, nxJ, 1150 x 49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
dim(x_mat)
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup_income_nativity_age))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 20000, 
                 burn = 10000, thin = 5, update = 1000, switch = 50, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_women_income_nativity_age"))
res_wrpc$runtime # 1.5 hours, K = 3


#================== Puerto Rican, INCOME (3-cat)-GENDER-AGE ====================

# Subset to Puerto Rican background: n = 2066
study_data <- comb_data %>%
  filter(BKGRD1_C6 %in% c(4))

# Exploratory proportions
# Income: 1= <=20k, 2 = 20-40k, 3= >40k
# Gender: 1=M, 2=F
# Age: 1=18-44, 2=45+
table(study_data$AGE_C3)  # few people in group 3
table(study_data$INCOME_C2, study_data$AGEGROUP_C2)
table(study_data$INCOME_NEW_C3, study_data$AGEGROUP_C2)
table(study_data$INCOME_NEW_C3, study_data$AGE_C3) 
table(study_data$INCOME_NEW_C3, study_data$GENDER)

### Create subgroup based on income, gender, age
subgroup_labels <- 
  c("incL-M-young", "incL-M-old", "incL-F-young", "incL-F-old",
    "incM-M-young", "incM-M-old", "incM-F-young", "incM-F-old",
    "incH-M-young", "incH-M-old", "incH-F-young", "incH-F-old")
study_data <- study_data %>%
  mutate(subgroup_income_gender_age = case_when(
    INCOME_NEW_C3 == 1 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 1,  # incL-M-young
    INCOME_NEW_C3 == 1 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 2,  # incL-M-old
    INCOME_NEW_C3 == 1 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 3,  # incL-F-young
    INCOME_NEW_C3 == 1 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 4,  # incL-F-old
    INCOME_NEW_C3 == 2 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 5,  # incM-M-young
    INCOME_NEW_C3 == 2 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 6,  # incM-M-old
    INCOME_NEW_C3 == 2 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 7,  # incM-F-young
    INCOME_NEW_C3 == 2 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 8, # incM-F-old
    INCOME_NEW_C3 == 3 & GENDER == 1 & AGEGROUP_C2 == 1 ~ 9,  # incH-M-young
    INCOME_NEW_C3 == 3 & GENDER == 1 & AGEGROUP_C2 == 2 ~ 10,  # incH-M-old
    INCOME_NEW_C3 == 3 & GENDER == 2 & AGEGROUP_C2 == 1 ~ 11,  # incH-F-young
    INCOME_NEW_C3 == 3 & GENDER == 2 & AGEGROUP_C2 == 2 ~ 12, # incH-F-old
    .default = NA
  ))
table(study_data$subgroup_income_gender_age)

# Check missingness
# Note: FPQ is already subsetted to complete cases
naniar::gg_miss_var(study_data %>% 
                      select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
                             INCOME_NEW_C3, GENDER, AGEGROUP_C2, 
                             subgroup_income_gender_age))

# Drop NAs
# Sample size after dropping all missingness: n = 1896
study_data_dropna <- study_data %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, 
         subgroup_income_gender_age, citrus_juice:soup_oth) %>% 
  drop_na()


### Run WRPC with subgroup 

# Categorical exposure matrix, nxJ, 1896 x 49
x_mat <- as.matrix(study_data_dropna %>% select(citrus_juice:soup_oth)) 
dim(x_mat)
# Stratum indicators, nx1
stratum_id <- c(as.numeric(study_data_dropna$STRAT))      
# Cluster indicators, nx1
cluster_id <- c(as.numeric(study_data_dropna$PSU_ID))                   
# Survey sampling weights, nx1
sampling_wt <- c(as.numeric(study_data_dropna$WEIGHT_FINAL_EXPANDED))   
# Subgroup variable, nx1
h_all <- c(as.numeric(study_data_dropna$subgroup_income_gender_age))

seed <- 1
res_wrpc <- wrpc(x_mat = x_mat, h_all = h_all, 
                 sampling_wt = sampling_wt, cluster_id = cluster_id, 
                 stratum_id = stratum_id, run_sampler = "both", K_max = 20, 
                 adapt_seed = seed, class_cutoff_global = 0.05, n_runs = 10000, 
                 burn = 5000, thin = 5, update = 1000, switch = 50, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_income_gender_age"))
res_wrpc$runtime # 1.3 hours, K = 5



#=================== Plotting results ==========================================

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
plot_wrpc_global_pattern_probs(res = res_wrpc, item_labels = item_labels, 
                          categ_labels = categ_labels,
                          num_rows = 7) 
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(res_wrpc$post_MCMC_out$dendrogram_global)
# Average posterior probability
mean(apply(res_wrpc$estimates$pred_global_class_probs, 1, max))
plot_wrpc_local_pattern_profiles(res = res_wrpc, 
                                 item_labels = item_labels, item_title = "Item",
                                 categ_labels = categ_labels,
                                 subgroup_labels = subgroup_labels,
                                 subgroup_title = "Subgroup") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_wrpc_allocation(res = res_wrpc, item_labels = item_labels, 
                     item_title = "Item",
                     subgroup_labels = subgroup_labels, 
                     subgroup_title = "Subgroup") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_wrpc_local_profiles_allocation(res = res_wrpc, 
                                    item_labels = item_labels, item_title = "Item",
                                    categ_labels = categ_labels,
                                    subgroup_labels = subgroup_labels,
                                    subgroup_title = "Subgroup") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# Proportion of each subpopulation following each global pattern
table(res_wrpc$estimates$c_all, res_wrpc$data_vars$h_all)
plot_wrpc_class_subgroup_dist(res = res_wrpc, subgroup_labels = subgroup_labels, 
                         subgroup_title = "Subgroup", normalize = TRUE)
plot_wrpc_class_subgroup_dist(res = res_wrpc, subgroup_labels = subgroup_labels, 
                         subgroup_title = "Subgroup", normalize = FALSE)

## plot trace of pi
K <- res_wrpc$estimates$K_red
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for pi"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  plot(res_wrpc$estimates$pi_red[, k], type = "l", ylim = c(0,1), xlab = "", ylab = "")
}

# Plot trace of nu
h <- 1
mult.fig(mfrow = c(4, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for nu, h = ", h, ", j = 1:4"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (j in 1:4) {
  plot(res_wrpc$estimates$nu_red[, h, j], type = "l", ylim = c(0,1), xlab = "", ylab = "")
}
h <- 4
mult.fig(mfrow = c(4, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for nu, h = ", h, ", j = 1:4"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (j in 1:4) {
  plot(res_wrpc$estimates$nu_red[, h, j], type = "l", ylim = c(0,1), xlab = "", ylab = "")
}

# Plot trace of theta global
j = 1
k = 1
mult.fig(mfrow = c(5, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for theta_global, j=", j, ", k=", k, ", r=1:5"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (r in 1:5) {
  plot(res_wrpc$estimates$theta_global_red[, j, k, r], type = "l", 
       ylim = c(0,1), xlab = "", ylab = "")
}

# Plot trace of theta local
j = 1
h = 4
mult.fig(mfrow = c(5, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Trace plots for theta_local, j=", j, ", h=", h, ", r=1:5"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (r in 1:5) {
  plot(res_wrpc$estimates$theta_local_red[, j, h, r], type = "l", 
       ylim = c(0,1), xlab = "", ylab = "")
}

