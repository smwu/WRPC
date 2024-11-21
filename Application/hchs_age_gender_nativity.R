#====================================================
# Analyzing HCHS data for Carribean background
# Patterns by age, gender, and nativity
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


#================= FPQ cleaning ================================================

# Select diet variables
fpq_data <- fpq_49categ %>% 
  select(ID, citrus_juice:soup_oth)

# Drop NAs: n=12738
fpq_dropna <- fpq_data %>% drop_na()

# Merge in derived variables, which include demographic and outcome variables.
# Merge in self-report variables
# Combined study data: 12761 x 521
study_data <- raw_derv_vars %>% 
  mutate(ID = as.integer(ID)) %>%
  left_join(raw_mhea_vars %>% mutate(ID = as.integer(ID)), by = join_by(ID)) %>%
  right_join(fpq_dropna %>% 
               mutate(across(citrus_juice:soup_oth, as.numeric)), by = join_by(ID)) 

# Check all 5-level foods
apply(study_data %>% select(citrus_juice:soup_oth), 2, 
      function(x) max(x, na.rm = TRUE))


### Check site x background numbers
# Site: B=Bronx, C=Chicago, M=Miami, S=San Diego
# Background: 0=Dominican, 1=Central or South American, 2=Cuban, 3=Mexican,
# 4=Puerto Rican, 5=More than one/Other heritage
table(study_data$CENTER, study_data$BKGRD1_C6)


#============ Variable conversions

# Variable conversions
study_data <- study_data %>%
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
                             .default = NA))

# Drop NAs
# Sample size after dropping all missingness: n = 12674
study_data_dropna <- study_data %>% 
  select(ID, STRAT, PSU_ID, WEIGHT_FINAL_EXPANDED, CENTER, BKGRD1_C6,
         AGE_C3, GENDER, YRSUS_C3, citrus_juice:soup_oth) %>% 
  drop_na()

### Subset to Carribean backgrounds (Dominican, Cuban, and Puerto Rican)
# n = 5168
study_data_dropna <- study_data_dropna %>%
  filter(BKGRD1_C6 %in% c(0, 2, 4))

### Create subgroup based on age and nativity
study_data_dropna <- study_data_dropna %>%
  mutate(subgroup = case_when(
    AGE_C3 == 1 & YRSUS_C3 == 1 ~ 1,  # young : us-born
    AGE_C3 == 1 & YRSUS_C3 == 2 ~ 2,  # young : accult
    AGE_C3 == 1 & YRSUS_C3 == 3 ~ 3,  # young : unaccult
    AGE_C3 == 2 & YRSUS_C3 == 1 ~ 4,  # mid : us-born
    AGE_C3 == 2 & YRSUS_C3 == 2 ~ 5,  # mid : accult
    AGE_C3 == 2 & YRSUS_C3 == 3 ~ 6,  # mid : unaccult
    AGE_C3 == 3 & YRSUS_C3 %in% c(1, 2) ~ 7,  # old : us-born or accult (us-born only 8)
    AGE_C3 == 3 & YRSUS_C3 == 3 ~ 8,  # old : unaccult
    .default = NA
  ))
table(study_data_dropna$subgroup)

#================= WRPC model with subgroups ===================================

### Run WRPC with subgroup 

# Categorical exposure matrix, nxJ, 5168 x 49
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
                 burn = 5000, thin = 5, update = 1000, switch = 50, save_res = TRUE, 
                 save_path = paste0(wd, res_dir, "HCHS_carrib_age_yrsus"))
res_wrpc$runtime # 3 hours



#======================
### Subset to Puerto Rican background
# n = 5168
study_data_dropna <- study_data_dropna %>%
  filter(BKGRD1_C6 %in% 4)