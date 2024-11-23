# **WRPC**: Weighted Robust Profile Clustering

A repo for running the Bayesian weighted robust profile clustering model on survey data. The following folders are included:

- `Model_Code`: Code containing main functions for running the WRPC model. Contains subfolder `NIMBLE` that contains implementation of WRPC in NIMBLE [https://r-nimble.org](https://r-nimble.org).
- `Application`: Code for applying WRPC to data from the Hispanic Community Health Study/Study of Latinos (HCHS/SOL)
- `Simulated_Data`: Toy simulated data

**Maintainer**: Stephanie M. Wu ([swu\@g.harvard.edu](mailto:swu@g.harvard.edu){.email})

We provide an example dataset from the National Health and Nutrition Examination Survey (NHANES) that includes multivariate categorical dietary intake data, categorical age data, and binary hypertension data for low-income women in the United States. Survey sampling weights and information on stratification and clustering are included to allow for adjustment for survey design when conducting estimation and inference.

<div id='id-section1'/>

## Example

``` r
#==== Setup ====
# Load libraries
library(baysc)  # Install at https://github.com/smwu/baysc
library(tidyverse)

# Source functions
source("Model_code/wrpc.R")                  # Main WRPC function
source("Model_Code/wrpc_mcmc_fns.R")         # MCMC functions
Rcpp::sourceCpp("Model_Code/wrpc_mcmc.cpp")  # MCMC functions Rcpp code
source("Model_Code/wrpc_utilities.R")        # Additional utility functions
source("Model_Code/wrpc_plotting_fns.R")     # Plotting functions

#==== Load and prepare data ====
# Load data example NHANES dataset (more information in `baysc` package)
load(file = "Model_Code/data_nhanes.rda")

# Categorical exposure matrix of food groups, nxJ
x_mat <- as.matrix(data_nhanes[, c(11:38)])

# Subpopulation indicators using age, nx1
h_all <- as.numeric(data_nhanes$age_cat)

# Binary outcome of hypertension
y_all <- data_nhanes$BP_flag

# Survey stratum indicators
stratum_id <- data_nhanes$stratum_id
# Survey cluster indicatos
cluster_id <- data_nhanes$cluster_id
# Survey sampling weights
sampling_wt <- data_nhanes$sample_wt

#==== Run WRPC model ====
res <- wrpc(x_mat = x_mat, h_all = h_all, sampling_wt = sampling_wt,
             cluster_id = cluster_id, stratum_id = stratum_id,
             run_sampler = "both", adapt_seed = 1, class_cutoff_global = 0.05,
             n_runs = 50, burn = 25, switch = 10, thin = 1, save_res = FALSE)

#==== Display results ====
# Global pattern profiles
plot_wrpc_global_pattern_profiles(res = res)
# Global pattern item consumption level probabilities
plot_wrpc_global_pattern_probs(res = res)
# Local pattern deviations for items allocated to local (nu < 0.5)
plot_wrpc_local_profiles_allocation(res = res)
# Full local pattern profiles
plot_wrpc_local_pattern_profiles(res = res)
# Global-local deviation
plot_wrpc_allocation(res = res)
# Distribution of global classes by subpopulation
plot_wrpc_class_subgroup_dist(res = res)


#==== Example for unweighted model without survey weights ====
# Leave stratum_id, cluster_id, and sampling_wt as NULL
res_unwt <- wrpc(x_mat = x_mat, h_all = h_all, run_sampler = "both",
                 adapt_seed = 1, class_cutoff_global = 0.05, n_runs = 50,
                 burn = 25, switch = 10, thin = 1, save_res = FALSE)
```

<div id='id-section2'/>

## Contributing and Getting Help

Please report bugs by opening an [issue](https://github.com/smwu/WRPC/issues/new/choose). If you wish to contribute, please make a pull request. If you have questions, you can open a [discussion thread](https://github.com/smwu/WRPC/discussions).
