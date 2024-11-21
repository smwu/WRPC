### Prep FPQ 50 item data

rm(list = ls())

library(tidyverse)  # data wrangling
library(baysc)      # bayesian survey clustering
library(haven)      # read sas file
library(naniar)     # missingness
library(sfsmisc)    # mult.fig
library(readxl)     # read excel file
library(data.table)

wd <- "~/Documents/Github/WRPC/"
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WRPC/"
code_dir <- "Model_Code/"
data_dir <- "Application/HCHS_Data/"
res_dir <- "Results/"

# Read in binned FPQ data
bin_fpq_vars <- read_sas(paste0(wd, data_dir, "fpq_bins.sas7bdat"))
# Read in raw FPE data: n=13866
raw_fpq_vars <- read_sas(paste0(wd, data_dir, "fpe_inv4.sas7bdat"))
# Raw variable labels
vars_labels = data.frame(label = sapply(raw_fpq_vars, function(x){attr(x,"label")}))
# Read in cleaned FPQ data
sol_data <- read.csv(paste0(wd, data_dir, "FPQ129_CVD07FEB2023.csv"))
# Select diet variables
fpq_data <- sol_data %>% 
  select(ID, STRAT, PSU_ID, CENTER, BKGRD1_C6, fp2:fp1)

# # Read in 50-item FPQ data: n = 16415
# fpq_52 <- read.csv(paste0(wd, data_dir, "FPQ52_CVD27MAY2020.csv"))
# # Read in labels
# fpq_52_labs <- readxl::read_xlsx(paste0(wd, data_dir, 
#                                         "NDSR_food_groups_with_Levels.xlsx"))
# # Get levels
# fpq_52_labs <- fpq_52_labs %>%
#   mutate(level_labs = case_when(
#     Levels == 4
#   ))
# 
# # Drop those with NA: n = 12772
# fpq_52_dropna <- fpq_52 %>%
#   drop_na(fp1:fp52)
# 
# # Get proportion in each food level
# apply(fpq_52_dropna %>% select(fp1:fp52), 2, 
#       function(x) round(prop.table(table(x, useNA = "always")), 3))
# 
# # Find 5-level foods
# which(apply(fpq_PR[, -c(1:5)], 2, 
#             function(x) max(x, na.rm = TRUE) == 5))
# # Convert 5-level foods to 4 levels
# study_data <- study_data %>%
#   mutate_at(c("fp2", "fp116", "fp7", "fp10", "fp11", "fp12", "fp13", "fp17", 
#               "fp53", "fp101"), ~ifelse(.x == 5, 4, .x))
# # Check no 5-level foods
# apply(study_data %>% select(fp2:fp1), 2, function(x) max(x, na.rm = TRUE))

#===========
# Binary foods: 115A, 115B, 115C, 115D
which(apply(raw_fpq_vars[, -c(1:3)], 2, 
            function(x) length(unique(x[!is.na(x)])) == 2))

# Frequencies equivalent in R
table_data <- apply(raw_fpq_vars %>% select(FPE115A:FPE115D), 2, 
      function(x) round(prop.table(table(x, useNA = "always")), 3))
print(table_data)

# Drop NAs: n = 12761
raw_fpq_vars_dropna <- raw_fpq_vars %>% drop_na(-FPE0A)
# Check NA harmonization
# Previous cleaned dataset
# Detected: some individuals with FPQ data but with missing Date of Completion.
# => Include these individuals in the dataset.
fpq_data <- fpq_data %>% drop_na()
missing <- setdiff(fpq_data$ID, as.numeric(raw_fpq_vars_dropna$ID))
missing_df <- fpq_data %>% filter(ID %in% missing)
raw_missing_df <- raw_fpq_vars %>% filter(ID %in% missing)

### Create Derived Variables
fpq_vars_derv <- raw_fpq_vars_dropna %>%
  mutate(
    # Fruit drinks
    FPE5F_reg = FPE5F * (1 - FPE5AP), # regular fruit drink
    FPE5F_diet = FPE5F * FPE5AP, # diet fruit drink
    # Soda
    FPE8F_reg = FPE8F * (1 - FPE8AP), # regular soda
    FPE8F_diet = FPE8F * (1 - FPE8AP), # diet soda
    # Cooked cereal
    FPE11F_oatmeal = FPE11F * FPE11AP, # oatmeal
    FPE11F_oth = FPE11F * (1 - FPE11AP), # other
    # Cold cereal
    FPE12F_whole = FPE12F * FPE12AP, # whole grain type
    FPE12F_not = FPE12F * (1 - FPE12AP), # not whole grain
    # Lettuce salad
    FPE38F_green = FPE38F * FPE38AP, # dark green leaf salad
    FPE38F_oth = FPE38F * (1 - FPE38AP), # other lettuce salad
    # Tortillas toacos
    FPE46F_corn = FPE46F * FPE46AP, # corn tortillas/tacos
    FPE46F_oth = FPE46F * (1 - FPE46AP), # other tortillas/tacos
    # Rice
    FPE53F_whole = FPE53F * FPE53AP, # whole grain rice
    FPE53F_oth = FPE53F * (1 - FPE53AP), # other rice
    # Sandwich bread
    FPE56F_white = FPE56F * FPE56AP, # white sandwich bread
    FPE56F_non = FPE56F * (1 - FPE56AP), # non-white sandwich bread
    # Bread/rolls
    FPE57F_white = FPE57F * FPE57AP, # white bread/rols
    FPE57F_non = FPE57F * (1 - FPE57AP), # non-white bread/rolls
    # Soup
    FPE80F_bean = FPE80F * FPE80AP, # bean soup
    FPE80F_cream = FPE80F * FPE80BP, # cream soup
    FPE80F_tomato = FPE80F * FPE80CP, # tomato/veg soup
    FPE80F_oth = pmax(0, FPE80F * (1 - FPE80AP - FPE80BP - FPE80CP)), # all other soup
    # Pizza
    FPE81F_pepp = FPE81F * FPE81AP, # pepperoni pizza
    FPE81F_non = FPE81F * (1 - FPE81AP), # non-pepperoni pizza
    # Pie
    FPE99F_fruit = FPE99F * FPE99AP, # fruit pie
    FPE99F_oth = FPE99F * (1 - FPE99AP) # other pies
  , .keep = "unused")  # drop used columns

### Create 49 food groups
fpq_vars_groups <- fpq_vars_derv %>%
  mutate(
    citrus_juice = FPE1F,
    fruit_juice_not_citrus = FPE2F + FPE3F + FPE4F,
    applesauce = FPE13F,
    citrus_fruit = FPE23F + FPE24F, 
    fruit_not_citrus = FPE14F + FPE15F + FPE16F + FPE17F + FPE18F + 
      FPE19F + FPE20F + FPE21F + FPE22F + FPE25F + FPE26F + FPE27F,
    avocado = FPE49F,
    fruit_fried = FPE51F,
    
    veg_dark_green = FPE28F + FPE29F + FPE34F + FPE38F_green,
    veg_yellow = FPE30F + FPE37F + FPE39F + FPE48F,
    tomato = FPE36F + FPE43F + FPE80F_tomato + FPE44F, 
    potato_white = FPE41F + FPE42F,
    potato_fried = FPE40F,
    veg_oth_starchy = FPE32F + FPE33F,
    legumes = FPE45F + FPE47F + FPE80F_bean,
    veg_oth = FPE31F + FPE35F + FPE38F_oth + FPE50F + FPE52F,
    
    grains_whole = FPE11F_oatmeal + FPE53F_whole + FPE12F_whole,
    grains_refined = FPE11F_oth + FPE53F_oth + FPE12F_not,
    bread_corn = FPE46F_corn + FPE83F,
    bread_refined = FPE46F_oth + FPE54F + FPE55F + FPE56F_white + FPE57F_white + 
      FPE84F,
    bread_oth = FPE56F_non + FPE57F_non,
    pizza = FPE81F_pepp + FPE81F_non,
    savory_snacks = FPE82F + FPE85F + FPE86F,
    sweet_baked_goods = FPE94F + FPE95F + FPE96F + FPE97F + FPE98F + 
      FPE99F_fruit + FPE99F_oth,
    
    beef = FPE64F + FPE65F + FPE67F + FPE68F + FPE69F,
    pork = FPE70F + FPE74F + FPE75F,
    poultry = FPE71F + FPE72F + FPE73F,
    cold_cuts = FPE59F + FPE60F + FPE61F + FPE62F,
    sausages_cured = FPE76F + FPE66F + FPE77F,
    fish = FPE63F + FPE78F + FPE79F,
    
    eggs = FPE102F,
    nuts = FPE87F,
    milk = FPE6F + FPE110F,
    cheese = FPE89F + FPE90F,
    yogurt = FPE88F,
    cream = FPE109F + FPE112F + FPE113F + FPE80F_cream,
    cream_nondairy = FPE108F,
    dairy_dessert = FPE91F + FPE92F + FPE93F,
    
    oil = FPE114F, # remove questions on oil types: FPE115A + FPE115B + FPE115C + FPE115D,
    sugar_sweetener = FPE58F + FPE106F + FPE111F,
    artif_sweetener = FPE107F,
    candy = FPE100F + FPE101F,
    
    bevs_sugar = FPE8F_reg + FPE5F_reg,
    bevs_diet = FPE8F_diet + FPE5F_diet,
    bevs_meal_repl = FPE7F,
    tea = FPE104F + FPE105F,
    coffee = FPE103F,
    beer = FPE9F,
    wine = FPE10F,
    soup_oth = FPE80F_oth,
  .keep = "unused") # drop used columns
# Drop the oil columns
fpq_vars_groups <- fpq_vars_groups %>%
  select(-c(FPE115A, FPE115B, FPE115C, FPE115D))

# Frequencies
table_data <- apply(fpq_vars_groups[, -c(1:3)], 2, 
                    function(x) round(prop.table(table(x, useNA = "always")), 3))
print(table_data)

### Consumption Level Binning
# Consumption levels: 1) Never, 2) 2-3 times per month or less, 
# 3) 1-2 times per week, 4) 3-6 times per week, 5) 1 time per day or more
totfoods <- ncol(fpq_vars_groups[, -c(1:3)])
fpq_categ <- fpq_vars_groups %>%
  mutate(across(citrus_juice:soup_oth, ~ case_when(
    . == 0 ~ 1, # No consumption
    . > 0 & . < (1 / 7) ~ 2,  # A few times per month (less than once per week)
    . >= (1 / 7) & . < (5 / 7) ~ 3, # A few (1-4) times per week
    . >= (5 / 7) & . <= 1 ~ 4, # Daily or almost daily
    . > 1 ~ 5, # More than daily
  )))

# Frequencies
table_data <- apply(fpq_categ[, -c(1:3)], 2, 
                    function(x) round(prop.table(table(x, useNA = "always")), 3))
print(table_data)

# ### Collapse upper bins with less than 5%. Do not collapse lower bins
# fpq_categ_bins <- fpq_categ %>%
#   mutate(
#     citrus_juice = ifelse(citrus_juice == 5, 4, citrus_juice),
#     applesauce = ifelse(applesauce %in% c(4, 5), 3, applesauce),
#     avocado = ifelse(avocado == 5, 4, avocado),
#     fruit_fried = ifelse(fruit_fried %in% c(4, 5), 3, fruit_fried),
#     potato_white = ifelse(potato_white %in% c(4, 5), 3, potato_white),
#     potato_fried = ifelse(potato_fried %in% c(4, 5), 3, potato_fried),
#     veg_oth_starchy = ifelse(veg_oth_starchy == 5, 4, veg_oth_starchy),
#     bread_oth = ifelse(bread_oth == 5, 4, bread_oth),
#     pizza = ifelse(pizza %in% c(4, 5), 3),
#     savory_
#     )


# Save final data
write.csv(fpq_categ, paste0(wd, data_dir, "fpq_49categ_2024Nov19.csv"), 
          row.names = FALSE)

