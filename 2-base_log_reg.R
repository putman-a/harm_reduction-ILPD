
# Custom functions --------------------------------------------------------

# ++++++++++++++++++++++++++++
# Calculate FNR
# ++++++++++++++++++++++++++++
# tab : a confusion table of predicted and actual values
FNR <- function(tab){
  TP <- tab[1,1]
  FP <- tab[1,2]
  FN <- tab[2,1]
  TN <- tab[2,2]
  FNR <- FN / (TP + FN)
  return(FNR)
} 


# +++++++++++++++++++++++++++++++++++++
# FNR diff function for boot strapping
# +++++++++++++++++++++++++++++++++++++
# data : a data frame with the predicted outcome, actual outcome, and admin sex
# index : an integer value used by the `boot` function to keep track of iterations
FNR_diff_sex_boot <- function(data, index) {
  # update index for boot
  data <- data[index, ]
  ## split by sex
  Fem <- data %>% filter(.,admin_sex == "Female")
  Male <- data %>% filter(.,admin_sex == "Male")
  ## female fnr
  F_tab <- table(Fem$pred_liver_disease, Fem$liver_disease)
  F_fnr <- FNR(F_tab)
  ## male fnr
  M_tab <- table(Male$pred_liver_disease, Male$liver_disease)
  M_fnr <- FNR(M_tab)
  ## fnr diff
  fnr_diff <- F_fnr - M_fnr
  return(fnr_diff)
}


# Set up ------------------------------------------------------------------
# load libraries
library(tidyverse)
library(boot)
library(rmoo)

# set seed
set.seed(42)

# import data
df <- read_rds(file = "data/cleaned_ILPD.RDS")


# Base logistic regression ------------------------------------------------
log_reg <- glm(
  formula = liver_disease ~ .,
  data = df,
  family = binomial(link = "logit")
)


## Calculate FNR by Admin. Sex ----------------------------------------------

# predict outcomes
log_reg_pred_probs <- predict(log_reg, type = "response")

## label predicted classifications
log_reg_pred_class <- factor(
  if_else(log_reg_pred_probs > 0.5, "nonLD", "LD"),
  levels = c("LD", "nonLD")
  )

## create pred vs. actual df, including sex for stratification
log_reg_results <- df %>% mutate(pred_liver_disease = log_reg_pred_class) %>%
  select(., admin_sex, liver_disease, pred_liver_disease)


### Base FNR  ---------------------------------------------------------------

# calculate for female
log_reg_results_F <- log_reg_results %>% filter(admin_sex == "Female")
log_reg_results_F_tab <- table(log_reg_results_F$pred_liver_disease, log_reg_results_F$liver_disease)
log_reg_results_F_FNR <- FNR(log_reg_results_F_tab)

# calculate for male
log_reg_results_M <- log_reg_results %>% filter(admin_sex == "Male")
log_reg_results_M_tab <- table(log_reg_results_M$pred_liver_disease, log_reg_results_M$liver_disease)
log_reg_results_M_FNR <- FNR(log_reg_results_M_tab)

# calculate diff
log_reg_FNR_diff <- round(log_reg_results_F_FNR - log_reg_results_M_FNR, 4)


### Bootstrap FNR  ---------------------------------------------------------
boot_FNR <- boot(data = log_reg_results, statistic = FNR_diff_sex_boot, R = 10000)
boot_FNR
boot_ci_FNR <- boot.ci(boot.out = boot_FNR, type = "norm")
boot_ci_FNR 


