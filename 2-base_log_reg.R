
# Custom functions --------------------------------------------------------

# +++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate accuracy of logistic regression prediction
# +++++++++++++++++++++++++++++++++++++++++++++++++++++
# tab : a confusion table of predicted and actual values
acc <- function(tab){
  TP <- tab[1,1]
  FP <- tab[1,2]
  FN <- tab[2,1]
  TN <- tab[2,2]
  acc <- (TP + TN) / (TP + FP + TN + FN)
  return(acc)
} 


# +++++++++++++++++++++++++++++++++++++
# accuracy function for boot strapping
# +++++++++++++++++++++++++++++++++++++
# data : a data frame with the predicted outcome and actual outcomes
# index : an integer value used by the `boot` function to keep track of iterations
acc_boot <- function(data, index) {
  ## update index for boot function
  data <- data[index, ]
  ## accuracy
  results_tab <- table(data$pred_liver_disease, data$liver_disease)
  results_acc <- acc(results_tab)
  ## return accuracy
  return(results_acc)
}


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


# ++++++++++++++++++++++++++++++++++++++++
# female FNR  function for boot strapping
# ++++++++++++++++++++++++++++++++++++++++
# data : a data frame with the predicted outcome, actual outcome, and admin sex
# index : an integer value used by the `boot` function to keep track of iterations
FNR_F_boot <- function(data, index) {
  ## update index for boot function
  data <- data[index, ]
  ## split by sex
  F_data <- data %>% filter(., admin_sex == "Female")
  ## female fnr
  F_tab <- table(F_data$pred_liver_disease, F_data$liver_disease)
  F_fnr <- FNR(F_tab)
  ## return fnr
  return(F_fnr)
}


# Set up ------------------------------------------------------------------
# load libraries
library(tidyverse)
library(boot)
# library(rmoo)

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
  if_else(log_reg_pred_probs > 0.5, "pred nonLD", "pred LD"),
  levels = c("pred LD", "pred nonLD")
  )

## create pred vs. actual df, including sex for stratification
log_reg_results <- df %>% mutate(pred_liver_disease = log_reg_pred_class) %>%
  select(., admin_sex, liver_disease, pred_liver_disease)

## export data frame with predicted results as RDS
write_rds(log_reg_results,file = "data/base_logreg_results.RDS")


### Base accuracy -----------------------------------------------------------
log_reg_results_tab <- table(log_reg_results$pred_liver_disease, log_reg_results$liver_disease)
log_reg_results_acc <- acc(log_reg_results_tab)

### Bootstrap 95% CI accuracy ------------------------------------------------

boot_acc <- boot(data = log_reg_results, statistic = acc_boot, R = 10000)
boot_acc
boot_ci_acc <- boot.ci(boot.out = boot_acc, type = "norm")
boot_ci_acc 


### Base FNR  ---------------------------------------------------------------

# overall FNR
log_reg_results_FNR <- FNR(log_reg_results_tab)

# calculate FNR for female
log_reg_results_F <- log_reg_results %>% filter(admin_sex == "Female")
log_reg_results_F_tab <- table(log_reg_results_F$pred_liver_disease, log_reg_results_F$liver_disease)
log_reg_results_F_FNR <- FNR(log_reg_results_F_tab)


### Bootstrap 95% CI FNR for females ----------------------------------------

boot_FNR <- boot(data = log_reg_results, statistic = FNR_F_boot, R = 10000)
boot_FNR
boot_ci_FNR <- boot.ci(boot.out = boot_FNR, type = "norm")
boot_ci_FNR 

# Results from base log reg -----------------------------------------------
base_logreg_results <- tibble(
  "Method" = c("Base Logisitic Regression"),
  "Optimization Objectives" = c("Negative Log-Likelihood"),
  "Log-Likelihood" = c(as.numeric(logLik(log_reg))),
  "Accuracy" = c(round(log_reg_results_acc,3)),
  "Accuracy Bootstrapped 95% CI" = c(
    paste0(
      "(", round(boot_ci_acc$normal[,2],3), ",", round(boot_ci_acc$normal[,3],3),")",
      collapse = ""
    )
  ),
  "FNR for Females" = c(round(log_reg_results_F_FNR,3)),
  "FNR Bootstrapped 95% CI" = c(
    paste0(
      "(", round(boot_ci_FNR$normal[,2],3), ",", round(boot_ci_FNR$normal[,3],3),")",
      collapse = ""
    )
  )
)
base_logreg_results
# save results
write_csv(base_logreg_results, file = "data/base_log_results_table.csv")
