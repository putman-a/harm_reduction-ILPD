
# Custom functions --------------------------------------------------------

# ++++++++++++++++++++++++++++
# Hmsic rcorr() matrix to table
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
cormat_to_table <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  tibble(
    "row variable" = rownames(cormat)[row(cormat)[ut]],
    "column variable" = rownames(cormat)[col(cormat)[ut]],
    "r"  = (cormat)[ut],
    "|r|" = abs((cormat)[ut]),
    "p value" = pmat[ut]
  )
}


# Set up -----------------------------------------------------------------

# load libraries
library(tidyverse)
library(mice)
# library(Hmisc) ## required but is called directly
library(rmoo)

# set seed
set.seed(42)

# import data
ilpd <- read_csv(file = "data/ILPD.csv")
ilpd_DD <- read_csv(file = "data/ILPD-DD.csv")


# Clean data --------------------------------------------------------------

# Tidy variable types for analysis
df <- ilpd %>% 
  mutate(
    admin_sex = as_factor(admin_sex),
    liver_disease = factor(
      if_else(
        liver_disease == 1,
        "LD",
        "nonLD"
        )
      )
    )


# MICE imputation ---------------------------------------------------------

# check for missing values
summary(df) # alkphos has 4 missing values

# identify rows with missing data
which(is.na(df$alkphos)) # rows 210, 242, 254, 313

# perform 5 multiple imputations using 50 iterations
imputed_df <- mice(df, m = 5, maxit = 50, method = 'pmm', seed = 42)

# identify the imputed values
imputed_values <-imputed_df$imp$alkphos

# average the imputed values per observation
averaged_imputted_values <- tibble(rowMeans(imputed_values, na.rm = TRUE))

# new df with missing values replace with average of the 5 imputed values
df_imp <- df
df_imp[210,"alkphos"] <- averaged_imputted_values[1,]
df_imp[242,"alkphos"] <- averaged_imputted_values[2,]
df_imp[254,"alkphos"] <- averaged_imputted_values[3,]
df_imp[313,"alkphos"] <- averaged_imputted_values[4,]


# Test for multicollinearity of predictor variables --------------------------

# test variables for highly correlated predictor variables (multicollinearity mitigation)

## create a df to work with cor() and Hmisc's rcorr()
df_cortest <- df_imp %>% mutate(
  admin_sex = if_else(admin_sex == "Female",1,0)
  ) %>% select(-c("liver_disease"))

## check correlations using cor()
cor_res <- cor(df_cortest)
round(cor_res,4)

## use Hmsic rcorr() to add p values for the r values
cor_res_pval <- Hmisc::rcorr(as.matrix(df_cortest))

## format corr() results into table
flattened_cor_res <- cormat_to_table(round(cor_res_pval$r,4), cor_res_pval$P)

## inter-predictor with absolute r values > 0.5
filtered_cor_res <- flattened_cor_res %>% filter(`|r|` >0.5)
filtered_cor_res

## Investigate theoretical basis for inclusion/exclusion ---------------------

### total vs. direct bilirubin
# if total bilirubin is elevated but most of the elevation is unconjugated 
# (indirect) bilirubin than hepatocellular disease is unlikely to be the cause
# therefore, we select direct bilirubin for the model

### albumin and a/g ratio
# a/g ratio is useful in differential diagnosis of liver and kidney disorders,
# however, is of limited clinical use without also knowing what the albumin 
# and globulin levels are separately
# therefore, we select albumin for the model

### SGPT and SGOT
# AKA AST and ALT, ALT(SGOT) is a more specific marker of hepatocellular injury
# therefore, we select SGOT for the model

### SGOT and alkphos
# AKA ALT and ALP, ALP(alkphos) is primarily indicative of cholestasis 
# (blocked bile ducts) whereas ALT(SGOT) is primarily indicative of 
# hepatocellular necrosis (liver cell death). While these can co-occur, the
# markers arise from distinct physiological mechanisms.
# therefore, we will use both variables in the model


# Create final dataset and export -----------------------------------------

# remove co-linear predictors
df_model <- df_imp %>%
  ## remove total bilirubin, a/g ratio, SGPT
  select(-c("tot_billirubin", "ag_ratio", "sgpt"))

# export df to RDS
write_rds(df_model,file = "data/cleaned_ILPD.RDS")
