
# set up -----------------------------------------------------------------

# load libraries
library(tidyverse)
library(rmoo)

# set seed
set.seed(42)

# import data
ilpd <- read_csv(file = "data/ILPD.csv")
ilpd_DD <- read_csv(file = "data/ILPD-DD.csv")


# clean data --------------------------------------------------------------

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


# test variables for colinearity
df_cortest <- df %>% select(-c("admin_sex","liver_disease"))
cor(df_cortest)
                            