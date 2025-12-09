
# Set up ------------------------------------------------------------------
# load libraries
library(tidyverse)
# library(boot)
# library(rmoo)

# set seed
set.seed(42)

# import data
df <- read_rds(file = "data/cleaned_ILPD.RDS")


# Prep data for use with `rmoo` -------------------------------------------
# separate outcome and predictors
y_prep <- as.numeric(if_else(df$liver_disease == "LD",1,0))
X_prep <- model.matrix(~ . - liver_disease, df)[, -1]

# scale predictors
X_scaled <- scale(X_prep)

# # use to extract values used to scale and centre, respectively
# x_center <- attr(X_scaled, "scaled:center")
# x_scale <- attr(X_scaled, "scaled:scale")

# un-scale admin sex variable
X_scaled[,"admin_sexMale"] <- if_else(X_scaled[,"admin_sexMale"] > 0,1,0)

df_prep <- data.frame(X_scaled, liver_disease = y_prep)
  


# Set up logistic regression function for MOO -----------------------------

# +++++++++++++++++++++++++++++++++++++
# Logistic regression fitness function 
# +++++++++++++++++++++++++++++++++++++
# x : a matrix input from `rmoo()`
logit_fitness <- function(x) {
  ## if x is a vector → turn into a 1-row matrix
  if (is.null(dim(x))) {
    x <- matrix(x, nrow =1)
  }
  
  ## prep data
  X <- as.matrix(df_prep[, setdiff(names(df_prep), "liver_disease")],drop=FALSE)
  y <- df_prep$liver_disease
  admin_sex <- df_prep$admin_sexMale
  
  ## add intercept column
  X_int <- cbind(1, X)  
  
  ## matrix multiply
  z <- X_int %*% t(x)
  
  ## logistic regression prediction 
  ## with overflow safety restricting possible coefs to p/m 10
  p_hat <- 1 / (1 + exp(-pmin(pmax(z, -10), 10)))
  
  ## Log-Likelihood calc
  LL <- sum(y * log(p_hat) + (1 - y) * log(1 - p_hat))
  
  ## female FNR calc
  FN_fem <- sum((p_hat < 0.5) & admin_sex == 0 & y == 1) # females pred no LD when LD is present
  pos_fem <- sum(admin_sex == 0 & y == 1) # females with LD
  FNR <- FN_fem / pos_fem # FNR calc
  
  # Return a 2-column matrix: (LL, FNR)
  return(c(-LL, FNR)) ## -LL because NGSA-II is a minimization function
}


# NSGA-II using `rmoo()` ----------------------------------------------------------

# number of predictors
p <- ncol(df_prep) - 1 # 7

# bounds for coefs
lower_bounds <- rep(-2, p+1)  
upper_bounds <- rep(2,  p+1)

## run nsga2
result <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  popSize = 300,
  pcrossover = 0.9,
  maxiter  = 10000,
  names = c('LL','FNR'),
  monitor = FALSE,
  seed = 42
)
result

# save result
write_rds(x = result, file = "saved_MOO_ouput/pop300max10000.RDS")

# put results in a data frame
result_df <- tibble(
  "Negative Log Likelihood" = -result@fitness[,1],
  "FNR for Females" = result@fitness[,2]
)

# plot
ggplot(result_df, aes(x = `Negative Log Likelihood`, y = `FNR for Females`)) +
  geom_point() +
  annotate("point", x=-290.3,y=0.207, colour = "blue",size=4, shape=17) +
  annotate("text",x=-290.4,y=.222, label="Standard MLE",colour="blue") +
  labs(
    title = "NSGA-II-derived Logistic Regression Performance",
    caption = "Population = 300, Max Iter. = 10,000"
    ) +
  ylim(0,.25) +
  scale_x_reverse(limits = c(-294,-290)) +
  theme_classic()
