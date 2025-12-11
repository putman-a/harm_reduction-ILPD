
# Set up ------------------------------------------------------------------

# load libraries
library(tidyverse)
# library(mice)
# library(boot)
# library(rmoo)

# set seed
set.seed(42)

# import data
df <- read_rds(file = "data/cleaned_ILPD.RDS")



## Utility functions ------------------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot trade-off between NLL and FNR in females from NSGA-II results
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# df : a tibble extracted from the NSGA-II object
# xlow : lower bound of the x-axis
# xhigh : upper bound of the x axis
MOO_plot <- function(df,pop,miter,xlow=-300,xhigh=-290){
  plot <- ggplot(df, aes(x = `Negative Log Likelihood`, y = `FNR for Females`)) +
    geom_point() +
    annotate(geom = "point", x=-290.3, y=0.207, fill="red", size=5, shape=23) +
    annotate(geom = "text", x=-290.3, y=0.207, label="Baseline", 
             hjust = -0.25,colour="red") +
    labs(subtitle = paste0("Population = ", pop, " Max Iter. = ", miter)) +
    ylim(0,.25) +
    scale_x_reverse(limits = c(xlow,xhigh)) +
    theme_classic()
  return(plot)
}


# ++++++++++++++++++++++++++++++++++++++++++++++
# Extract results from nsga2 object into tibble
# ++++++++++++++++++++++++++++++++++++++++++++++
# nsga2 : an nsga2 class object
MOO_results <- function(nsga2){
  tib <- tibble(
    "Negative Log Likelihood" = -nsga2@fitness[,1],
    "FNR for Females" = nsga2@fitness[,2]
  )
  return(tib)
}


## Prep data for use with `rmoo` ------------------------------------------
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


# number of predictors
p <- ncol(df_prep) - 1 # 7

# bounds for coefs
lower_bounds <- rep(-2, p+1)  
upper_bounds <- rep(2,  p+1)


# Testing parameters ------------------------------------------------------

# All of the following test use the following parameters:
#
# Fitness function = logistic regression coefficient and intercept creation
# Optimization target 1 = negative log likelihood (inverted as log likelihood 
#   so minimization decreases the value)
# Optimization target 2 = FNR for females (minimize number of females who have 
#   liver disease being incorrectly predicted as not having liver disease)
# Algorithm = "NSGA-II"
# Lower bound for coefficients = -2
# Upper bound for coefficients = 2
# Crossover probability = 0.9
# Crossover = simulated binary crossover
# Mutation = polynomial mutation
# 


# Max iter 100 ------------------------------------------------------------

## Pop 40 -----------------------------------------------------------------
max100pop40 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 40,
  maxiter  = 100,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max100pop40, file = "saved_MOO_ouput/pop40max100.RDS")

# extract results
max100pop40_results <- MOO_results(max100pop40)

# plot results
MOO_plot(max100pop40_results,40,100)


## Pop 100 ----------------------------------------------------------------
max100pop100 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 100,
  maxiter  = 100,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max100pop100, file = "saved_MOO_ouput/pop100max100.RDS")

# extract results
max100pop100_results <- MOO_results(max100pop100)

# plot results
MOO_plot(max100pop100_results,100,100)


## Pop 200 ----------------------------------------------------------------
max100pop200 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 200,
  maxiter  = 100,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max100pop200, file = "saved_MOO_ouput/pop200max100.RDS")

# extract results
max100pop200_results <- MOO_results(max100pop200)

# plot results
MOO_plot(max100pop200_results,200,100)


## Pop 300 ----------------------------------------------------------------
max100pop300 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 300,
  maxiter  = 100,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max100pop300, file = "saved_MOO_ouput/pop300max100.RDS")

# extract results
max100pop300_results <- MOO_results(max100pop300)

# plot results
MOO_plot(max100pop300_results,300,100)


# Max iter 1000 -----------------------------------------------------------

## Pop 40 -----------------------------------------------------------------
max1000pop40 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 40,
  maxiter  = 1000,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max1000pop40, file = "saved_MOO_ouput/pop40max1000.RDS")

# extract results
max1000pop40_results <- MOO_results(max1000pop40)

# plot results
MOO_plot(max1000pop40_results,40,1000)


## Pop 100 ----------------------------------------------------------------
max1000pop100 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 100,
  maxiter  = 1000,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max1000pop100, file = "saved_MOO_ouput/pop100max1000.RDS")

# extract results
max1000pop100_results <- MOO_results(max1000pop100)

# plot results
MOO_plot(max1000pop100_results,100,1000)


## Pop 200 ----------------------------------------------------------------
max1000pop200 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 200,
  maxiter  = 1000,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max1000pop200, file = "saved_MOO_ouput/pop200max1000.RDS")

# extract results
max1000pop200_results <- MOO_results(max1000pop200)

# plot results
MOO_plot(max1000pop200_results,200,1000)


## Pop 300 ----------------------------------------------------------------
max1000pop300 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 300,
  maxiter  = 1000,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max1000pop300, file = "saved_MOO_ouput/pop300max1000.RDS")

# extract results
max1000pop300_results <- MOO_results(max1000pop300)

# plot results
MOO_plot(max1000pop300_results,300,1000)


# Max iter 10000 ----------------------------------------------------------

## Pop 40 -----------------------------------------------------------------
max10000pop40 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 40,
  maxiter  = 10000,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max10000pop40, file = "saved_MOO_ouput/pop40max10000.RDS")

# extract results
max10000pop40_results <- MOO_results(max10000pop40)

# plot results
MOO_plot(max10000pop40_results,40,10000)


## Pop 100 ----------------------------------------------------------------
max10000pop100 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 100,
  maxiter  = 10000,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max10000pop100, file = "saved_MOO_ouput/pop100max10000.RDS")

# extract results
max10000pop100_results <- MOO_results(max10000pop100)

# plot results
MOO_plot(max10000pop100_results,100,10000)


## Pop 200 ----------------------------------------------------------------
max10000pop200 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 200,
  maxiter  = 10000,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max10000pop200, file = "saved_MOO_ouput/pop200max10000.RDS")

# extract results
max10000pop200_results <- MOO_results(max10000pop200)

# plot results
MOO_plot(max10000pop200_results,200,10000)


## Pop 300 ----------------------------------------------------------------
max10000pop300 <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  pcrossover = 0.9,
  popSize = 300,
  maxiter  = 10000,
  monitor = FALSE,
  seed = 42
)

# save result
write_rds(x = max10000pop300, file = "saved_MOO_ouput/pop300max10000.RDS")

# extract results
max10000pop300_results <- MOO_results(max10000pop300)

# plot results
MOO_plot(max10000pop300_results,300,10000)

