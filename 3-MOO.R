library(rmoo)




# prep data
y_prep <- as.numeric(if_else(df$liver_disease == "LD",1,0))
X_prep <- model.matrix(~ . - liver_disease, df)[, -1]

# scale predictors
X_scaled <- scale(X_prep)
x_center <- attr(X_scaled, "scaled:center")
x_scale <- attr(X_scaled, "scaled:scale")

# unscale sex variable
X_scaled[,"admin_sexMale"] <- if_else(X_scaled[,"admin_sexMale"] > 0,1,0)

df_prep <- data.frame(X_scaled, liver_disease = y_prep)
  
########################

logit_fitness <- function(x) {
  # If x is a vector → turn into a 1-row matrix
  if (is.null(dim(x))) {
    x <- matrix(x, nrow =1)
  }
  
  # x = coefficient vector (beta0, beta1, beta2, ...)
  
  ## prep data
  X <- as.matrix(df_prep[, setdiff(names(df_prep), "liver_disease")],,drop=FALSE)
  y <- df_prep$liver_disease
  
  # add intercept column
  X_aug <- cbind(1, X)  
  
  ## matrix multiply
  z <- X_aug %*% t(x)
  
  ## logistic transform (plus overflow safety)
  p_hat <- 1 / (1 + exp(-pmin(pmax(z, -20), 20)))
  
  ## RMSE calc
  rmse <- sqrt(mean((p_hat - y)^2))
  ## FNR calc
  fnr <- sum((p_hat < 0.5) & y == 1) / sum(y == 1)
  
  # Return a 2-column matrix: (RMSE, FNR)
  return(c(rmse, fnr))
}


########################

p <- ncol(df_prep) - 1              # number of predictors

# bounds for coefs
lower_bounds <- rep(0, p+1)  
upper_bounds <- rep(1,  p+1)

## run nsga2
result <- rmoo::rmoo(
  fitness  = logit_fitness,
  type  = "real-valued",
  algorithm = "NSGA-II",
  lower = lower_bounds,
  upper = upper_bounds,
  nObj = 2,
  popSize = 100,
  pcrossover = 0.9,
  maxiter  = 10000,
  names = c("RMSE", "FNR"),
  monitor = FALSE,
  seed = 42
)

