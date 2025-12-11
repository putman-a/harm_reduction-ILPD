
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
  


# Logistic fitness function for MOO -----------------------------

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
  popSize = 300,  # selected from testing
  pcrossover = 0.9,
  maxiter  = 10000,  # selected from testing
  names = c('LL','FNR'),
  monitor = FALSE,
  seed = 42
)
result

# # save result
# write_rds(x = result, file = "saved_MOO_ouput/pop300max10000.RDS")

# load results
result <- read_rds(file = "saved_MOO_ouput/pop300max10000.RDS")

# put results in a data frame
result_df <- MOO_results(result)

# plot MOO suggestions
MOO_plot(result_df,pop = 300,miter = 10000) +
  labs(title = "NSGA-II")



# plot for paper ----------------------------------------------------------

# define polygon coordinates for the CI shape
x_points <- c(-317.7, -290.3, -271.7,-290.3)
y_points <- c(.207, .124, .207, .290)

# rename column
plot_df <- result_df %>% mutate("Female FNR" = `FNR for Females`)

# plot
ggplot(plot_df, aes(x = `Negative Log Likelihood`, y = `Female FNR`)) +
  ## scatterplot
  geom_point() +
  ## draw polygon representing the 95%CI
  annotate(
    geom = "polygon", 
    x=x_points,
    y=y_points, 
    fill="purple",
    alpha=.15
    ) +
  ## add the SO point
  annotate(
    geom = "point", 
    x=-290.3, 
    y=0.207,
    colour="purple", 
    size=5, 
    shape=18
    ) +
  ## label for SO point
  annotate(
    geom = "text", 
    x=-303, 
    y=.29, 
    label="SO (NLL)",
    size=3.5
    ) +
  ## arrow from label to SO point
  annotate(
    geom = "curve", 
    x = -299.5, 
    y = .289, 
    xend =-290.4, 
    yend = 0.215,
    curvature = -.3, 
    arrow = arrow(length = unit(2, "mm"))
  ) +
  ## label for 95%CI area
  annotate(
    geom = "text", 
    x=-312, 
    y=0.25, 
    label="SO 95%CI",
    size=3.5
    ) +
  ## arrow from label to shaded area
  annotate(
    geom = "curve", 
    x = -312, 
    y = .242, 
    xend =-307, 
    yend = 0.20,
    curvature = .3, 
    arrow = arrow(length = unit(2, "mm")),
  ) +
  ## label for NSGA-II points
  annotate(
    geom = "text", 
    x=-307, 
    y=0.02, 
    label="NSGA-II",
    size=3.5
    ) +
  ## arrow from label to NSGA-II points
  annotate(
    geom = "curve", 
    x = -304, 
    y = .02, 
    xend =-292, 
    yend = 0.05,
    curvature = .2, 
    arrow = arrow(length = unit(2, "mm"))
  ) +
  ## subtitle
  labs(subtitle = "NSGA-II: Population = 300, Max. Iter. = 10000") +
  ## y and x limits
  ylim(0,.3) +
  scale_x_reverse(limits = c(-320,-270)) +
  ## classic theme
  theme_classic()

# save plot as png
ggsave(, filename = "saved_MOO_ouput/MOOplotwCI.png", device = "png",scale = 3,
       width = 600, height = 600,units = "px")


# NSGA-II population summary ----------------------------------------------


## All results -------------------------------------------------------------

# create tibble of outputs
NSGA_result_summary <- tibble(
  Measure = c(
    "Mean(SD)",
    "Median(IQR)",
    "Min:Max"
  ),
  NLL = c(
    ## mean(sd)
    paste0(
      round(mean(result_df$`Negative Log Likelihood`),1),
      "(",
      round(sd(result_df$`Negative Log Likelihood`),2),
      ")"
    ),
    ## median(IQR)
    paste0(
      round(median(result_df$`Negative Log Likelihood`),1),
      "(",
      round(IQR(result_df$`Negative Log Likelihood`),2),
      ")"
    ),
    ## min:max
    paste0(
      round(min(result_df$`Negative Log Likelihood`),1),
      ":",
      round(max(result_df$`Negative Log Likelihood`),1)
    )
  ),
  "FNR in Females" = c(
    ## mean(sd)
    paste0(
      round(mean(result_df$`FNR for Females`),4),
      "(",
      round(sd(result_df$`FNR for Females`),5),
      ")"
    ),
    ## median(IQR)
    paste0(
      round(median(result_df$`FNR for Females`),4),
      "(",
      round(IQR(result_df$`FNR for Females`),5),
      ")"
    ),
    ## min:max
    paste0(
      round(min(result_df$`FNR for Females`),4),
      ":",
      round(max(result_df$`FNR for Females`),4)
    )
  )
)


## Weak accept -------------------------------------------------------------

# select weak accept
weak <- result_df %>%
  filter(
    `FNR for Females` < 0.124 | `Negative Log Likelihood` > -271.7
  )

## add to results
# create tibble of outputs
NSGA_result_plus_weak <- NSGA_result_summary %>% 
  mutate(
  Measure_weak = c(
    "Mean(SD)",
    "Median(IQR)",
    "Min:Max"
  ),
  NLL_weak = c(
    ## mean(sd)
    paste0(
      round(mean(weak$`Negative Log Likelihood`),1),
      "(",
      round(sd(weak$`Negative Log Likelihood`),2),
      ")"
    ),
    ## median(IQR)
    paste0(
      round(median(weak$`Negative Log Likelihood`),1),
      "(",
      round(IQR(weak$`Negative Log Likelihood`),2),
      ")"
    ),
    ## min:max
    paste0(
      round(min(weak$`Negative Log Likelihood`),1),
      ":",
      round(max(weak$`Negative Log Likelihood`),1)
    )
  ),
  "FNR in Females_weak" = c(
    ## mean(sd)
    paste0(
      round(mean(weak$`FNR for Females`),4),
      "(",
      round(sd(weak$`FNR for Females`),5),
      ")"
    ),
    ## median(IQR)
    paste0(
      round(median(weak$`FNR for Females`),4),
      "(",
      round(IQR(weak$`FNR for Females`),5),
      ")"
    ),
    ## min:max
    paste0(
      round(min(weak$`FNR for Females`),4),
      ":",
      round(max(weak$`FNR for Females`),4)
    )
  )
)


# Moderate accept ---------------------------------------------------------

# select moderate accept
moderate <- result_df %>%
  filter(
    (`FNR for Females` < 0.124 & `Negative Log Likelihood` > -271.7) | 
      (`FNR for Females` < (0.124 - 0.083) ) |
      (`Negative Log Likelihood` > (-271.7 + 18.6))
  )

## add to results
# create tibble of outputs
NSGA_result_plus_mod <- NSGA_result_plus_weak %>% 
  mutate(
    Measure_mod = c(
      "Mean(SD)",
      "Median(IQR)",
      "Min:Max"
    ),
    NLL_mod = c(
      ## mean(sd)
      paste0(
        round(mean(moderate$`Negative Log Likelihood`),1),
        "(",
        round(sd(moderate$`Negative Log Likelihood`),2),
        ")"
      ),
      ## median(IQR)
      paste0(
        round(median(moderate$`Negative Log Likelihood`),1),
        "(",
        round(IQR(moderate$`Negative Log Likelihood`),2),
        ")"
      ),
      ## min:max
      paste0(
        round(min(moderate$`Negative Log Likelihood`),1),
        ":",
        round(max(moderate$`Negative Log Likelihood`),1)
      )
    ),
    "FNR in Females_mod" = c(
      ## mean(sd)
      paste0(
        round(mean(moderate$`FNR for Females`),4),
        "(",
        round(sd(moderate$`FNR for Females`),5),
        ")"
      ),
      ## median(IQR)
      paste0(
        round(median(moderate$`FNR for Females`),4),
        "(",
        round(IQR(moderate$`FNR for Females`),5),
        ")"
      ),
      ## min:max
      paste0(
        round(min(moderate$`FNR for Females`),4),
        ":",
        round(max(moderate$`FNR for Females`),4)
      )
    )
  )


## Strong accept -----------------------------------------------------------

# select strong accept
strong <- result_df %>%
  filter(
      (`FNR for Females` < (0.124 - 0.083) ) &
      (`Negative Log Likelihood` > (-271.7 + 18.6))
  )

# no models meet this criteria!!!


# Save results ------------------------------------------------------------

write_csv(NSGA_result_plus_mod, file = "data/NSGA2_summary_table.csv")
