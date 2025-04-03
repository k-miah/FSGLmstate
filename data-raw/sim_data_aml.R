## code to prepare `sim_data_aml` dataset goes here

require(dplyr)


#################################
# Data generation: sim_data_aml #
#################################

# Data generation based on cause-specific hazards [Beyersmann et al., 2009;
# Beyersmann et al., 2012, Section 8.2 (pp. 172ff)]:


## generate_dat - R-function to simulate multi-state data (in wide format) of a
##                9-state model with q=8 transitions for AML based on
##                cause-specific hazards [Beyersmann et al., 2009]
##
## Input: n           [numeric]: Number of subjects
##        p           [numeric]: Number of covariates
##        X.design  [character]: Character string indicating whether covariates
##                               are drawn from a normal ('norm'; mu=0, sd=1) or
##                               binomial ('bin'; p=0.5) distribution
##        pars           [list]: List with variables
##                               - beta_12,...,beta_69 [transition-specific regression parameters]
##                               - k_12,...,k_69 [transition-specific baseline hazards]
##       seed         [numeric]: Single value for random number generation
##
## Output: dat [data frame]: Simulated mstate data in wide format


generate_dat <- function(n, p, X.design = 'norm', pars, seed = 2024){

  #set.seed(seed)
  p.total <- p

  if(X.design == 'norm') X <- matrix(rnorm(n*p.total, 0, 1), ncol = p.total, nrow = n)
  if(X.design == 'bin') X <- matrix(rbinom(n*p.total, 1, 0.5), ncol = p.total, nrow = n)


  # Transition-specific regression coefficients:
  beta_12 <- pars$beta_12
  beta_13 <- pars$beta_13
  beta_24 <- pars$beta_24
  beta_25 <- pars$beta_25
  beta_46 <- pars$beta_46
  beta_47 <- pars$beta_47
  beta_68 <- pars$beta_68
  beta_69 <- pars$beta_69

  # Transition-specific baseline hazards:
  k12 <- pars$k_12
  k13 <- pars$k_13
  k24 <- pars$k_24
  k25 <- pars$k_25
  k46 <- pars$k_46
  k47 <- pars$k_47
  k68 <- pars$k_68
  k69 <- pars$k_69

  dat <- data.frame(X)

  dat <- dat %>%
    mutate(
      hr_12 = exp(X %*% beta_12),               # Individual HR
      hr_13 = exp(X %*% beta_13),
      lambda_12 = hr_12 * k12,                  # Individual hazard
      lambda_13 = hr_13 * k13,
      time_123 = rexp(n, lambda_12+lambda_13),  # Time-to-event drawn from an exp-distribution

      hr_24 = exp(X %*% beta_24),
      hr_25 = exp(X %*% beta_25),
      lambda_24 = hr_24 * k24,
      lambda_25 = hr_25 * k25,
      time_245 = rexp(n, lambda_24+lambda_25),

      hr_46 = exp(X %*% beta_46),
      hr_47 = exp(X %*% beta_47),
      lambda_46 = hr_46 * k46,
      lambda_47 = hr_47 * k47,
      time_467 = rexp(n, lambda_46+lambda_47),

      hr_68 = exp(X %*% beta_68),
      hr_69 = exp(X %*% beta_69),
      lambda_68 = hr_68 * k68,
      lambda_69 = hr_69 * k69,
      time_689 = rexp(n, lambda_68+lambda_69)
    )

  dat <- dat %>%
    mutate(
      # Competing events: CR1 & Death_noCR
      start1 = rep(0, n), # All patients start from active disease set at time = zero
      stop1 = time_123,
      x1 = runif(n),
      CR1 = 1*(x1 <= lambda_12/(lambda_12+lambda_13)),
      Death_noCR = 1-CR1,

      # Competing events: Relapse1 & NRM_CR
      start2 = if_else(CR1 == 1, time_123, NA_real_),
      stop2 = if_else(CR1 == 1, time_123 + time_245, NA_real_),
      x2 = runif(n),
      Relapse1 = if_else(CR1 == 1, 1*(x2 <= lambda_24/(lambda_24+lambda_25)), NA_real_),
      NRM_CR = 1-Relapse1,

      # Competing events: CR2 & RM
      start3 = if_else(CR1 == 1 & Relapse1 == 1, time_123 + time_245, NA_real_),
      stop3 = if_else(CR1 == 1 & Relapse1 == 1, time_123 + time_245 + time_467, NA_real_),
      x3 = runif(n),
      CR2 = if_else(CR1 == 1 & Relapse1 == 1, 1*(x3 <= lambda_46/(lambda_46+lambda_47)), NA_real_),
      RM = 1-CR2,

      # Competing events: Relapse2 & Death_CR2
      start4 = if_else(CR1 == 1 & Relapse1 == 1 & CR2 == 1, time_123 + time_245 + time_467, NA_real_),
      stop4 = if_else(CR1 == 1 & Relapse1 == 1 & CR2 == 1, time_123 + time_245 + time_467 + time_689, NA_real_),
      x4 = runif(n),
      Relapse2 = if_else(CR1 == 1 & Relapse1 == 1 & CR2 == 1, 1*(x4 <= lambda_68/(lambda_68+lambda_69)), NA_real_),
      Death_CR2 = 1-Relapse2
    )

  return(dat)
}


# Parameter settings:

q <- 8 #transitions
p <- 2 #covariates

# Transition-specific regression effects and baseline hazards
pars <- list(beta_12 = c(1.5, rep(0, p-1)), beta_13 = rep(0, p),
             beta_24 = c(1.2, rep(0, p-1)), beta_25 = c(-0.8, rep(0, p-1)),
             beta_46 = rep(0, p), beta_47 = rep(0, p),
             beta_68 = c(1.2, rep(0, p-1)), beta_69 = c(-0.8, rep(0, p-1)),
             k_12 = 0.05, k_13 = 0.05, k_24 = 0.05, k_25 = 0.05,
             k_46 = 0.05, k_47 = 0.05, k_68 = 0.05, k_69 = 0.05)
#ensure enough events for last transition: k_69 = 0.2 or increase N

# Stacked regression parameter
beta <- c(1.5, 0, 1.2, -0.8, 0, 0, 1.2, -0.8, rep(0, q))


# Single simulated dataset: only binary covariates
MS <- generate_dat(n = 500, p = p, pars = pars, X.design = 'bin')

sim_data_aml <- subset(MS, select = c('stop1', 'stop2', 'stop3', 'stop4',
                                      'CR1', 'Death_noCR', 'Relapse1', 'NRM_CR',
                                      'CR2', 'RM', 'Relapse2', 'Death_CR2',
                                      'X1', 'X2'))


# Save sim_data_aml:
usethis::use_data(sim_data_aml, overwrite = TRUE)
