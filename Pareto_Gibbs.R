#==========================
# Pareto Example for Gibbs
#==========================

install.packages("dplyr")
library(dplyr) # for data frame row selection
library(coda) # for N_eff


#-------------------------------------
# Gibbs Sampler for Pareto likelihood
#-------------------------------------
# This is a function for sampling two parameters alpha, c in Pareto example. 

# Likelihood: pareto(alpha, c | Data) 
# Prior: improper prior I(alpha, c > 0)
# Posterior: No closed form to sample directly

# Using Gibbs to sample from conditional distribution
  # alpha|c, Data ~ Gamma(alpha| n+1, Sum_i log(Data_i) - n * log c)
  # c|alpha, Data ~ Mono(c|n * alpha + 1, x_min)


# Arg: 
  # alpha_start: starting value for parameter alpha 
  # c_start: starting value for parameter c
  # n.sims: total number of simulations
  # Data: data in vector form

# Return:
  # a matrix n.sims * 2
  # the 1st col is the Gibbs samples (MCMC samples) for par alpha
  # the 2nd col is the Gibbs sampes (MCMC samples) for par c

Pareto_Gibbs <- function(alpha_start, c_start, n.sims, Data) {
  n = length(Data)
  sum_log = sum(log(Data))
  Dt_min = min(Data)
  
  alpha_c_post <- matrix(NA, nrow = n.sims, ncol = 2) 
  alpha_c_post[1, ] <- c(alpha_start, c_start)
  
  for (i in 2:n.sims) {
    # sample alpha ~ Gamma(n+1, sum_i log(xi) - n*log(c))
    alpha_c_post[i, 1] <- rgamma(1, shape = n+1, 
                                 rate = sum_log - n*log(alpha_c_post[i-1, 2]))
    
    # sample c ~ Mono(n*alpha + 1, xmin)
    u <- runif(1)
    alpha_c_post[i, 2] <- Dt_min * u^{1/ (n * alpha_c_post[i, 1] + 1)}
  }
  return(alpha_c_post)
}



#-----------------------------
# Data: 
# World city population
#-----------------------------
Wrld_city <- read.csv("world_city_pareto.csv", header = T)
head(Wrld_city)


Wrld_city_eff <- filter(Wrld_city, Wrld_city$pop2023 > 750000)

Pop2023 <- Wrld_city_eff$pop2023

min(Pop2023)
# [1] 750097

# mean of Pareto distribution 
  # alpha * xmin / (alpha - 1)

mean(Pop2023)
# [1] 2654390

find_alpha <- function(xmin, b) {
  alpha = b / (b - xmin)
  return(alpha)
}

find_alpha(xmin = 750097, b = 2654390)
# [1] 1.393898

quantile(Pop2023)
#     0%        25%        50%        75%       100% 
#    41.0   971703.2  1364431.0  2526683.2  37194104.0

# c = 900000

## test on ParetoGibbs
a_c_post <- Pareto_Gibbs(alpha_start = 1, c_start = 900000, n.sims = 1e4, Data = Pop2023)
alpha_post <- a_c_post[, 1]
c_post <- a_c_post[, 2]


#----------
# Plots:
# Checking for convergence, MCMC samples independence
#----------
## check for convergence or not 
par(mar = c(4.5, 4.5, 6, 2))
plot(1: 1e4, alpha_post, type = "l",
     xlab = "MCMC Iterations", 
     ylab = expression(paste("posterior of ", alpha)),
     #main = bquote(mu ~ ":" ~ .(mu2) ~ mm),
     main = bquote(atop("Posterior of "~alpha, 
                        atop("lag-1 acf:"~ .(round(acf_a_lag1, 2)),
                             "N_eff:"~ .(round(N_eff, 2))))),
     cex = 0.55)


## check for independence / correlation of MCMC samples
  # acf of MCMC samples, the smaller, the more indpendent
  # N_eff: if similar to n.sims, indipendent 

acf_alpha <- acf(alpha_post)
acf_a_lag1 <- acf_alpha$acf[2]
# [1] 0.01980889

N_eff <- effectiveSize(alpha_post)
# 9610.556 



#------------
# Gibbs sampler including plot
#------------

Pareto_Gibbs_plt <- function(alpha_start, c_start, n.sims, Data, plot = T) {
  n = length(Data)
  sum_log = sum(log(Data))
  Dt_min = min(Data)
  
  alpha_c_post = matrix(NA, nrow = n.sims, ncol = 2)
  alpha_c_post[1, ] = c(alpha_start, c_start)
  
  for (i in 2:n.sims) {
    # sample alpha ~ gamma(n+1, sumlogxi - n*log(c))
    alpha_c_post[i, 1] <-
      rgamma(1, shape = n+1, rate = sum_log - n*log(alpha_c_post[i-1, 2]))
    
    # sample c ~ Mono(n*alpha + 1, xmin)
      # 1. u~ U(0, 1)
      # 2. c = F_inv (u)
    u <- runif(1)
    alpha_c_post[i, 2] <- Dt_min * u^{1/(n*alpha_c_post[i, 1] + 1)}
  }
  
  # lag-1 acf
  lag1_acf_alpha <- acf(alpha_c_post[, 1])$acf[2]
  lag1_acf_c <- acf(alpha_c_post[, 2])$acf[2]
  
  # effective sample size for MCMC samples
  library(coda)
  N_eff_alpha <- effectiveSize(alpha_c_post[, 1])
  N_eff_c <- effectiveSize(alpha_c_post[, 2])
  
  if (plot) {
    # trace plot for checking convergence
    par(mfrow = c(1, 2), mar = c(4.5, 4.5, 6, 3))
    plot(1:n.sims, alpha_c_post[, 1], type = "l",
         xlab = "MCMC iterations", 
         ylab = expression(paste("Posterior of ", alpha)),
         main = bquote(atop("Posterior of " ~ alpha,
                            atop("lag-1 acf:" ~ .(lag1_acf_alpha),
                                 "N_eff:" ~ .(N_eff_alpha)))),
         cex = 0.55)
    
    plot(1:n.sims, alpha_c_post[, 2], type = "l",
         xlab = "MCMC iterations", 
         ylab = expression(paste("Posterior of ", c)),
         main = bquote(atop("Posterior of " ~ c,
                            atop("lag-1 acf:" ~ .(lag1_acf_c),
                                 "N_eff:" ~ .(N_eff_c)))),
         cex = 0.55)
  }
  return(alpha_c_post)
}

alf_c_post <- Pareto_Gibbs_plt(alpha_start = 1, c_start = 800000, n.sims = 1e4, 
                 Data = Pop2023, plot = T)

mean(alf_c_post[, 1])
# [1] 1.163881

mean(alf_c_post[, 2])
# [1] 749286
