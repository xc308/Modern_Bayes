#==============================================
# MonteCarlo_Importance_Sampling_Compact_Spport
#==============================================

# evaluate I = int_{-inf}^{inf} exp(-x^4) dx
# Method : Monte Carlo Integration
  # E_f[h(x)] = int h(x) f(x) dx
  # where f(x) is easy to sample from 

# to construct an f(x) easy to sample from 
# let y = sqrt(2) x^2
# by substitution derivation
# I = 2^{-5/4} int_0^{inf} sqrt(2*pi / y) 2*phi(y) dy
  # 2*phi(y) is P(|Y| < y) the density of std.N(0, 1) 
    # truncated from postive Real line
  # computation can 1. sample from N(0, 1)
                    # 2. take abs of these sample values


#-------------------
# mc sampling target
#-------------------
# mc1: truncated (folded) normal sampling distribution
mc1 <- function(x) {
  sqrt(2 * pi / abs(x))
}

draws <- rnorm(1e5)
values1 <- sapply(draws, mc1)
values1 <- values1 * 2^{-5/4}
hist(values1, xlab = paste0("sample values from ", expression(f(x))),
     main = "Samples from folded N(0, 1)")


# mc2: importance sampling from N(0, 1) 
mc2 <- function(x) {
  exp(-x^4) / dnorm(x)
}

draws2 <- rnorm(1e5)
values2 <- sapply(draws2, mc2)
hist(values2, xlab = paste0("sample values from ", expression(f(x))),
     main = "Samples from N(0, 1)")


# mc3: uniform sampling distribution

# since when x = 10, exp^{-10^4} = 0
# means exp^{-x^4} decreases to 0 quickly
# so can truncate the tail and turn I into a compact support
# I = int_{-10}^{10} exp(-x^4)
  # = int_{-10}^{10} 20*exp(-x^4) / (10-(-10)) dx
  # which is an average of 20*exp(-x^4) wrt unif(-10, 10)

# so sample from Unif(-10, 10)
# evaluate 20*exp(-x^4) and take average

mc3 <- function(x) {
  20 * exp(-x^4)
}

draws3 <- runif(1e5, -10, 10)
values3 <- sapply(draws3, mc3)

hist(values3, xlab = paste0("sample values of ", expression(f(x))),
     main = "Samples from Unif(-10, 10)")


#----------------------------------
# Root Mean Square Errors of samples
#----------------------------------

RMSE <- function(values) {
  sd(values) / sqrt(length(values))
}


#--------------------------------------------------
# Means of sampled values from f(x) using 3 methods
#--------------------------------------------------
# Method 0, 
  # closed from derived using substituion u = x^4 and
  # Gamma(x |a, b)

# true value of I 
m0 <- gamma(1/4) / 2 
# [1] 1.812805

# Method 1,
  # use Monte Carlo integration
  # let y = sqrt(2)*x^2
  # substitue, 
  # f(x) = 2phi(x), a folded N(0,1), only positive real line

m1 <- mean(values1) 
# values1 <- values1 * 2^{-5/4}


# Method 2,
  # importance sampling
  # sample from f(x) = N(0, 1)
  # evaluate exp^{-x^4} / dnorm(x)

m2 <- mean(values2)


# Method 3,
  # truncated support or compact support
  # I = int_{-10}^{10} exp^{-x^4} 
    # = int_{-10}^{10} 20*exp^{-x^4} / (10 - (-10))
  # sample from f(x) = Unif(-10, 10)
  # evaluate 20*exp^{-x^4}

m3 <- mean(values3)

MEANS <- c(m0, m1, m2, m3)

RMSEs <- c(NA, RMSE(values1), RMSE(values2), RMSE(values3))

mc.df <- data.frame(Means = MEANS, RMSEs = RMSEs, 
           row.names = c("True value", "Folded N(0, 1)",
                         "Full N(0, 1)", 
                         "Truncated Unif(-10, 10)"))

library(xtable)
xtable(mc.df)

## ideally, histogram of sampled values using each method 
# is centered around the true mean 1.813
# but the sampled values scattered in different ways






