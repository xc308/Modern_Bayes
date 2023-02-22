#=======================
# Two-components Mixture
#=======================

install.packages("mice")
library(mice)  # selfreport data set containing measured heights
library(dplyr)

install.packages("LaplacesDemon")
library(LaplacesDemon)


#-----------
# Data set
#-----------
head(selfreport)

self_krul <- filter(selfreport, src == "krul")
Heights <- self_krul[, "hm"]

str(Heights) 
# num [1:1257]
# heights containing M/F

Heightsex <- select(self_krul, c("sex", "hm"))
H_F <- filter(Heightsex, sex == "Female") # 695 obs
H_M <- filter(Heightsex, sex == "Male") # 562 obs

mean(H_F[, 2]) # [1] 167.998
mean(H_M[, 2]) # [1] 181.4133
mean(H_M[, 2]) - mean(H_F[, 2]) # [1] 13.41536

round(mean(Heights), 0)
# [1] 173.9959
sqrt(var(Heights))
# [1] 10.46119


#------------------------------------------
# Gibss sampler for two-components mixture
#------------------------------------------

rbern(1, 0.5) # generate 1 sample from p^x (1-p)^(1-x), x = 0,1
              # generate 0 or 1 according probability 0.5

## Herirachical Model
# Xi|mu, Z ~iid N(mu_zi, lambda^{-1}), 
  # Z1, ... Zn ~ iid Bern(Pi)
    # mu ~ N(m, l^{-1})
    # Pi ~ Beta(a, b)

## Prior settings hyperparameters
a = b = 1 # a/a+b = 0.5, prior on propotion Pi is 50%:50%
m = round(mean(Heights), 0) # 174

sigma = round(sqrt(var(Heights)), 0)  # 10
lambda = 1/(sigma^2)

s = round(mean(H_M[, 2]) - mean(H_F[, 2]), 0) # 13
l = 1/(s^2)


## Starting values setting for desired unkowns
Pi = .5   # starting proportion for two components is 50:50
mu0 = mu1 = m 
Z = rbern(n = length(Heights), prob = Pi)


#------------------------------------------
## Gibbs Sampler for Two-components Mixture
#------------------------------------------
Mixture_Gibbs <- function(Pi_start, mu0_start, mu1_start, Z_start, n.sim, burn.in, Data) {
  Z = Z_start
  n = length(Data)
  
  Res <- matrix(NA, nrow = n.sim, ncol = (1 + 2))
  for (i in 1:n.sim) {
    n1 = sum(Z)
    n0 = n - n1
    
    # update Pi
    Pi_pst <- rbeta(1, n1 + a, n0 + b)
    
    # update mu0, mu1
    X <- data.frame(cbind(Z, Data))
    X0 <- X[Z == 0, ]
    X1 <-  X[Z == 1, ]
    
    L0 <- n0 * lambda + l
    M0 <- (m * l + lambda * sum(X0[, 2])) / L0
    
    mu0_pst <- rnorm(1, M0, sd = sqrt(1/L0))
    
    L1 <- l + n1 * lambda
    M1 <- (m * l + lambda * sum(X1[, 2])) / L1
    
    mu1_pst <- rnorm(1, M1, sd = sqrt(1/L1))
    
    # update Z
    for (j in 1:n) {
      D1 = dnorm(Data[j], mu1_pst, sqrt(1/lambda))
      D0 = dnorm(Data[j], mu0_pst, sqrt(1/lambda))
      
      alpha1 = Pi_pst * D1
      alpha0 = Pi_pst * D0
      
      prob = alpha1 / (alpha1 + alpha0)
      
      Z[j] <- rbern(1, prob = prob)
    }
    
    ## collect posteriors 
    Res[i, ] <- c(Pi_pst, mu0_pst, mu1_pst)
  }
  return(Res[burn.in:n.sim, ])
}


#-------------
# Test on data
#-------------

## Starting values setting for desired unkowns
Pi = .5   # starting proportion for two components is 50:50
mu0 = mu1 = m 
Z = rbern(n = length(Heights), prob = Pi)


Res <- Mixture_Gibbs(Pi_start = Pi, mu0_start = mu0, mu1_start = mu1, 
                     Z_start = Z, n.sim = 1e3, burn.in = 1, Data = Heights)



#-----------
# Plots
#-----------
## for convergence
  # trace plot
  # running average

## for samples independence
  # lag-1 acf
  # N_eff 
# in coda package

## estimated density for each of desired parameters 
# posterior mean
# HPD in coda for Credible intervals















#-------
# Try 
#-------
# try for sum_{i:Zi == 0} Xi in M0
str(Heights)
# num [1:1257]

length(Z)
x <- data.frame(cbind(Z, Heights))
x0 <- x[Z==0, ]
x1 <- x[Z==1, ]
head(x0)
str(x0)
# 'data.frame':	649 obs. o
# > head(x)
#  Z Heights
#1 0   172.3
#2 0   169.4
#3 0   165.9
head(x0[, "Heights"])
str(x0[, "Heights"])
# num [1:649] 
sum(x0[, "Heights"])

X <- data.frame(cbind(Z, Heights))
X0 <- X[Z == 0, ]
X1 <-  X[Z == 1, ]
head(X0)

sum(X0[, 2])





# try density of norm when feed in 1 data point as quantile
dheights <- dnorm(Heights[1], 174, 10)
str(dheights)
# num [1:1257]
z <- rbern(1, dheights)

