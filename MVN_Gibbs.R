#===============
# The MVN Model
#===============

# yi: p*1 vector
# yi | theta, S0 ~ MVN(theta, S0^{-1})

# likelihood: 
  # y1, ..., yn |theta, Sigma ~ MVN(theta, Sigma)

# priors:
  # theta ~ MVN(mu_0, Lambda_0)
  # Sigma ~ InvWishart(nu_0, S_0) // a mv form of gamma

# joint posterior:
  # (theta, Sigma | y1,...,yn) unknown & untractable


## Consider full conditional
# if Sigma is known, 
  # then likelihood:  y1,..., yn | theta ~ MVN(theta, Sigma)
      #prior on theta:       theta ~ MVN(mu_0, Lambda_0)
          # posterior: theta | y1,...,yn, Sigma ~ MVN(mu_n, L_n)

# if theta is known, 
  # then likelihood: y1,..., yn | Sigma ~ MVN(theta, Sigma)
      # prior on Sigma: Sigma ~ InvWishart(nu_0, S0^{-1})
          # posterior: Sigma | y1,...,yn, theta ~ InvWishart(nu_n, S_n^{-1})


# Given starting value for theta and Sigma, 
# Gibbs sample from above full conditionals 


#-------------------------------
# Basic computational knowledges
#-------------------------------

## MVN distributions:
#--------------------

install.packages("mvtnorm")
library(mvtnorm)

## Simulate ##
  # simulate one MVN random vector
rmvnorm(1, mean = c(0, 2), sigma = diag(2))
#         [,1]     [,2]
# [1,] 1.691275 1.615846
# each row corresponds to one sample
# above is one row, ie one sample

str(rmvnorm(1, mean = c(0, 2), sigma = diag(2)))
# num [1, 1:2] -1.5 3.49
str(c(rmvnorm(1, mean = c(0, 2), sigma = diag(2))))
# num [1:2] -0.16 1.86


str(rmvnorm(3, mean = c(0, 2), sigma = diag(2)))
# num [1:3, 1:2]
rmvnorm(3, mean = c(0, 2), sigma = diag(2))
#           [,1]     [,2]
#[1,]  0.5435996 1.720621
#[2,]  1.4296302 1.273166
#[3,] -1.8272226 2.414615



## Evaluate ##
  # evaluate the mvnormal density at a sigle quantile (value) 
    # if bivariate quantile x = (x1, x2) = c(1,1)
dmvnorm(x = c(1, 1), mean = c(0, 2), sigma = diag(2))
# [1] 0.05854983


## Wishart distribution:
#------------------------

# stats package contains functions for evaluating 
                # and simulating from a Wishart density

library(stats)

## Simulate ##
# can simulate a single Wishart distriburted matrix 
nu_0 <- 2
sigma_0 <- diag(2)
rWishart(1, df = nu_0, Sigma = sigma_0)[, , 1]
# Sigma:	
# positive definite (p×p) “scale” matrix, 
  # the matrix parameter of the distribution

# Value:
  # a numeric array of dimension p×p×n, 
  # where each R[,,i] is a positive definite matrix,
  # a realization of the Wishart distribution 


#            [,1]      [,2]
# [1,]  0.02454644 -0.060081
# [2,] -0.06008100  3.588225


rWishart_mat <- rWishart(3, df = nu_0, Sigma = sigma_0)

rWishart_mat
# , , 1

#        [,1]      [,2]
#[1,] 4.0827158 0.8191552
#[2,] 0.8191552 0.4904935

#, , 2

#        [,1]      [,2]
#[1,] 5.1570517 0.9274856
#[2,] 0.9274856 0.8977809

#, , 3

#         [,1]       [,2]
#[1,]  0.09478939 -0.5087137
#[2,] -0.50871374  4.0967433


#==================================================
# setting parameters in the prior -- hyper-priors
#==================================================

# theta ~ MVN(mu_0, Lambda_0)
# Sigma ~ InvWishart(nu_0, S_0^{-1})

# mu_0: the exam is designed to have average scores of ~ 50 out of 100
  # mu_0 = [mu_01, mu_02] = [50, 50]

# Lambda_0:
  # since score has to be in [0, 100]
    # so set variance or std to ensure majoirty of data fall into [0,100]
    # by normality, 
        # 68% data fall in 50 +/- std, 
        # 95% fall into 50 +/- 2std
    # so set std = 25, then P(theta_j NOT in [0, 100]) = 0.05
    # var of Lambda_0 = [Lambda_01, Lambda_02]^T = [25^2, 25^2]
    
    # for off-diag Lambda_012 and Lambda_021, 
      # due to the 2 tests are similar content, set correlation = 0.5
      # i.e., Lambda_012 / Lambda_01 = 0.5
      # Lambda_012 = 312.5 = Lambda_021

    # so diag(Lambda_0) = [Lambda_01, Lambda_02]^T = [25^2, 25^2]
    # off-diag(Lambda_0) = Lambda_012 = 312.5 = Lambda_021


# S0:
  # set S0 = Lambda_0

# nu_0:
  # to set Sigma centerend around Lambda_0 so nu_0 = p+2 = 4
  # E[Sigma] = (1/nu_0 - p - 1) * S0 = S0

  # nu_0 = p+2 = 4


#-----------
# Data set
#-----------



structure()
# returns the given object with further attributes set.
head(Y)

dim(Y)
# [1] 22  2

# get ybar for each test, columnwise
ybar <- apply(Y, 2, mean)
# pretest posttest 
# 47.18182 53.86364
str(ybar)

Sigma <- cov(Y)
#          pretest posttest
# pretest  182.1558 148.4069
# posttest 148.4069 243.6472


#------------------
# set hyper-priors
#------------------

mu_0 <- c(50, 50)
Lamb_0 <- matrix(c(625, 312.5, 312.5, 625), nrow = 2, ncol = 2)
S_0 <- Lamb_0
nu_0 <- 4

str(mu_0)
#num [1:2] 50 50

#--------------
# Gibbs Sampler
#--------------

# This function is a Gibbs sampler for theta and Sigma
  # in MVN(theta, Sigma)

# Arg:
  # theta_start
  # Sigma_start



MVN_Gibbs <- function(theta_start, Sigma_start, n.sim, Data) {
  
  ## check dimension of theta_start and Sigma_start
  #if (dim(Sigma_start) != c(ncol(Data), ncol(Data))) {
    #print("Stop! Sigma_start matrix dimension is not right!")
  #}
  
  n <- nrow(Data)
  Sigma <- Sigma_start
  
  THETA <- NULL
  SIGMA <- NULL
  
  for (i in 1:n.sim) {
    
    # sample theta|data, Sigma
    Lamb_n = solve(solve(Lamb_0) + n * solve(Sigma)) # Sigma_inv
    mu_n = Lamb_n %*% (solve(Lamb_0) %*% mu_0 + n * solve(Sigma) %*% ybar) 
    
    theta <- rmvnorm(1, mean = mu_n, sigma = Lamb_n)  # num [1, 1:2] df
    theta <- c(theta) # num [1:2] vector
    
    
    # sample Sigma|data, theta ~ InvWishart(nu_0, S0^{-1})
    nu_n = nu_0 + n
    S_theta = (t(Data) - theta) %*% t(t(Data) - theta) # 2*2
    Sn_inv = solve(S_0 + S_theta)
    
    Sigma = solve(rWishart(1, df = nu_n, Sigma = Sn_inv)[, , 1]) # solve Wishart for InvWishart
    
  
    # collect theta and Sigma 
    THETA <- rbind(THETA, theta)  # df 
    SIGMA <- rbind(SIGMA, c(Sigma))
  }
  return(list(theta = THETA, Sigma = SIGMA))
}


Sigma_start <- matrix(c(150, 100, 100, 200), nrow = 2)
RES <- MVN_Gibbs(theta_start = c(50, 50), Sigma_start = Sigma_start, 
          n.sim = 500, Data = Y)

RES$theta
RES$Sigma


Sigma <- cov(Y)
#          pretest posttest
# pretest  182.1558 148.4069
# posttest 148.4069 243.6472



str(rbind(1, 2))







































