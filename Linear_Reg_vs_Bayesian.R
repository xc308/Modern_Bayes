#================================
# Regression Models vs Bayesian 
#================================

#--------
# Model
#-------

# for running group: 
  # E[y | x] = beta1 + beta3 * X3

# for exercise group:
  # E[y | x] = (beta1 + beta2) + (beta2 + beta4) * X3



#-------
# Data
#-------

# running is 0, 1 is aerobic
x1<-c(0,0,0,0,0,0,1,1,1,1,1,1)

# age
x2<-c(23,22,22,25,27,20,31,23,27,28,22,24)

# change in maximal oxygen uptake
y<-c(-0.87,-10.74,-3.27,-1.97,7.50,
     -7.25,17.05,4.96,10.40,11.05,0.26,2.51)

x3 <- x2  # x3 = age
# [1] 23 22 22 25 27 20 31 23 27 28 22 24

x2 <- x1  # x2 is a 2 level indicator variable
# [1] 0 0 0 0 0 0 1 1 1 1 1 1

# index of each person
x1 <- seq(1:length(x2))
# [1]  1  2  3  4  5  6  7  8  9 10 11 12

x1 <- rep(1, 12)

# x4 = x2*x3
x4 <- x2 * x3
# [1]  0  0  0  0  0  0 31 23 27 28 22 24


# design matrix 
X <- cbind(x1, x2, x3, x4)
head(X)


#-------------------------------
# OLS estimates of unknown beta
#-------------------------------

# beta_ols = (X^T X)^{-1} X^t y 
beta_ols <- solve(t(X) %*% X) %*% (t(X) %*% y)
#           [,1]
#x1 -51.2939459
#x2  13.1070904
#x3   2.0947027
#x4  -0.3182438



# alternative way using lm model fitting
fit.ols <- lm(y ~ X[, 1] + X[, 2] + X[, 3] + X[, 4])
fit.ols$coefficients
# (Intercept)      X[, 1]      X[, 2]      X[, 3]    X[, 4] 
# -51.2939459          NA  13.1070904   2.0947027   -0.3182438


sum(fit.ols$residuals^2)
#  68.33981 same as below SSR_beta


# for running group: 
# E[y | x] = beta1 + beta3 * X3
          # = -51.2939459 +  2.0947027 X3


# for exercise group:
# E[y | x] = (beta1 + beta2) + (beta3 + beta4) * X3
          # = (-51.2939459 + 13.1070904) + (2.0947027-0.3182438) X3


#------------------------------
## SSR sum of squares residuals
#------------------------------
## SSR(beta) = sum_{i = 1}^n (yi - beta^t %*% xi)^2
  ## (y - X beta)^t (y - X beta)

SSR_beta <- t(y - X %*% beta_ols) %*% (y - X %*% beta_ols)


#-------
# unbiased estimator for sigma2
#-------
n = length(y)
p = ncol(X)
sigma2_ols <- 1/(n-p-1) * SSR_beta
# 9.762831
sigma2_ols <- c(sigma2_ols)


#-----
# the sampling variance of beta_ols: var[beta_ols]
# (t(X)X)_{-1} sigma2
#-----
XtX <- solve(t(X) %*% X) 
Var_beta_ols <- sigma2_ols * XtX
# [1:4, 1:4]


#----------------------------------------
# Bayesian analysis of linear regression
#----------------------------------------

## unknown parameters of interest
  # beta
  # sigma2

## use full conditionals of each of them + Gibbs

# likelihood: y | X, beta, sigma2 ~ MVN(m_beta_ols, v_beta_ols) a function of beta
    # prior: beta ~ MVN(beta0, Sigma0)
        # posterior: 
          # beta | y, X, sigma2 ~ MVN(m_beta, V_beta)
              # var_beta_inv = (t(X) %*% X) /sigma2
              # V_beta = (Sigma0_inv + var_beta_inv)^{-1}
              # m_beta = (beta0 * Sigma0_inv + (t(X)%*%y) / sigm2) * V_beta

# likelihood: y | X, beta, sigma2 ~ MVN(;)
    # prior : gam = 1/sigma2 ~ Gamma(a, b)
        # posterior:
            # sigma2 | y, X, beta ~ InvGamma(A, B)
                # A = a + n / 2
                # B = (b + SSE) / 2


#-----
# Parameters in the prior
#-----
# beta0 = beta_ols    # ols estimate of beta
# Sigma0_inv = {var[beta_ols]}^{-1} / n # unit informaion
            #= 1/n * t(X)%*%X / sigma2
            # sigma2 = sigma2_ols plug in

# a = 1
# b = sigma2_ols 


beta0 = beta_ols
Sigma0_inv = 1/n * (t(X)%*%X) /sigma2_ols  # [1:4, 1:4]

a = 1
b = sigma2_ols


#--------
# Gibbs sampler
#--------
library(mvtnorm)

Gibbs_LinearReg <- function(sigma2_start , n.sim, Data = y) {
  
  n = length(Data)
  
  sigma2 = sigma2_start
  #sigma2 = sigma2_ols
  
  BETA <- SIGMA2 <- NULL
  for (i in 1:n.sim) {
    
    # sample beta | y, X, sigma2 ~ MVN(m_beta, V_beta)
    V_beta = solve(Sigma0_inv + t(X) %*% X / sigma2) # 4*4
    m_beta = V_beta %*% (Sigma0_inv %*% beta0 + t(X) %*% y / sigma2)  # 4*1
    
    beta = rmvnorm(1, mean = m_beta, sigma = V_beta)
    # str(beta)  num [1, 1:4]
    beta = t(beta) # 4*1
    
    # sample sigma2 | y, X, beta ~ InvGamma(A = a + n/2, B = (b + SSE_beta)/2)
    
    A = a + n / 2
    SSE_beta = t(y - X %*% beta) %*% (y - X %*% beta) # scalar
    B = (b + SSE_beta) / 2
    
    sigma2 = 1/rgamma(1, shape = A, rate = B)
    
    
    # update the unknowns for this run
    BETA <- rbind(BETA, t(beta)) # 1*4 
    SIGMA2 <- rbind(SIGMA2, sigma2)
  }
  return(list(beta = BETA, sigma2 = SIGMA2))
}



#-------------
# Test on data
#-------------

Res <- Gibbs_LinearReg(sigma2_start = sigma2_ols, n.sim = 1e3, Data = y)

Res_beta <- Res$beta
Res_sigma2 <- Res$sigma2


#------
# Plots
#------

plot(1:1e3, Res_sigma2, type = "l", 
     main = bquote("Trace plot of "~ beta))


head(Res_beta)
par(mfrow = c(1, 4))
for(idx in 1:4) {
  plot(1:1e3, Res_beta[, idx], type = "l",
       main = bquote("Trace plot of "~ beta[.(idx)]))
}



#---------------
# Summary stats
#---------------

beta_pst_mean <- apply(Res_beta, 2, mean)
sigma2_pst_mea <- mean(Res_sigma2)

beta_CI <- apply(Res_beta, 2, function(x) quantile(x, c(0.025, 0.975)))
str(t(beta_CI))

sigma2_CI <- quantile(Res_sigma2, c(0.025, 0.975))

beta_stats <- cbind(beta_pst_mean, t(beta_CI))

sigma2_stats <- cbind(sigma2_pst_mea, t(sigma2_CI))

Sum_stats <- rbind(beta_stats, sigma2_stats)
rownames(Sum_stats) <- c("beta1", "beta2", "beta3", "beta4", "sigma2")
colnames(Sum_stats) <- c("pst_mean", "2.5%", "97.5%")


library(xtable)
xtable(Sum_stats)


Sum_stats_ols <- cbind(Sum_stats, rbind(beta_ols, sigma2_ols))
colnames(Sum_stats_ols)[4] <- c("ols")




#---------------
# g-prior model
#---------------

# likelihood:
  # y|X, beta, sigma2 ~ MVN((t(X)X)_inv t(X)y, t(X)X)_inv sigma2 )

# g-prior on beta:
  # beta ~ MVN(0, g*sigma2 (t(X)X)_inv)

# prior on simga2 or gamma = 1/sigma2
  # gamma ~ Gamma(nu0/2, nu0 sigma0^2 / 2)

# marginal posterior density of beta:
  # beta|y, X, sigma2 ~ MVN(m_beta, V_beta)
     # m_beta = g/(g+1) (t(X)X)^{-1} t(X)y
     # V_beta = g/(g+1) (t(X)X)^{-1} sigma2

# marginal posterior density of sigma2:
  # sigma2|y, X, beta ~ InvGamma(a, b)
    # a = 1/2(nu0 + n)
    # b = 1/2(nu0 sigma02 + SSR_g)

    # SSR_g = t(y)(I - g/1+g X(t(X)X)^{-1}t(X))y


# parameter in the priors:
  # g = n = length(y)
  # nu0 = 1
  # sigma0^2 = sigma_ols^2 = 1/n-p SSR(beta) = 8.54

    # n = length(y)
    # p = ncol(X)
    # sigma2 <- 1/(n-p) * SSR_beta
      # 8.542477


# joint posterior distribution (beta, sigma2 | y, X) 
  # can be made with Monte Carlo approximation 
  # not necessarily with Gibbs

# 1. Sample 1/sigma2 ~ Gamma(a, b)
# 2. sample beta ~ MVN(m_beta, V_beta)


y
X
# parameters in prior
g <- n <- length(y)
p <- ncol(X)
nu0 <- 1
sigma2_0 <- sigma2_ols <- 1/(n-p) * SSR_beta

S <- 1e3  # number of MC simulations

H_g <- g/(g+1) * X %*% solve(t(X) %*% X) %*% t(X)
str(H_g) # num [1:12, 1:12]
SSR_g <- t(y) %*% (diag(1, n) -  H_g) %*% y
# 123.9618


# sample from sigma2 ~ IG(a, b)
a <- (nu0 + n) / 2 
b <- (nu0 * sigma2_0 + SSR_g) / 2

sigma2 <- 1/rgamma(1, shape = a, rate = b)


# sample beta ~ MVN(m_beta, V_beta)
m_beta <- g/(g+1) * solve(t(X) %*% X) %*% (t(X) %*% y)
#           [,1]
#x1 -47.3482578
#x2  12.0988527
#x3   1.9335717
#x4  -0.2937635

V_beta <- g/(1+g) * sigma2 * solve(t(X) %*% X)
str(V_beta)
# num [1:4, 1:4] 

library(mvtnorm)
beta <- rmvnorm(1, mean = m_beta, sigma = V_beta)


U <- rmvnorm(1, mean = rep(0, length(m_beta)), 
             sigma = diag(1, nrow = length(m_beta)))
str(U)
# num [1, 1:4]

# U %*% chol(V_beta)
#            x1        x2        x3        x4
# [1,] -5.438595 -5.626524 0.2455127 0.1795727

beta_alternative <- m_beta + t(U %*% chol(V_beta) ) 
























