#========================
# Metropolis algo in GLM
#========================

# Using know Normal model and normal prior with known variance

# y1,...,yn | theta ~ iid N(theta, sigma2); sigma2 known
  # theta ~ N(mu_0, tau_0), 
    # theta | y1,...,yn ~ N(mu_n, tau_n)
      
      # tau_n = (1/tau_0 + n/sigma2)^{-1}  # variance
      # mu_n = (mu_0 * 1/tau_0 + ybar * n/sigma2) * tau_n


# so the posterior theta | y1,...,yn can be computable
  # p(theta | y1,...,yn) = dnorm(mean = mu_n, sd = sqrt(tau_n))


y <- c(9.37, 10.18, 9.16, 11.60, 10.33)
ybar <- mean(y)
n <- length(y)
#mean(y) # [1] 10.128
#var(y) # [1] 0.93047

sigma2 <- 1
tau_0 <- 100
mu_0 <- 5

tau_n <- (1/tau_0 + n/sigma2)^{-1}
sqrt(tau_n)
mu_n <- (mu_0 * 1/tau_0 + ybar * n/sigma2) * tau_n

# posterior density of theta | y1, .., yn
p(theta|y1, ..., yn) = dnorm(x, mean = mu_n, sd = sqrt(tau_n))


## Now assume we don't know the form of the posterior p(theta|y1, ..., yn)
  # so has to resort to the Metropolis algo
delta2 <- 2

Metropolis_NN <- function(theta_0, n.sim, delta2, Data) {
  set.seed(1)
  
  THETA <- NULL
  theta <- theta_0
  count <- 0
  for (i in 1:n.sim) {
    
    # sample theta_star ~ N(theta^(s), delta2)
    theta_star <- rnorm(1, mean = theta, sd = sqrt(delta2))
    
    
    # compute logr 
    logr = sum(dnorm(Data, theta_star, sqrt(sigma2), log = T) - 
      dnorm(Data, theta, sqrt(sigma2), log = T)) +
      dnorm(theta_star, mu_0, sqrt(tau_0), log = T) - 
      dnorm(theta, mu_0, sqrt(tau_0), log = T)
    
    # u ~ unif(0, 1)
    u <- runif(1)
    if(log(u) < logr) {
      # update theta using theta_star
      theta <- theta_star
      count <- count + 1
    }
  
    # update
    THETA <- rbind(THETA, theta)
  }
  AC_ratio <- count / n.sim
  return(list(theta = THETA, AC_ratio = AC_ratio))
}

Theta_pst <- Metropolis_NN(theta_0 = 0, n.sim = 1e4, Data = y)

par(mfrow = c(1,1))
plot(1:1e3, Theta_pst, type = "l")

hist(Theta_pst[10:1e4], freq = F, breaks = 1e2)
lines(density(Theta_pst), col = "red")


#------------------
# Selecting delta2: ideal AC_ratio : 20%-50%
#------------------

# delta2 = 1/32, 1/2, 2, 32, 64

delta <- c(1/32, 1/2, 2, 32, 64)
sapply(delta, function(x) Metropolis_NN(theta_0 = 0, n.sim = 1000, 
                                    delta2 = x, Data = y)$AC_ratio)


# [1] 0.816 0.579 0.349 0.098 0.068
# want AC_ratio = 20% - 50% 
  # when delta2 = 32, 64, AC_ratio is too low with more than 90% RJ
  # when delta2 = 1/32, 1/2 AC_ratio is too high with most of them AC
  # both of these two scenarios end up with highly correlated MC values
    # either correlated with theta_start or theta_{s-1}

  # so choose delta2 = 2 with AC_ratio = 35% 




#==============================
# Metropolis for Poisson model
#==============================
# a sample of 52 female song sparrows was studied over
# a summer and their reproductive activities were recorded

# The boxplot indicates that two-year-old birds in this population
# had the highest median reproductive success, with the number of offspring
# declining beyond two years of age

# Aim: 
  # to understand the relationship between age and reproductive success,
  # or to make population forecasts for this group of birds.

# Model construction:
  # The number of off-spring of each bird is non-neg {0, 1, 2, ...}
  # Y_i conditional on Age_i ~ Poi(theta)
    # so model log(E[Y_i | Age_i]) = b1 + b2*Age_i + b3*Age_i^2


#-----
# Data set
#-----
head(sparrows)
# 52 * 2
#     fledged age
#[1,]       3   3
#[2,]       1   3
#[3,]       1   1
#[4,]       2   1

## non-bayesian method 
sparrow_dat <- as.data.frame(sparrows)

fit.mle <- glm(fledged ~ age + age^2, data = sparrow_dat, family = "poisson")
fit.mle$coefficients
summary(fit.mle)

n <- length(y)
p <- ncol(X)

SSE <- sum(fit.mle$residuals^2) # 28.8
sigma2_ols <- 1/(n - p - 1) * SSE  # [1] 0.5997504


# variance of Beta use ols variance of Beta:
  # sigma2 * (t(X)X)_{inv}
  # sigma2 use sigma2_ols plug in
  

#----------
# prepare for Bayesian
#----------

y <- sparrows[, 1]

# design matrix
X <- cbind(rep(1, length(y)), sparrows[, -1], sparrows[, -1]^2)
X  #52 * 3
t(X) %*% X

# data set for Bayesian
yX <- cbind(y, X)
colnames(yX) <- c("fledged", "intercept", "age", "age2")
head(yX)

n <- length(y)
p <- ncol(X)



#----------------------
# Model for Metropolis
#----------------------

# the offspring for each bird i Yi 
  # Yi|x_i ~ Poisson (exp^{b*x_i})

# prior:
  # beta = (b1, b2, b3) ~iid N(0, 100)

# posterior :
  # beta | y, X has no closed form 

## resort to metropolis MC algo


#--------------------------------------
# Metropolis sampler for Poisson model
#--------------------------------------
library(mvtnorm)


# sigma2_hat = var(log(y + 1/2)) sample variance on log scale
  # + 1/2 to avoid log 0 = inft, as count y has 0's
Var_beta_ols <- var(log(y + 1/2)) * solve(t(X) %*% X) 
beta_0 <- 0  # starting beta element, will be replicate to a vector of 3

prior.m.b <- rep(0, ncol(X)) # mean of prior of beta
prior.sd.b <- rep(10, ncol(X)) # sd of prior of beta


Metro_Poisson <- function(beta_0, n.sim, Data = y, DM = X, 
                          burn_in = 1000, Thinning = T, thin = 10) {
  Beta <- rep(beta_0, ncol(X))
  #str(Beta) # num [1:3]
  
  BETA <- NULL
  counts <- 0
  
  for (i in 1:n.sim) {
    # proposal:
    beta_star <- rmvnorm(1, mean = Beta, 
                         sigma = Var_beta_ols) # [1, 1:3]
    beta_star <- t(beta_star)
    
    #L <- X %*% t(beta_star) # 52*1
    #sum(dpois(y, lambda = exp(L), log = T))
    
    # calculate AC ratio
    logr <- sum(dpois(y, exp(X %*% beta_star), log = T)) - 
      sum(dpois(y, exp(X %*% Beta), log = T)) + 
      sum(dnorm(beta_star, prior.m.b, prior.sd.b, log = T)) - 
      sum(dnorm(Beta, prior.m.b, prior.sd.b, log = T))
    
    
    # AC or RJ
    u <- runif(1)
    if (log(u) < logr) {
      Beta <- beta_star  # 3*1
      counts <- counts + 1
    }
    
    # update 
    BETA <- rbind(BETA, t(Beta)) # 1*3
  }
  
  # AC ratio
  AC_ratio <- counts / n.sim
  
  if (Thinning) {
    return(list(beta_pst = BETA[seq(burn_in, n.sim, by = thin), ], AC_ratio = AC_ratio))
  } else{
    return(list(beta_pst = BETA[burn_in:n.sim, ], AC_ratio = AC_ratio))
  }
}



#-----------
# Test on data
#-----------
Res <- Metro_Poisson(beta_0 = 0, n.sim = 500, Data = y, DM = X)

Res$AC_ratio # 0.456
# so the proposal is acceptable
# now re-run with longer chain for stationary

Res <- Metro_Poisson(beta_0 = 0, n.sim = 1e4, Data = y, DM = X, burn_in = 1000)
Res$AC_ratio

beta_pst <- Res$beta_pst
str(beta_pst)

#---------
# Plot
#--------
par(mfrow = c(1, 3))
for(idx in 1:3) {
  plot(1000:1e4, beta_pst[, idx], type = "l",
       main = bquote("Trace plot of " ~ beta[.(idx)]))
  
}


par(mfrow = c(1, 3))
for (idx in 1:3) {
  acf(beta_pst[, idx])$acf[2]
}



lag1_acf <- acf(beta_pst[, 2])$acf[2]

Res2 <- Metro_Poisson(beta_0 = 0, n.sim = 1e4, Data = y, DM = X, burn_in = 1000, 
                      Thinning = T, thin = 10)
Res2$AC_ratio # [1] 0.4272

beta_pst2 <- Res2$beta_pst

par(mfrow = c(1, 3))
for (idx in 1:3) {
  acf(beta_pst2[, idx])$acf[2]
}

sapply(1:3, function(x) acf(beta_pst2[, x], plot = F)$acf[2])
# [1] 0.1817538 0.2270815  0.2549073


library(coda)
effectiveSize(as.mcmc(beta_pst2[, x]))

sapply(1:3, function(x) effectiveSize(as.mcmc(beta_pst2[, x])))
# from the remaining 1000 samples after 1 tenth thinning
  # the effective sample sizes for each beta are:
#    var1     var1     var1 
# 609.5962 391.2855 356.7926 

## after thinning, we have nearly 1000 independent samples


## posterior marginal densities
par(mfrow = c(1, 3))
for(idx in 1:3) {
  plot(density(beta_pst2[, idx]), col = "red")
}


plot(density(beta_pst2[, 1]), col = "red")




#------------
# AGAIN
#------------

load("sparrows.RData")
head(sparrows)
#      fledged age
#[1,]       3   3
#[2,]       1   3
#[3,]       1   1

str(sparrows)
# num [1:52, 1:2]
# each bird give birth only once, either in age 1 or 2 or...

y <- sparrows[, 1]
Age <- sparrows[, -1]

# design matrix X
X <- cbind(rep(1, length(y)), sparrows[, -1], sparrows[, -1]^2)
head(X)

yX <- cbind(y, X)

colnames(yX) <- c("ofsp", "intercept", "age", "age2")
yX_dat <- as.data.frame(yX)
head(yX)
#     ofsp intercept age age2
#[1,]    3         1   3    9
#[2,]    1         1   3    9


## non-bayesian glm fitting for beta
mod <- glm(y ~ age + age2, data = yX_dat, family = "poisson")
summary(mod)

# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)  
#(Intercept)  0.27662    0.44219   0.626   0.5316  
#age          0.68174    0.33850   2.014   0.0440 *
#age2        -0.13451    0.05786  -2.325   0.0201 *
  ---

# can be used for comparison with Bayesian method result


  
#-----------------
# Bayesian version
#-----------------

# likelihood: 
  # yi | xi ~ Poisson (theta_x) 
      # theta_x is the true # of off-sp
      # theta_x = E[yi|xi]

   # yi is non-neg, so use log link
    # log(yi) = b1 + b2*age + b3*age2
    # yi = exp(beta^t X)
          # beta = (b1, b2, b3),
          # X = (1, age, age2)

# beta is unknown and of interest
# prior:
  # beta ~ iid N(b.prior.m, b.prior.sd)

  # b.prior.m = 0
  # b.prior.sd = 10


# posterior: 
  # beta | y, X under N likelihood and Poi prior 
    # no closed form 
    # resort to Metropolis


#--------------------------------
# Metropolis sampler for Poisson
#--------------------------------

# proposal:
    # Beta_star ~ MVN(Beta^(s), Sigma)

      # Beta^(s) = Beta_0
      # Sigma approximated by Var[beta_ols] = sigma2 * t(X)%*%X
          # sigma2 here is approximated var[log(y + 1/2)] 
          # while usually sigma2 is approximated by 1/(n-p-1) * SSE
          # but here beta is modeled when y is on log scale


# calculate AC ratio r:
  # logr <- sum(dpoi(y, exp(X %*% Beta_star), log = T)) - 
            # sum(dpoi(y, exp(X %*% Beta), log = T)) + 
            # sum(dnorm(Beta_star, rep(b.prior.m, 3), rep(b.prior.sd, 3), log = T)) - 
            # sum(dnorm(Beta, rep(b.prior.m, 3), rep(b.prior.sd, 3), log = T))



# AC or RJ
  # u ~ runif(1)
  # if (log(u) < logr) {
      # Beta <- Beta_star
      # Cnts <- Cnts + 1  # for AC ratio calculation 20% - 50%
  #}


library(mvtnorm)

Beta_0 <- 0
Metropolis_Poi <- function(Beta_0, b.prior.m = 0, b.prior.sd = 10, n.sim, Data = y, DM = X, burn_in) {
  #request(mvtnorm)
  Beta_0 <- rep(Beta_0, ncol(X))
  
  Beta <- Beta_0 # 1:3
  Sig_prop <- var(log(y + 1/2)) * solve(t(X) %*% X)
  
  #b.prior.m <- 0
  #b.prior.sd <- 10
  Cnts <- 0
  BETA <- NULL
  for (i in 1:n.sim) {
    # proposal:
    Beta_star <- rmvnorm(1, mean = Beta, sigma = Sig_prop) # [1, 1:3]
    Beta_star <- t(Beta_star) # 3*1
    
    
    # AC ratio 
    logr <- sum(dpois(y, exp(X %*% Beta_star), log = T)) -
      sum(dpois(y, exp(X %*% Beta), log = T)) +
      sum(dnorm(Beta_star, mean = rep(b.prior.m, ncol(X)), 
                sd = rep(b.prior.sd, ncol(X)), log = T)) -
      sum(dnorm(Beta, mean = rep(b.prior.m, ncol(X)), 
                sd = rep(b.prior.sd, ncol(X)), log = T))
    
    
    # AC or RJ
    u <- runif(1)
    if (log(u) < logr) {
      Beta <- Beta_star
      Cnts <- Cnts + 1
    }
    
    # update desired Beta
    BETA <- rbind(BETA, t(Beta)) # 1:3
  }
  
  # calculate AC ratio for proposal tunning, 20%-50% AC ratio
  AC_ratio <- Cnts / n.sim
  
  return(list(Beta_pst = BETA[(burn_in+1):n.sim, ], AC_ratio = AC_ratio))
}

#----------
try <- rbind(t(y[1:3]), 1)
try[1:3] # as a vector element selection
try[1:2, ] # select first two rows
#----------



Res <- Metropolis_Poi(Beta_0 = 0, b.prior.m = 0, b.prior.sd = 10, n.sim = 500, 
               Data = y, DM = X, burn_in = 0)


Res$AC_ratio # [1] 0.426







#--------
# Data set（depreached）
#--------
set.seed(1)

## try make sure age 2 has the largest, age 6 smallest offsp
springs <- round(rnorm(6, 3, sqrt(2)), 0)
Offsp <- vector("numeric", length = 6)

Offsp[2] <- max(springs)
Offsp[6] <- min(springs)

idx1 <- which(springs == max(springs))[1]
idx2 <- which(springs == min(springs))[1]

Offsp[-c(2, 6)] <- springs[-c(idx1, idx2)]
#---------------------------------------


OFFSP <- NULL
for (i in 1:52) {
  springs <- round(rnorm(6, 3, sqrt(2)), 0)
  
  Offsp <- vector("numeric", length = 6)
  
  Offsp[2] <- max(springs)
  Offsp[6] <- min(springs)
  idx1 <- which(springs == max(springs))[1]
  idx2 <- which(springs == min(springs))[1]
  Offsp[-c(2, 6)] <- springs[-c(idx1, idx2)]
  
  OFFSP <- rbind(OFFSP, abs(Offsp))
}

OFFSP # 20 * 6
y <- c(t(OFFSP))
Age <- replicate(52, seq(1:6))
Age <- c(Age)
ID <- rep(1:52, each = 6)

Offspring_dat <- data.frame(ID = ID, Age = Age, Ofsp = y)
head(Offspring_dat, 15)


boxplot(Offspring_dat$Ofsp ~ Offspring_dat$Age)
# age 2 has the highest offspring number
# age 6 has the least


