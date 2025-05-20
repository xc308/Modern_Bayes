#==========================================
# M-H sampler for Tumour Location Detection
#==========================================

Aim:
     - Construct a Bayesian generalised linear mixed effects model to detect tumour location
     - Sampling methods: hybrid Metropolis - Gibbs


load("tumorLocation.RData")

Y <- tumorLocation
head(Y)
str(Y) # num [1:21, 1:20]

# create location
loc <- seq(1, 20, 1)
xs <- loc / 20 


#----------------------
# Exploratory plots
#----------------------

par(mfrow = c(1, 2))

## Each mouse curve at each location and average (over 21 mice)
plot(c(0, 1), range(Y), xlab = "Location",
     ylab = "Tumor counts", type = "n")

for (i in 1:nrow(Y)) {
  lines(xs, Y[i, ], col = "gray")
}

Avg_y <- apply(Y, 2, mean)
lines(xs, Avg_y, lwd = 3)


##--------
# To find the representative polynomial for log(Y_xj)
##--------
# Y_xj : TM counts at location x for jth mouse
# log(Y_xj): to ensure counts must be non-neg

# to model log(Y_xj)
  # use the relationship between log(Avg_Y_x) and location x
  # as representative

log_avg_y <- log(apply(Y, 2, mean))
plot(xs, log_avg_y, type = "l")

# use polynomial to fit to this log avg y curve
# 1st construct design matrix X

X <- cbind(rep(1, length(xs)), poly(xs, 4, raw = T))
# evaluate raw xs polynomials rather than orthogonal to constant

#try <- poly(xs, 4)
#head(try)
#try2 <- poly(xs, 4, raw = T)
#head(try2)

head(X)
#          1      2        3          4
# [1,] 1 0.05 0.0025 0.000125 0.00000625
# [2,] 1 0.10 0.0100 0.001000 0.00010000
# [3,] 1 0.15 0.0225 0.003375 0.00050625


# fit linear regression polynomial to the log avg y curve

fit2 <- lm(log_avg_y ~ -1 + X[, 1:3])$coef # quadratic
fit3 <- lm(log_avg_y ~ -1 + X[, 1:4])$coef # cubic
fit4 <- lm(log_avg_y ~ -1 + X[, 1:5])$coef # quatic

# get the corresponding y values under the fitted coef
y2 <- X[, 1:3] %*% fit2 
y3 <- X[, 1:4] %*% fit3
y4 <- X[, 1:5] %*% fit4

plot(c(0, 1), range(c(log_avg_y, y2, y3, y4)),
     xlab = "Location", ylab = "log avg TM counts",
     type = "n")

lines(xs, log_avg_y, lwd = 2)

lines(xs, y2, col = "gray")
points(xs, y2, pch = "2", col = "gray")

lines(xs, y3, col = "blue")
points(xs, y3, pch = "3", col = "blue")

lines(xs, y4, col = "red")
points(xs, y4, pch = "4", col = "red")

## Conclusion:
  # 4th degrees polynomial fits the log avg y counts best
  # so use poly(4) as the model for log(Y_xj)
  # Y_xj ~ iid Poisson(exp(Y_xj))


#=======================================
# M-H sampler for desired theta, Sigma, beta_j, j = 1:m (m = 21)
#=======================================

#----------------------
# parameters for prior
#----------------------
BETA <- NULL
for (i in 1:nrow(Y)) {
  betaj <- lm(log(Y[i, ]) ~ -1 + X)$coef
  BETA <- rbind(BETA, betaj)
}

# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# NA/NaN/Inf in 'y'

# So, log(Y[i, ] + 1/20) to avoid Y[i, ] very close to 0

range(Y) # [1]  0 17 Y has 0 elements 

BETA <- NULL
for (i in 1:nrow(Y)) {
  betaj <- lm(log(Y[i, ] + 1/20) ~ -1 + X)$coef
  BETA <- rbind(BETA, betaj)
}

str(BETA )
# num [1:21, 1:5] 


mu_0 <- apply(BETA, 2, mean)
L_0 <- S_0 <- cov(BETA)

p <- ncol(X)
eta_0 <- p + 2


#-------
# starting values
#-------

Sigma <- cov(BETA)
#Sigma_0 <- cov(BETA)

#------------------------------------------
# logmvnorm function for Metropolis sampler
#------------------------------------------

logmvnorm <- function(Beta, theta, Sigma) {
  dSigma <- det(Sigma)
  iSigma <- solve(Sigma)
  
  -0.5 * (length(Beta) * log(2 * pi) + log(dSigma)) -
    0.5 * t(Beta - theta) %*% iSigma %*% (Beta - theta)
}


# verify: correct
logmvnorm(Beta = BETA[1, ], theta = mu_0, Sigma = Sigma_0) # -24.35743
dmvnorm(BETA[1, ], mu_0, Sigma_0, log = T) # [1] -24.35743



#------------
# M-H sampler
#------------
library(mvtnorm)
library(MCMCpack) # for riwish(v, S)

n.sim <- 1e4
iL_0 <- solve(L_0)
m <- nrow(Y)
cnts <- 0

THETA_pst <- NULL
SIGMA_pst <- NULL

for (i in 1:n.sim) {
  
  # sample theta ~ MVN(theta_m, L_m)
  iSigma <- solve(Sigma)
  beta_bar <- apply(BETA, 2, mean)
  
  L_m = solve(iL_0 + m * iSigma)
  mu_m = L_m %*% (iL_0 %*% mu_0 + m * iSigma %*% beta_bar)
  
  theta <- rmvnorm(1, mean = mu_m, sigma = L_m) 
  theta <- t(theta) # 5 * 1
  
  
  # sample Sigma ~ InvWishart(eta_m, S_m)
  eta_m <- eta_0 + m
  S_theta <- (t(BETA) - c(theta)) %*% t(t(BETA) - c(theta))
  S_m <- S_0 + S_theta
  
  Sigma <- riwish(v = eta_m, S = S_m)
  
  
  # Metropolis sample BETA for each of mouse j
  for (j in 1:m) {
    # 1. propose betaj_star
    Betaj_start <- rmvnorm(1, mean = BETA[j, ], sigma = 0.5*Sigma)
    Betaj_start <- t(Betaj_start)
    
    # 2. calculate log(r)
    logr <- sum( dpois(Y[j, ], lambda = exp(X %*% Betaj_start), log = T) -
                   dpois(Y[j, ], lambda = exp(X %*% BETA[j, ]), log = T)) +
      dmvnorm(c(Betaj_start), c(theta), Sigma, log = T) -
      dmvnorm(BETA[j, ], c(theta), Sigma, log = T)
      # dmvnorm requires vector input
    
    
    # 3. AC or RJ
    if (log(runif(1)) < logr) {
      BETA[j, ] <- Betaj_start
      cnts <- cnts + 1
    }
  }
  
  # Print AC ratio 20% - 50%
  cat(i, cnts / (m * i), "\n")
  
  # Update and collect posterior results every 10 samples 
  if (i %% 10 == 0) {
    THETA_pst <- rbind(THETA_pst, t(theta))
    SIGMA_pst <- rbind(SIGMA_pst, c(Sigma))
  }
}


#---------
# Results
#---------
# AC ratio is 0.3045857 
THETA_pst
SIGMA_pst


#-------------------------------------------
# Quantiles of exp(thetaX), exp(betaX), Y|x
#-------------------------------------------

str(THETA_pst)
# num [1:1000, 1:5]
str(X) # num [1:20, 1:5]

## quantile of exp(thetaX)
thetaX <- THETA_pst %*% t(X) 
#str(thetaX) # num [1:1000, 1:20]
exp_thetaX <- exp(thetaX)
exp_thetaX_CI <- apply(exp_thetaX, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
t(exp_thetaX_CI)

## quantile of exp(betaX)
# 1. get posterior BETA

str(SIGMA_pst) # num [1:1000, 1:25]

BETA_pst <- NULL
for(i in 1:1000) {
  BETAj_pst <- rmvnorm(1, mean = THETA_pst[j, ], 
                       sigma = matrix(SIGMA_pst[j, ], 5, 5))
  BETA_pst <- rbind(BETA_pst, BETAj_pst)
}

str(BETA_pst)
# num [1:1000, 1:5]
BETAX <- BETA_pst %*% t(X)
exp_BETAX <- exp(BETAX)
str(exp_BETAX) # num [1:1000, 1:20]

BETAX_CI <- apply(exp_BETAX, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
t(BETAX_CI )


# predictive distribution of Y|x
# log(y_x.) = BETAX
# Y_xj ~ iid Poisson(exp(BETAX))
      # str(exp_BETAX) # num [1:1000, 1:20]
Y <- rpois(1, lambda = exp_BETAX)

Y <- NULL
for (j in 1:1000) {
  Yj <- rpois(dim(exp_BETAX)[2], lambda = exp_BETAX[j, ])
  Y <- rbind(Y, Yj)
}
str(Y)
# int [1:1000, 1:20]

Y_CI <- apply(Y, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
t(Y_CI)


#----------------------------------
# Plot above 3 quantities quatiles 
#----------------------------------

par(mfrow = c(1, 3))

plot(c(0, 1), range(c(exp_thetaX_CI[1, ], exp_thetaX_CI[2, ], exp_thetaX_CI[3, ])),
     xlab = "Location", 
     ylab = expression(exp(theta*X)),
     type = "n")

lines(xs, exp_thetaX_CI[1,])
lines(xs, exp_thetaX_CI[2,], lwd = 2)
lines(xs, exp_thetaX_CI[3,])


plot(c(0, 1), range(c(BETAX_CI[1, ], BETAX_CI[2, ], BETAX_CI[3,])),
     xlab = "Location", 
     ylab = expression(exp(beta*X)),
     type = "n")

lines(xs, BETAX_CI[1, ])
lines(xs, BETAX_CI[2, ], lwd = 2)
lines(xs, BETAX_CI[3, ])


plot(c(0, 1), range(c(Y_CI[1, ], Y_CI[2, ], Y_CI[3,])),
     xlab = "Location", 
     ylab = "Y|x",
     type = "n")

lines(xs, Y_CI[1, ])
lines(xs, Y_CI[2, ], lwd = 2)
lines(xs, Y_CI[3, ])

# Conclusion
  # the 1st panel shows the variance in the fix 
      # but unknown theta
  
  # 2nd panel shows the variance and covariance in the 
      # betaj in addition to the above theta
      # to reflect the cross-mouse heterogeneity

  # 3rd panel shows the fluction of actual observed
      # TM counts around it's true expected TM counts
      # in addition to the above two sources of variance

  # if we want to predict the TM counts of a new mouse
    # need to use the 3rd CI band
  
  # if just want to use assess the uncertainty of theta
    # then just need to use the 1st CI band. 












