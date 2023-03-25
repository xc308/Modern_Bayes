#==================================
# Hoff Chp11 GL Mix effects Model
#==================================
# Analysis of tumor location data

load("tumorLocation.RData")
head(tumorLocation)

Y <- tumorLocation
m <- nrow(Y)

Loc <- seq(1, 20, by = 1)
xs <- Loc / 20


plot(c(0, 1), range(Y), xlab = "location", ylab = "number of tumors",
     type = "n")
for (j in 1:m) {
  lines(xs, Y[j, ], col = "gray", type = "l")
  # joint the correponding points
}

# average of 21 mice at each location
AVG_Loci <- apply(Y, 2, mean)
lines(xs, AVG_Loci, lwd = 3)



#------------------------
# To fit the log average
#------------------------
# log average tumor counts is the lambda rate in the Poission
# our model object
# want to model the log avg TM counts = polynomial of certain degrees

  # from the above plot, we know the shape of the density of the log avg
  # try polynomial of different degrees 



### 1st, Log average of TM counts for each location
log_Y_avg <- log(apply(Y, 2, mean))


# 2nd, construct design matrix (DM)
X <- cbind(rep(1, ncol(Y)), poly(xs, 4, raw = T))
  # up to 4 degrees of freedom of xs, ie xs, xs^2, xs^3, xs^4
  # raw = T means the xs, xs^2, xs^3, xs^4 are not orthogonal to intercept
head(X, 3)


# 3rd, fit a lm using above log_Y_avg and X
fit2 <- lm(log_Y_avg ~ -1 + X[, 1:3]) # fit quadratic
fit3 <- lm(log_Y_avg ~ -1 + X[, 1:4]) # fit cubic
fit4 <- lm(log_Y_avg ~ -1 + X[, 1:5]) # fit quatic


## To see which line fit the density of log avg of TM counts best
  # need to get the density of each of these curve at each location

y2 <- X[, 1:3] %*% fit2$coefficients
y3 <- X[, 1:4] %*% fit3$coefficients
y4 <- X[, 1:5] %*% fit4$coefficients

plot(c(0, 1), range(c(log_Y_avg, y2, y3, y4)), 
     xlab = "location", ylab = "log avg tumor counts",
     type = "n")
lines(xs, log_Y_avg, lwd = 3, type = "l")

lines(xs, y2, col = "gray")
points(xs, y2, pch = "2", col = "gray")

lines(xs, y3, col = "blue")
points(xs, y3, pch = "3", col = "blue")

lines(xs, y4, col = "red", lwd = 2)
points(xs, y4, pch = "4", col = "red")

## the polynomial of degree 4 is the best fit to the density
  # of log avg of TM counts


X <- cbind(rep(1, ncol(Y)), poly(xs, 4, raw = T))
# raw = T, each of polynomials are not orthogonal to constant 1 col

fit2 <- lm(log_Y_avg ~ -1 + X[, 1:3]) # quadratic
fit3 <- lm(log_Y_avg ~ -1 + X[, 1:4]) # cubic
fit4 <- lm(log_Y_avg ~ -1 + X[, 1:5]) # quatic


y2 <- X[, 1:3] %*% fit2$coefficients
y3 <- X[, 1:4] %*% fit3$coefficients
y4 <- X[, 1:5] %*% fit4$coefficients

# X[, 1:3]: 20*3, 
# fit2$coefficients : 3*1

# X[, 1:4]: 20*4
# fit3$coefficients: 4*1

# str(fit4$coefficients) # Named num [1:5]

plot(xs, AVG_Loci,
     ylim = range(c(AVG_Loci, y2, y3, y4)),
     type = "n", xlab = "locations", 
     ylab = "log avg number of tumors",
     lwd = 3)

lines(xs, AVG_Loci, lwd = 3)
lines(xs, y2)
points(xs, y2, pch = "2", col = "black")

lines(xs, y3, col = "gray")
points(xs, y3, pch = "3", col = "red")

lines(xs, y4, col = "blue")
points(xs, y4, pch = "4", col = "blue")

## 
# poly with degrees of 4 is the best fit to the 
  # log avg of number of tumors


#---------------------
# Paramters on priors
#---------------------

# get the OLS estimates for beta
str(Y) #  num [1:21, 1:20]
m <- nrow(Y)
str(X) # num [1:20, 1:5]
p <- ncol(X)


BETA <- NULL
for (j in 1:m) {
  BETA <- rbind(BETA, lm(log(Y[j, ] + 1/20) ~ -1 + X)$coef)
}
str(BETA)
# num [1:21, 1:5]]


mu_0 <- apply(BETA, 2, mean)
# 1*5

L_0 <- S_0 <- cov(BETA)
# 5*5

eta_0 <- p + 2


#----
# starting values
#-----

theta_0 <- apply(BETA, 2, mean)
Sigma_0 <- cov(BETA)
Beta_0 <- BETA 



#-------------
# M-H samplers
#-------------
library(mvtnorm)
library(MCMCpack) # for Invwishart riwish(v, S)


Mix_GLM_MH_sampler <- function(theta_0, Sigma_0, BETA_0, Data = Y, DM = X, n.sim) {
  
  m <- nrow(Y)
  theta <- theta_0
  Sigma <- Sigma_0
  Beta <- BETA_0  # 21*5
  
  THETA <- NULL
  SIGMA <- NULL
  
  cnts <- 0
  iL_0 <- solve(L_0)
  
  for (i in 1:n.sim) {
    
    # update theta ~ MVN(mu_m, Sigma_m)
    iSigma <- solve(Sigma)
    beta_bar <- apply(Beta, 2, mean)
    
    L_m <- solve(iL_0 + m * iSigma)
    mu_m <- L_m %*% (iL_0 %*% mu_0 + m * iSigma %*% beta_bar)
    
    theta <- rmvnorm(1, mean = mu_m, sigma = L_m) # [1, 1:5]
    theta <- t(theta)
    
    
    # update Sigma ~ InvWishart(eta_m, S_m)
    eta_m = eta_0 + m
    S_theta <- (t(Beta) - c(theta)) %*% t(t(Beta) - c(theta))
    S_m = S_0 + S_theta
    
    Sigma <- riwish(v = eta_m, S = S_m) #1:5, 1:5

    
    ## update Beta for each j 1:m, m = 21
    iSigma <- solve(Sigma); dSigma <- det(Sigma)
    for (j in 1:m) {
      # 1. propose
      betaj_star <- rmvnorm(1, mean = Beta[j, ], 0.5 * Sigma)
      betaj_star <- t(betaj_star)
      
      # 2. logr
      #logr <- sum(dpois(Y[j, ], lambda = exp(X %*% betaj_star), log = T) -
      #  dpois(Y[j, ], lambda = exp(X %*% Beta[j, ]), log = T)) + 
      #  logmvnorm(betaj_star, theta, Sigma) - 
      #  logmvnorm(Beta[j, ], theta, Sigma)
      
    
      #logr <- sum(dpois(Y[j, ], lambda = exp(X %*% betaj_star), log = T) -
                   # dpois(Y[j, ], lambda = exp(X %*% Beta[j, ]), log = T)) + 
        #dmvnorm(c(betaj_star), c(theta), Sigma, log = T) - 
        #dmvnorm(Beta[j, ], c(theta), Sigma, log = T)
      
       logr <- sum(dpois(Y[j, ], lambda = exp(X %*% betaj_star), log = T) -
                    dpois(Y[j, ], lambda = exp(X %*% Beta[j, ]), log = T)) + 
        ldmvnorm(t(betaj_star), theta, Sigma, iSigma=iSigma, dSigma=dSigma) -
        ldmvnorm( t(Beta[j,]),theta,Sigma, iSigma=iSigma, dSigma=dSigma )
       
       
      # 3. AC or RJ
      if (log(runif(1)) < logr) {
        Beta[j, ] <- betaj_star
        cnts <- cnts + 1
      }
    }
    
    
    # AC ratio 
    AC_ratio = cnts / (m * i)
    cat(i, AC_ratio, "\n")
    
    
    # collect the updates for this run if i %% 10 == 0
    if (i %% 10 == 0) {
      THETA <- rbind(THETA, t(theta))
      SIGMA <- rbind(SIGMA, c(Sigma))
    }
    
  }
  return(list(Theta_pst = THETA, Sigma_pst = SIGMA))
}



Res_100 <- Mix_GLM_MH_sampler(theta_0 = theta_0, Sigma_0 = Sigma_0, 
                              BETA_0 = BETA, Data = Y, DM = X, n.sim = 1000)



theta_pst <- Res_100$Theta_pst # 100 * 5
Sigma_pst <- Res_100$Sigma_pst
Beta_pst <- Res_100$Beta_pst # 2100 *5
head(Beta_pst)


plot(1:2100, Beta_pst[, 1], type = "l")
plot(1:2100, Beta_pst[, 2], type = "l")

plot(1:100, theta_pst[, 1], type = "l")

library(coda)
effectiveSize(theta_pst)
#     var1      var2      var3       var4      var5 
#   100.00000 100.00000 100.00000  100.00000  72.74689 
# almost independt MCMC samples for all 5 theta components

effectiveSize(Sigma_pst)
# almost all 25 element of the 5*5 Sigma are indepents MCMC sampels

effectiveSize(Beta_pst)


Res_5e4 <- Mix_GLM_MH_sampler(theta_0 = theta_0, Sigma_0 = Sigma_0, 
                              BETA_0 = BETA, Data = Y, DM = X, n.sim = 5e4)

theta_pst <- Res_5e4$Theta_pst
str(theta_pst) # num [1:5000, 1:5]

Sigma_pst <- Res_5e4$Sigma_pst
str(Sigma_pst) # num [1:5000, 1:25]

Beta_pst <- Res_5e4$Beta_pst
str(Beta_pst)
# num [1:105000, 1:5]


effectiveSize(theta_pst)
# var1 var2 var3 var4 var5 
#5000 5000 5000 5000 5000 
# Monte Carlo standard error is approximately the same 
# as posterior of standard error


effectiveSize(Sigma_pst)
# the majority sample size of elements of Sigma_pst are 5000
# few elements have negative sample correlations ending 
# larger sample size than 5000 (at 5200)

effectiveSize(Beta_pst)
#         X         X1         X2        X3         X4 
#   20175.835 117790.678  12054.060   3872.510   2162.569 
# beta_3 beta_4 are positively correlated with effN around 3000-4000
# beta_1, 2, 3 are heavily negatively correlated
# 



#------------
# Quantiles fo exp(theta x), exp(beta x), Y|x
#------------

theta_pst <- Res_5e4$Theta_pst
str(theta_pst) # num [1:5000, 1:5]

str(X) # num [1:20, 1:5] 


exp_thetaX <- exp(theta_pst %*% t(X)) # 5000 * 20
exp_thetaX_CI <- apply(exp_thetaX, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
t(exp_thetaX_CI)
#             2.5%       50%      97.5%
#[1,] 0.008255019 0.1726579  4.3349878
#[2,] 0.036299381 0.1756825  0.8462361
#[3,] 0.052343186 0.1806182  0.6019891

str(xs) #num [1:20]


par(mfrow = c(1, 3))
#plot(xs, t(exp_thetaX_CI)[, 1], type = "l")

plot(c(0, 1), range(c((t(exp_thetaX_CI)[, 1]), 
                (t(exp_thetaX_CI)[, 2]),
                (t(exp_thetaX_CI)[, 3]))),
     xlab = "location", 
     ylab = "number of turmors", 
     type = "n")

for (j in 1:3) {
  lines(xs, t(exp_thetaX_CI)[, j])
}



#---
# exp(X beta)
#---

# 1st get beta from the theta_pst and Sigma_pst


theta_pst <- Res_5e4$Theta_pst
str(theta_pst) # num [1:5000, 1:5]

Sigma_pst <- Res_5e4$Sigma_pst
str(Sigma_pst) # num [1:5000, 1:25]

EXP_BETAX <- NULL
for (i in 1:5000) {
  beta <- rmvnorm(1, mean = theta_pst[i, ],
                  sigma = matrix(Sigma_pst[i, ], 5, 5))
  # 5000 * 5
  
  EXP_BETAX <- rbind(EXP_BETAX, exp(beta %*% t(X))) 
  # 5000 * 20
}

str(EXP_BETAX)
# num [1:5000, 1:20]


EXP_BETAX_CI <- apply(EXP_BETAX, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
t(EXP_BETAX_CI)
str(EXP_BETAX_CI) # num [1:3, 1:20]
range(EXP_BETAX_CI[1, ])



plot(c(0, 1), range(c(EXP_BETAX_CI[1, ],
                      EXP_BETAX_CI[2, ])),
     xlab = 'location', ylab = "", type = "n")

for(i in 1:2) {
  lines(xs, EXP_BETAX_CI[i, ])
}


yXT.pp<-matrix( rpois(prod(dim(eXB.post)),eXB.post),
                dim(eXB.post)[1],dim(eXB.post)[2] )


#-----
# y ~ Poisson(exp(xb))
#-----

str(EXP_BETAX)
# num [1:5000, 1:20]

Y_pst <- matrix(rpois(n = prod(dim(EXP_BETAX)),lambda = EXP_BETAX),
       dim(EXP_BETAX)[1], dim(EXP_BETAX)[2])

str(Y_pst)
# num [1:5000, 1:20] 

Y_pst_CI <- apply(Y_pst, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
t(Y_pst_CI)

