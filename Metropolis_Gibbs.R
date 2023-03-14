#=======================================
# Hoff 10.5 Combing Metropolis and Gibbs
#=======================================

load("icecore.RData")

head(icecore)
#     year   co2  tmp
#1 -419095 277.6 0.94
#2 -416872 286.9 1.36

# data include 200 values of temperature measured at roughly equal time intervals
# time between consecutive measurements being approximately 2,000 years
# For each value of temperature there is a CO2 concentration value 
    # corresponding to a date that
    # is roughly 1,000 years previous to the temperature value, on average.
# Temperature is recorded in terms of its difference from 
    # current present temperature in degrees Celsius,

# CO2 concentration is recorded in parts per million by volume.


max((icecore[, 3] - mean(icecore[, 3])) / sd(icecore[, 3])) # [1] 2.867939
min((icecore[, 3] - mean(icecore[, 3])) / sd(icecore[, 3])) # 
# [1] -1.664384

# less variation than on original scale
min(icecore[, 3])
# [1] -9.39

Mx <- c()
Mn <- c()
for (i in 2:3) {
  Mx[i-1]<- max((icecore[, i] - mean(icecore[, i])) / sd(icecore[, i]))
  Mn[i-1] <- min((icecore[, i] - mean(icecore[, i])) / sd(icecore[, i]))
}

Mx # [1] 2.089509 2.867939
Mn # [1] -1.604304 -1.664384



# mgp
# The margin line (in mex units) for the axis title, 
# axis labels and axis line. Note that mgp[1] affects title 
# whereas mgp[2:3] affect axis. The default is c(3, 1, 0).

pdf("fig10_7.pdf",family="Times") 
# font change from "Helctive to Times"

par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
layout(matrix( c(1,1,2),nrow=1,ncol=3))
layout(matrix( c(1,1,2),nrow=1,ncol=3))
# a matrix object specifying the location of the next NN figures on the output device
# Each value in the matrix must be 0 or a positive integer.
# If N is the largest positive integer in the matrix,
# then the integers {1,…,N−1} must also appear at least once in the matrix.

par(mfrow = c(1, 2))

plot(icecore[,1], (icecore[, 3] - mean(icecore[, 3])) / sd(icecore[, 3]),
     type = "l", col = "black", ylim = c(-2.5, 3), 
     xlab = "year", ylab = "standardized measurement")

lines(icecore[, 1], (icecore[, 2] - mean(icecore[, 2])) / sd(icecore[, 2]),
      type = "l", col = "grey")

legend("topright", legend = c("temp", expression(CO[2])), 
       bty = "n", lwd = c(2, 2), col = c("black", "grey"))

plot(icecore[, 2], icecore[, 3], 
     xlab = expression(paste(CO[2], " (ppmv)")), 
     ylab = "temperature difference (deg C)")

## Remark:
  # plot indicates that the temporal history of 
    # temperature and CO2 follow very similar patterns.
  # 2nd plot indicates that CO2 concentration at a given 
    # time point is predictive of temperature 
      # following that time point.
  # One way to quantify this is by fitting a linear regression model 
    # for temperature (Y ) as a function of CO2 (x).

# Ordinary least squares regression
  # gives an estimated model of ˆE[Y |x] = −23.02+0.08x with a nominal standard
  # error of 0.0038 for the slope term.

# The validity of this standard error relies
    # on the error terms in the regression model being 
    # independent and identically distributed, and
    # standard confidence intervals further rely on the errors being
    # normally distributed.




## Residual correlation visualisation
lmfit <- lm(icecore$tmp ~ icecore$co2)
hist(lmfit$residuals)
acf(lmfit$residuals, xlab = "lag", ci.col = "gray")
acf(lmfit$residuals, plot = F)$acf[2]
# [1] 0.5165516


summary(lmfit)
# Coefficients:
#           Estimate Std. Error t value
#(Intercept) -23.024140   0.879543  -26.18
#icecore$co2   0.079853   0.003834   20.83

# the std error is nominal as it dependes on the assumption
  # the the error being iid, 
# and the CI dependends on the error being normal

# the 1st plot indicate no obvious deviation from normality
# the 2nd plot shows non-trivial auto-correlation 0.52 betw
  # two residuals at time consecutives

# such postive correlation generally means
  # there's less information in the data and 
  # less evidence for relationship btw two variables 
    # than is assumed by OLS regression analysis


S <- 1e4
thin <- c(1,(1:1000)*(S/1000))
head(thin) # [1]  1 10 20 30 40 50
tail(thin) # [1]  9950  9960  9970  9980  9990 10000

thin2 <- seq(1, 1e4, by = 10)
head(thin2) # [1]  1 11 21 31 41 51
tail(thin2)

thin3 <- seq(1, 1e4, by = 9)
head(thin3) # [1]  1 10 19 28 37 46


head(icecore)


#-------------------------------
# How to code correlation C_rho
#-------------------------------

# C_rho = matrix(c(1, rho, rho^2, rho^3, 
                #  rho, 1, rho^2, rho^3, 
                #  rho^2, rho, 1, rho^3, 
                #  rho^3, rho^2, rho, 1))

# 1st create rho
  # lag-1 acf
n <- 5
lmfit <- lm(icecore$tmp ~ icecore$co2)
summary(lmfit)
phi <- acf(lmfit$residuals, plot = F)$acf[2]

# 2nd create exponent matrix 
DY <- abs(outer( (1:n),(1:n) ,"-"))

# 3rd create the desired correlation matrix C_rho
phi^DY



#---------------------------
# How to code log(det(C_rho))
#---------------------------

determinant(phi^DY, log = T)
# $modulus
#[1] -1.241486
#attr(,"logarithm")
#[1] TRUE

#$sign
#[1] 1

#attr(,"class")
#[1] "det"

determinant(phi.p^DY,log=TRUE)$mod


#----------------------------
# Metropolis - Gibbs Sampler
#----------------------------

# likelihood of data;
  # Y|X, beta, sigma2, rho ~ MVN(X beta, sigma2 C_rho)

# likehood of data as a function of beta:
  # beta | y, X, sigma2, rho ~ MVN((t(X)X)^{-1} t(X)y, sigma2 C_rho t(X)X)^{-1})
  

# desired parameters beta, sigma2, rho
# prior of each of them:

  # beta ~ MVN(beta0, Sigma0)
  # sigma2 ~ InvGamma(a, b)
  # rho ~ unif(1)

# posterior:
  # beta | y, X, sigma2, rho ~ MVN(beta_n, Sigma_n)
      # Sigma_n <- solve(Sigma0_inv + t(X) %*% solve(C_rho) %*% X / sigma2)
      # beta_n <- Sigma_n %*% [Sigma0_inv %*% beta0 + t(X) %*% solve(C_rho) %*% y / sigma2]

  # sigma2 | y, X, beta, rho ~ InvGamma(An, Bn)
      # An <- a + n/2
      # Bn <- b + 1/2 SSR_rho
        # SSR_rho <- t(y - X beta) %*% solve(C_rho) %*% (y - X beta)


# rho: 
  # no closed form of posterior distribution
  # Metropolis: 
    # 1. propose rho_star ~ Unif(rho - delta, rho + delta)

    # 2. r = p(rho_star|y, X, beta, sigma2) / p(rho|y, X, beta, sigma2)
      # = p(y|X, beta, sigma2, rho_star) * p(rho_star) 
        #/p(y|X, beta, sigma, rho) * p(rho)

    # pdf density for p(y|X, beta, sigma2, rho_star) 
      # propto det(sigma2 * C_rho_star)^{-1/2} 
                  # exp^{1/2 t(y - X beta) (sigma2 C_rho_star) (y - X beta) }
    
    # pdf density for p(y|X, beta, sigma2, rho)
      # propto det(sigma2 * C_rho)^{-1/2}
                  # exp^{1/2 t(y - X beta) (sigma2 C_rho) (y - X beta)}

    # logr = -.5 * [log(det(sigma2 C_rho_star)) - log(det(sigma2 C_rho))]
        # -.5 * trace[(y - X beta) %*% t(y - X beta) %*% (solve(C_rho_star) - solve(C_rho))] / sigma2


  # 3. u ~ Unif(1)
      # log(u) < logr AC, rho <- rho_star 
    


#------------------------------------------
# data, design matrix, correlation matrix
#------------------------------------------

y <- icecore$tmp
X <- cbind(1, icecore$co2)

n <- length(y)
DY <- abs(outer(1:n, 1:n, "-"))
str(DY) # int [1:200, 1:200]


# to get rho need lag-1 acf of residual of ols fit
# C_rho <- rho^DY

lmfit <- lm(y ~ icecore$co2)
summary(lmfit)
rho <- acf(lmfit$residuals, plot = F)$acf[2]

C_rho <- rho^DY # num [1:200, 1:200]


#---------------
# starting value
#---------------

# beta: beta_ols
beta_ols <- lmfit$coefficients
str(beta_ols) # Named num [1:2] -23.0241 0.0799

# sigma2: sigma2_ols
sum_lmfit <- summary(lmfit)
sigma2_ols <- sum_lmfit$sigma^2  # [1] 2.349845

# rho: lag-1 acf of residual from ols fit
rho_ols <- acf(lmfit$residuals, plot = F)$acf[2]
# [1] 0.5165516


#---------------------
# Parameters in prior
#---------------------

beta0 <- rep(0, ncol(X))
#Sigma0 <- diag(1000, nrow = ncol(X))
Sigma0_inv <- diag(1/1000, nrow = ncol(X))

a = 1/2
b = 1/4

delta = 0.1


#----------------------------------------------------
# Metropolis-Gibbs sampler for error-correlated data
#----------------------------------------------------
library(mvtnorm)

Metro_Gibbs_er_corr <- function(beta0, sigma2_0, rho_0, Data = y, DM = X, n.sim, burn.in, 
         Thinning = T, thin) {
  
  n <- length(y)
  sigma2 <- sigma2_0
  rho <- rho_0
  
  Cnts <- 0
  OUTPUT <- NULL
  
  for (i in 1:n.sim) {
    
    # sample beta ~ MVN(beta_n, Sigma_n)
    C_rho <- rho^DY
    Sigma_n <- solve(Sigma0_inv + t(X) %*% solve(C_rho) %*% X / sigma2)
    beta_n <- Sigma_n %*% (Sigma0_inv %*% beta0 + t(X) %*% C_rho %*% y / sigma2)
    
    beta <- rmvnorm(1, mean = beta_n, sigma = Sigma_n)
    beta <- t(beta) # 2*1
    
    
    # sample sigma2 ~ InvGamma(An, Bn)
    An <- a + n/2
    SSR_r <- t(y - X %*% beta) %*% solve(C_rho) %*% (y - X %*% beta)
    Bn <- b + SSR_r / 2
    
    sigma2 <- 1/rgamma(1, shape = An, rate = Bn)
    
    
    # sample rho using Metropolis
    # propose: rho must be in (0, 1)
    rho_star <- abs(runif(1, rho - delta, rho + delta))
    rho_star <- min(rho_star, 2 - rho_star) 
    
  
    # calculate r or logr
    C_rho_star <- rho_star^DY
    logr <- -1/2 * (determinant(C_rho_star, log = T)$mod - 
                      determinant(C_rho, log = T)$mod) -
      1/2 * sum(diag((y - X %*% beta) %*% t(y - X %*% beta) %*% 
                       (solve(C_rho_star) - solve(C_rho)))) / sigma2
    
    
    # AC or RJ
    u <- runif(1)
    if (log(u) < logr) {
      rho <- rho_star
      Cnts <- Cnts + 1
    }
    
    # Collect 
    OUTPUT <- rbind(OUTPUT, cbind(t(beta), sigma2, rho))
    
    # print message
    cat(i, Cnts/i, beta, sigma2, rho,"\n") 
  }
  
  # AC ratio
  AC_ratio <- Cnts / n.sim
  
  # return
  if (Thinning) {
    thin.idx <- c(1, (1:(n.sim / thin)) * thin)
    return(list(Pst_results = OUTPUT[thin.idx, ], AC_ratio = AC_ratio))
  } else {
    return(list(Pst_results = OUTPUT, AC_ratio = AC_ratio))
  }
}


Res <- Metro_Gibbs_er_corr(beta0 = beta_ols, sigma2_0 = sigma2_ols, rho_0 = rho_ols, 
                    Data = y, DM = X, n.sim = 1e3, Thinning = F)


Res$AC_ratio # [1] 0.032
Res$AC_ratio 

thin = 10; n.sim <- 1e4; burn.in <- 1000
c(1, (1:(n.sim / thin)) * thin)

S <- 1000
odens<-S/1000

10000%%odens==0

c(burn.in, (1:(n.sim / thin)) * thin + burn.in)



