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
load("icecore.RData")

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
b = 1/2

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
    beta_n <- Sigma_n %*% (Sigma0_inv %*% beta0 + t(X) %*% solve(C_rho) %*% y / sigma2)
    
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
    
    # Collect with or w/o thinning
    if (Thinning) {
      if (i %% thin == 0) { # i is a multiple of thin, e.g.25, collect
        OUTPUT <- rbind(OUTPUT, cbind(t(beta), sigma2, rho))
        cat(i, Cnts/i,"\n")
      }
    } else {
      OUTPUT <- rbind(OUTPUT, cbind(t(beta), sigma2, rho))
      cat(i, Cnts/i,"\n")
    }
  }
  
  # AC ratio
  AC_ratio <- Cnts / n.sim
  
  # return
  return(list(Pst_results = OUTPUT, AC_ratio = AC_ratio))
  
  
  #if (Thinning) {
    #thin.idx <- c(1, (1:(n.sim / thin)) * thin)
    #return(list(Pst_results = OUTPUT[thin.idx, ], AC_ratio = AC_ratio))
  #} else {
   # return(list(Pst_results = OUTPUT, AC_ratio = AC_ratio))
  #}
  
}


#---------------------------
# n.sim = 1000 w/o thinning
#---------------------------
Res <- Metro_Gibbs_er_corr(beta0 = beta_ols, sigma2_0 = sigma2_ols, 
                           rho_0 = rho_ols, 
                    Data = y, DM = X, n.sim = 1e3, Thinning = F)


Res$AC_ratio  # [1] 0.258 Acceptable
Results_1000 <- Res$Pst_results

pst_rho <- Results_1000[, 4]

par(mfrow = c(1, 2))
plot(1:1e3, pst_rho, type = "l", 
     xlab = "scan", ylab = expression(rho),
     main = bquote("Trace plot of " ~ rho))

acf(pst_rho, ci.col = "gray")
acf(pst_rho)$acf[2]
# [1] 0.9616772

library(coda)
effectiveSize(pst_rho)
# 19.51615  out of 1000

## simulated samples are highly correlated 
  # with lag-1 auto correlation 0.96
  # eff size 19 out of 1000


apply(Results_1000, 2, effectiveSize)
#                       sigma2       rho 
# 40.640169 40.695030  8.752703 19.516152


#-------------------------------------
# n.sim = 1e4, Thinning = T, thin = 25
#-------------------------------------

Res_1e4 <- Metro_Gibbs_er_corr(beta0 = beta_ols, sigma2_0 = sigma2_ols, 
                               rho_0 = rho_ols, Data = y, DM = X, 
                               n.sim = 1e4, Thinning = T, thin = 25)


Res_1e4$AC_ratio #[1] 0.2817 
pst_rho_1e4 <- Res_1e4$Pst_results[, 4]


plot(1:400, pst_rho_1e4, xlab = "scan/25", 
     ylab = expression(rho), type = "l")

acf(pst_rho_1e4, xlab = "lag/25") # lag 25th

# to get lag-1 of samples after thinning
acf(pst_rho_1e4, plot = F, )$acf[2] 
# [1] 0.4488623


#------
# Plots
#------
library(coda)

pdf("Trace_ACF_MH.pdf")
layout(mat = matrix(c(1,2), nrow = 1, ncol = 2))

pdf("Trace_ACF_MH_2.pdf")
layout(mat = matrix(c(1,2, 3, 4), nrow = 2, ncol = 2))

pdf("Trace_ACF_MH_3.pdf")
layout(mat = matrix(c(1,2, 3, 4), nrow = 2, ncol = 2, byrow = T))

for (idx in 1:2) {
  plot(100:400, Res_1e4$Pst_results[100:400, idx], 
       xlab = "scan/25", 
       ylab = bquote(beta[.(idx)] ~ "| Data"),
       type = "l", 
       main = bquote("Trace plot of " ~ beta[.(idx)]))

  ACF <- acf(Res_1e4$Pst_results[, idx], plot = F)$acf[2]
  EFF <- effectiveSize(Res_1e4$Pst_results[, idx])
  
  acf(Res_1e4$Pst_results[, idx], xlab = "lag/25",
      ci.col = "gray", 
      main = bquote(atop("ACF of " ~ beta[.(idx)],
                         atop("lag-1 acf: " ~ .(round(ACF, 2)),
                         "Effsize: " ~ .(round(EFF, 2))))),
      cex = 0.35)
}
dev.off()



#--------------------
# summary statistics
#--------------------

# posterior mean
pst_mean <- apply(Res_1e4$Pst_results, 2, mean)
pst_mean[2]

# posterior quantiles
CI <- apply(Res_1e4$Pst_results, 2, function(x) quantile(x, c(0.025, 0.975)))
CI_t <- t(CI)
CI_t[2,]

Sum_Stats <- cbind(pst_mean, CI_t)
rownames(Sum_Stats) <- c("beta_1", "beta_2", "sigma2", "rho")
Sum_Stats <- round(Sum_Stats, 3)

library(xtable)
xtable(Sum_Stats)


#---------
# Remarks
#---------

# posterior of beta-2 (slop) is 0.029
# still indicating positive relationship between 
  # temp and CO2, 
  # but not as strong as the one from OLS (0.07985291)

#beta_ols
# (Intercept)  icecore$co2 
#-23.02414043   0.07985291 

pdf("Compare_OLS_GLM.pdf", family = "Times")
layout(mat = matrix(c(1, 2), nrow = 1, ncol = 2))
plot(density(Res_1e4$Pst_results[, 2], adj = 2), # adjust  for smooth 
     xlab = expression(beta[2]), 
     ylab = "posterior marginal density", 
     main = "")
abline(v = pst_mean[2], lwd = 1.5, col = "red")
abline(v = CI_t[2,], lwd = 1, col = "gray")


plot(icecore$co2, icecore$tmp, xlab = "CO2", ylab = "Temp")
abline(a = beta_ols[1], b = beta_ols[2], col = "gray", lwd = 2)
abline(a = pst_mean[1], b = pst_mean[2], col = "black", lwd = 2)
legend(180, 2.5, legend = c("OLS", "BayesGLM"),
       col = c("gray", "black"), bty = "n", lwd = c(2, 2),
       cex = 0.85)
dev.off()

## GLS put a weight to data according to temporal correlation
  # the further temporal apart, the less correlation
  # so some data point although with high Temp (y) values
  # but is temporally far away, so less weighted
  # weaker slope

## While OLS treated different data with the same weight
  # easier to be influenced by those high-value y. 





