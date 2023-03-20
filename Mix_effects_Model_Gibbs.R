#==========================================
# Hoff Chp11 Linear & GL Mix effects Model
#==========================================

load("nelsSES.RData")
str(nels)
# 'data.frame':	1993 obs. of  4 variables:
#$ sch_id       : num  1011 1011 1011 1011 1011 ...
#$ sch_freelunch: num  6 6 6 6 6 6 6 6 6 6 ...
#$ stu_ses      : num  -0.0599 1.0517 -0.8635 -0.7966 -1.6135 ...
#$ stu_mathscore: num  52.1 57.6 66.4 44.7 40.6 ...


#-----------------------------------------------------------
# construct outcome and regressor variables for each school
#-----------------------------------------------------------
idx <- sort(unique(nels$sch_id))
m <- length(idx) # 100 schools

# for each school
Y <- list(); X <- list()
N <- c()
for (j in 1:m) {
  Y[[j]] <- nels[nels$sch_id == idx[j], 4]
  N[j] <- sum(nels$sch_id == idx[j])
  
  xj <- nels[nels$sch_id == idx[j], 3]
  xj <- xj - mean(xj)  # centered covariates for easy interepretation of Y
  X[[j]] <- cbind(rep(1, N[j]), xj)
}


N 
head(Y, 2)
head(X, 2)


#-------------------------
# Fit OLS for each school
#-------------------------

BETA.LS <- NULL
SIG2.LS <- NULL
for (j in 1:m) {
  lmfit <- lm(Y[[j]] ~ -1 + X[[j]])
  
  BETA.LS <- rbind(BETA.LS, lmfit$coef)
  SIG2.LS <- rbind(SIG2.LS, summary(lmfit)$sigma^2)
}

head(BETA.LS, 3)
head(SIG2.LS)


## for parameters in prior on theta, mu_0, 
BETA.mean <- apply(BETA.LS, 2, mean)
# 48.129449  2.373533 
BETA.var <- apply(BETA.LS, 2, var)
# 30.99446 10.09956 


## for parameters in prior on SIGMA, Lamb_O, S_0
BETA.cov <- cov(BETA.LS)
#             X[[j]]  X[[j]]xj
#X[[j]]   30.994456  2.195044
#X[[j]]xj  2.195044 10.099561


# for starting value of sigma2
SIGMA.var <- apply(SIG2.LS, 2, mean)
# [1] 75.90852


#-----------
# Set-up
#-----------

p <- dim(X[[1]])[2]

mu_0 <- apply(BETA.LS, 2, mean)
# 48.129449  2.373533 

L_0 <- cov(BETA.LS)
# 30.994456  2.195044
# 2.195044 10.099561

S_0 <- cov(BETA.LS)

eta_0 <- p + 2

a <- b <- 1


## Remark:
# the parameters in the prior distribution (hyper-prior)
  # uses ols estimates reflect the info is unbiased but 
  # weak. 

# For example, 
  # 1. the quantile of theta[, 2] slope under this 
  # hyper-prior setting, the prior of theta[, 2] slope has
  # 95% quantile of (-3.041309, 9.116819)
  # which is very large

theta_prior <- rmvnorm(100, mean = mu_0, sigma = L_0)
quantile(theta_prior[, 2], c(0.025, 0.975))
#      2.5%     97.5% 
# -3.041309  9.116819 


  # 2. prior distribution on Sigma ~ InvWishart(eta_0, S_0)
    # set S_0 = sample covariance cov(BETA.LS)
    # set eta_0 = p+2
    # to ensures E[Sigma] = 1/(eta_0 - p - 1) * S_0 = S_0 
      # which is the sample covariance matrix of the OLS estimates

  # 3. prior dist on sig2 ~ InvGamma(1, 1)
    # is essentially within school sample variance 
    # with very vague mean = a/b = 1



#-----------------------------------------
# Starting value of parameters of interest
#-----------------------------------------

beta_0 <- BETA.LS
#[1:100, 1:2]

theta_0 <- apply(BETA.LS, 2, mean) # 48.129449  2.373533 

Sigma_0 <- cov(BETA.LS) 
# 30.994456  2.195044
# 2.195044 10.099561

sig2_0 <- apply(SIG2.LS, 2, mean) # 75.90852


#--------------
# Gibbs sampler
#--------------
library(mvtnorm) # for rmvnrom
library(MCMCpack) # for inverse wishart
#riwish(v, S)

Gibbs_Mixeff <- function(beta_0, theta_0, Sigma_0, sig2_0, n.sim, 
                         Data = Y, DM = X, Thinning = T, thin = 10) {
  
  #idx <- sort(unique(nels$sch_id))
  #m <- length(idx) # 100 schools
  
  
  Sigma <- Sigma_0
  sig2 <- sig2_0
  theta <- theta_0
  
  BETA <- NULL
  THETA <- NULL
  SIGMA <- NULL
  SIG2 <- NULL
  for (i in 1:n.sim) {
    
    # update beta ~ MVN(mu_n, SIGMA_n) for each school
    Sigma_inv <- solve(Sigma)
    beta <- matrix(NA, nrow = m, ncol = 2)
    for (j in 1:m) {
      
      Sigma_n <- solve(Sigma_inv + t(X[[j]]) %*% X[[j]] / sig2)
      mu_n <- Sigma_n %*% (Sigma_inv %*% theta + t(X[[j]]) %*% Y[[j]] / sig2)
      
      #Sigma_n <- solve(Sigma_inv + t(X[[1]]) %*% X[[1]] / sig2)
      #mu_n <- Sigma_n %*% (Sigma_inv %*% theta + t(X[[1]]) %*% Y[[1]] / sig2)
      
      beta[j, ] <- rmvnorm(1, mean = mu_n, sigma = Sigma_n) # 1:2
      #head(beta)
    }
    
    
    # update theta ~ MVN(mu_n, L_n)
    L_0_inv <- solve(L_0)
    beta_bar <- apply(beta, 2, mean)
    
    L_n <- solve(L_0_inv + m * Sigma_inv)
    mu_theta_n <- L_n %*% (L_0_inv %*% mu_0 + m * Sigma_inv %*% beta_bar)
    
    theta <- rmvnorm(1, mean = mu_theta_n, sigma = L_n)
    theta <- t(theta)
    
    
    # update Sigma ~ InvWishart(eta_n, Sn_inv)
    eta_n <- eta_0 + m 
    
    S_theta <- (t(beta) - c(theta)) %*% t((t(beta) - c(theta)))
    Sn <- S_0 + S_theta
    
    Sigma <- riwish(v = eta_n, S = Sn) # 2*2
    
    
    # update sig2 ~ InvGamma(An, Bn)
    An <- a + sum(N) / 2
    
    SSR_j <- c()
    for (j in 1:m) {
      SSR_j[j] <- t(Y[[j]] - X[[j]] %*% beta[j, ]) %*% (Y[[j]] - X[[j]] %*% beta[j, ])
    }  
    SSR <- sum(SSR_j)
    Bn <- b + SSR / 2 
    
    sig2 <- 1/rgamma(1, shape = An, rate = Bn)
    
    
    # collect desired parameters of this iteration i without thinnng
    #BETA <- rbind(BETA, beta)       # nj *2
    #THETA <- rbind(THETA, t(theta)) # 1*2
    #SIGMA <- rbind(SIGMA, c(Sigma)) # 1*4
    #SIG2 <- rbind(SIG2, sig2)       # 1
    
    
    # return with thinning
    if (Thinning) {
      if (i %% thin == 0) {
        BETA = rbind(BETA, beta) 
        THETA = rbind(THETA, t(theta))
        SIGMA = rbind(SIGMA, c(Sigma))
        SIG2 = rbind(SIG2, sig2)
      }
    }
  }
  # return results after thinning
  return(list(Beta_pst = BETA, 
              Theta_pst = THETA, 
              Sigma_pst = SIGMA, 
              sig2_pst = SIG2))
}  


#-----------
# Test
#-----------

Res_500 <- Gibbs_Mixeff(beta_0 = beta_0, theta_0 = theta_0, Sigma_0 = Sigma_0, 
             sig2_0 = sig2_0, n.sim = 500, Data = Y, DM = X)

head(Res_500$Beta_pst)
head(Res_500$theta_pst)
head(Res_500$Sigma_pst)
head(Res_500$sig2_pst)


#---------
# Plots
#---------
par(mfrow = c(1, 2))

for (j in 1:2) {
  plot(1:50000, Res_500$Beta_pst[, j], type = "l", 
       xlab = "iterations", 
       ylab = bquote(beta[.(j-1)] ~ "| Data"),
       main = "")
}


#------------------------------
# check samples autocorrelation
#------------------------------
library(coda)
# lag-1
apply(Res_500$Beta_pst, 2, function(x) acf(x)$acf[2])
# [1] 0.06889570 0.02063625

# lag-10
apply(Res_500$Beta_pst, 2, function(x) acf(x)$acf[11])
# [1] 0.0792868 0.0805535


str(Res_500$theta_pst)
# num [1:500, 1:2] 
apply(Res_500$theta_pst, 2, function(x) acf(x)$acf[2])
# [1] 0.1362201 0.6024050
apply(Res_500$theta_pst, 2, function(x) acf(x)$acf[11])
# [1]  0.003052509 -0.006553851


apply(Res_500$Sigma_pst, 2, function(x) acf(x)$acf[2])
# [1] 0.3338999 0.7087971 0.7087971 0.7688635
apply(Res_500$Sigma_pst, 2, function(x) acf(x)$acf[11])
# [1] -0.03187579  0.13013764  0.13013764 0.16968983


apply(Res_500$sig2_pst, 2, function(x) acf(x)$acf[2])
# [1] 0.1649147
apply(Res_500$sig2_pst, 2, function(x) acf(x)$acf[11])
# [1] 0.04057485


## lag-1 auto-correlation of theta, esp. the 2nd theta, 
# and 4 Sigma_pst all have relatively high auto-correlations

## while lag-10 AC correlation of every sequence is low
## so thinning is indicated, and thin = 10. 


#----------------------
# Results with thinning
#----------------------

Res_thin10 <- Gibbs_Mixeff(beta_0 = beta_0, theta_0 = theta_0, Sigma_0 = Sigma_0,
             sig2_0 = sig2_0, n.sim = 5000, Data = Y, DM = X, 
             Thinning = T, thin = 10)

# str(Res_thin10$Beta_pst) # num [1:50000, 1:2]
# str(Res_thin10$sig2_pst) # num [1:500, 1]



#----------
# Plots
#----------

par(mfrow = c(1, 2))

for (j in 1:2) {
  plot(1:50000, Res_thin10$Beta_pst[, j], type = "l", 
       xlab = "iterations (scans) / 10", 
       ylab = bquote(beta[.(j-1)] ~ "|Data"),
       main = "")
}


for (j in 1:2) {
  plot(1:500, Res_thin10$Theta_pst[, j], type = "l", 
       xlab = "iterations (scans) / 10", 
       ylab = bquote(theta[.(j)] ~ "| Data"),
       main = "")
}


par(mfrow = c(2, 2))
for (j in 1:4) {
  plot(1:500, Res_thin10$Sigma_pst[, j], type = "l", 
       xlab = "iterations (scans) / 10", 
       ylab = bquote(Sigma[.(j)] ~ "| Data"),
       main = "")
}

par(mfrow = c(1, 1))
plot(1:500, Res_thin10$sig2_pst, type = "l", 
     xlab = "iterations (scans) / 10",
     ylab = bquote(sigma^2 ~ "| Data"),
     main = "")


## Density
par(mfrow = c(1, 2))
pdf("posterior density.pdf")
plot(density(Res_thin10$Theta_pst[, 2], adj = 2), 
     xlab = "slope parameter", 
     ylab = "Posterior density", 
     xlim = range(Res_thin10$Beta_pst[, 2]), 
     lwd = 2, main = "")

lines(density(Res_thin10$Beta_pst[, 2], adj = 2), 
      col = "gray", lwd = 2)

legend(-5, 1, legend = c(expression(theta[2]), expression(tilde(beta)[2])),
       lwd = c(2, 2), col = c("black", "gray"), 
       bty = "n")
dev.off()


# theta_pst[, 2]
  # posterior density of the slope for any one school

# beta[, 2]
  # posterior predictive distribution of slope



#--------------------
# Quantile of theta 2
#--------------------
quantile(Res_thin10$Theta_pst[, 2], c(.025, 0.975))
#     2.5%    97.5% 
# 1.817438 2.895650 

# compare to prior quantile:
quantile(theta_prior[, 2], c(0.025, 0.975))
#      2.5%     97.5% 
# -3.041309  9.116819 

# indicates a strong alteration of data information
# theta[2] 95% does not include 0 or negative values
  # indicates the majority school-level average slope
  # is positive


#--------------------------------
# Posterior distribution of beta
#--------------------------------

# samples from this distribution can be sample from 
  # MVN(theta_pst, Sigma_pst)
  # for each scan

# this posterior predictive distribution is more spread out
  # then theta's
  # reflecting heterogenity in slope across schools

# Pr(beta[, 2] < 0 | y, X) = 

sum(Res_thin10$Beta_pst[, 2] < 0) / 50000
# [1] 0.07858
# small but not negligible



head(Res_thin10$Beta_pst)
str(Res_thin10$Beta_pst)
# num [1:50000, 1:2]


#Beta_pst_mean <- Res_thin10$Beta_pst / 500 
# 500 samples for each school
# get the posterior mean of Beta
#str(Beta_pst_mean)
# num [1:50000, 1:2]
range(Beta_pst_mean[, 1])
range(Beta_pst_mean[, 2])

range(Res_thin10$Beta_pst[, 1])
range(Res_thin10$Beta_pst[, 2])

## Averaging 500 samplers for each of 100 schools
pst_avg <- NULL
for (i in 1:100) {
  total <- matrix(0, 1, 2)
  for (j in 0:4) {
    total = total + Res_thin10$Beta_pst[j * 100 + i, ]
  }
  pst_avg <- rbind(pst_mean, total/500)
}
str(pst_avg)
#  num [1:100, 1:2]
head(pst_avg)
range(pst_avg[, 1])
range(pst_avg[, 2])



plot(range(nels[, 3]), range(nels[, 4]), 
     xlab = "SES", ylab = "Math Score", 
     type = "n")
for (j in 1:m) {
  abline(a = Res_thin10$Beta_pst[j, 1], b = Res_thin10$Beta_pst[j, 2], 
         col = "gray")
}

abline(a = Res_thin10$Theta_pst[, 1], 
       b = Res_thin10$Theta_pst[, 2], 
       lwd = 2)

# compare to the original the school-wise regression lines
# the posterior regression line scattered more closely
# to each other 
# reflecting the heirarchical model is able to share
# information across groups and shrinks extreme regression 
# line to the across-group average

# in particular, hardly any regression line is negative
# when sharing info across groups. 









head(Y, 2)
head(X, 2)

#riwish(3, matrix(c(1,.3,.3,1),2,2))

matrix(c(1,.3,.3,1, 1,.3, .3,1),2,4) - matrix(c(0.1, 0.2), 2, 1)

t(beta[1, ])
c(0.1, 0.2)
X[[1]]
