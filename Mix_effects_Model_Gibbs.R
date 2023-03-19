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

Gibbs_Mixeff <- function(beta_0, theta_0, Sigma_0, sig2_0, n.sim, Data = Y, DM = X) {
  
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
    Sn_inv <- solve(Sn)
    
    Sigma <- riwish(v = eta_n, S = Sn_inv) # 2*2
    
    
    # update sig2 ~ InvGamma(An, Bn)
    An <- a + sum(N) / 2
    
    SSR_j <- c()
    for (j in 1:m) {
      SSR_j[j] <- t(Y[[j]] - X[[j]] %*% beta[j, ]) %*% (Y[[j]] - X[[j]] %*% beta[j, ])
    }  
    SSR <- sum(SSR_j)
    Bn <- b + SSR / 2 
    
    sig2 <- 1/rgamma(1, shape = An, rate = Bn)
    
    
    # collect desired parameters of this iteration i
    BETA <- rbind(BETA, beta)       # nj *2
    THETA <- rbind(THETA, t(theta)) # 1*2
    SIGMA <- rbind(SIGMA, c(Sigma)) # 1*4
    SIG2 <- rbind(SIG2, sig2)       # 1
  }
  
  # return
  return(list(Beta_pst = BETA, theta_pst = THETA, 
              Sigma_pst = SIGMA, sig2_pst = SIG2))
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


head(Y, 2)
head(X, 2)

#riwish(3, matrix(c(1,.3,.3,1),2,2))

matrix(c(1,.3,.3,1, 1,.3, .3,1),2,4) - matrix(c(0.1, 0.2), 2, 1)

t(beta[1, ])
c(0.1, 0.2)
X[[1]]
