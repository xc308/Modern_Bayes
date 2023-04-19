#==========================
# Ordinal probit regression
#==========================

# social mobility data

load("socmob.RData")
head(socmob)


yDEG <- socmob$DEGREE
unique(yDEG)
# [1]  1  0  3  4  2 NA

yCHILD <- socmob$CHILDREN
unique(yCHILD)
# [1]  3  1  2  0  4  5 NA  8  6  7

yPDEG <- socmob$PDEGREE 

unique(yPDEG)
yPDEG <- 1*(yPDEG > 2) # change logical into numerical


#------------------------------------------
# Setting up for Ordinal probit regression
#------------------------------------------

# design matrix X
X <- cbind(yCHILD, yPDEG, yCHILD * yPDEG)
y <- yDEG + 1

# remove rows with NA either in X or y
keep <- (1:length(y))[!is.na(apply(cbind(y, X), 1, sum))] # 959
X <- X[keep, ]
y <- y[keep]


# dim
n <- dim(X)[1]
p <- dim(X)[2]


# latent process z1,...,zn of y1,...,yn
  # by yi = g(zi) = N(zi; 0, 1)
z <- qnorm(rank(y, ties.method = "random") / (n + 1))


# K: catogory number of values of y
K <- length(sort(unique(y)))
# 5

K - 1
# [1]  4


#--------------------------
# parameters in the priors
#--------------------------

# in g1, ... g_{k-1}|. ~iid  N_{ak, bk}(mu_k, sigma_k^2)
mu_k = 0; sigma_k = 100

# g1, .., g_{k-1}
mu = rep(0, (K-1))
sigma = rep(100, (K-1))


#----------------------------------------
# starting values of unknowns of interest
#----------------------------------------

z_start <- qnorm(rank(y, ties.method = "random") / (n + 1))
# 1:959

beta_start <- rep(0, p) # p = 3


#--------------
# Gibbs Sampler
#--------------
iXX <- solve(t(X) %*% X) # [1:3,1:3]
Sig <- iXX * (n / (n+1))
choSig <- chol(Sig)


g <- rep(NA, K-1)
BETA <- matrix(NA, 1000, p)
Z <- matrix(NA, 1000, n)

ac <- 0
S <- 25000

z <- z_start

for (s in 1:S) {
  
  # update g1,...,g_{k-1} |. ~ N_(ak, bk)(muk, sigmak)
  for (k in 1:(K-1)) {
 
    ak = max(z[y == k])
    bk = min(z[y == (k + 1)])
    
    u = runif(1, pnorm((ak - mu[k]) / sigma[k]), 
              pnorm((bk - mu[k]) / sigma[k]))
    
    g[k] = mu[k] + sigma[k] * qnorm(u)
  }
  
  
  # update beta |. ~ N(mu_n, Sig_n)
  mu_n = Sig %*% (t(X) %*% z) 
  beta = ChoSig %*% rnorm(p) + mu_n  # 3*1
  
  
  # update z1,...,zn|. ~ N_{a, b}(t(beta)%*%X, 1)
  mu_z = X %*% beta  # n * 1
  
  a = c(-Inf, g)[match(y-1, 0:(K-1))] # n*1
  #head(y-1) # [1] 1 0 1 3 3 4
  #0:(K-1) # [1] 0 1 2 3 4
  #match(y-1, 0:(K-1)) : 2, 1, 2, 4, 5 the corresponding K location of y-1
  # indexing c(-Inf, g) with index 2, 1, 2, 4, 5
  
  b = c(g, Inf)[match(y, 1:K)] # n*1
  #y == match(y, 1:K) # all true
  
  # pnorm(Inf) # [1] 1 so even when b = Inf, a = -Inf,
  # pnorm(-Inf) # [1] 0
  
  u = runif(n, pnorm(a - mu_z), pnorm(b - mu_z)) # n*1
  z = mu_z + 1 * qnorm(u) # n*1
  
  
  
  # help mixing -- Metropolis sample of z and g
  # propose
  c = rnorm(1, 0, 0.1)
  z_star = z + c; g_star = g + c
  
  # logr
  logr = sum(dnorm(z_star, mu_z, 1, log = T) - dnorm(z, mu_z, 1, log = T)) + 
    sum(dnorm(g_star, mu, sigma, log = T) - dnorm(g, mu, sigma, log = T))
  
  # AC
  if (log(runif(1)) < logr) {
    z = z_star
    g = g_star
    
    ac = ac + 1
  }
  
  
  # collect, every 25th (S/1000) sample 
  if (s %% (S/1000) == 0) {
    
    cat(s, ac/s, "\n")
    
    BETA[s/(S/1000), ] = beta
    Z[s/(S/1000), ] = z
  }
}


#--------
# Results
#--------

# 25000 0.38496 

str(BETA) # num [1:1000, 1:3]
beta_pst <- apply(BETA, 2, mean)
# [1] -0.02166479  0.82628568  0.07401959

BETA_CI <- apply(BETA, 2, quantile, prob = c(0.25, 0.5, 0.975))
t(BETA_CI)


#----------------
# Regression Line
#----------------

# X1 = CHILD; X2 = PDEG; X3 = CHILD*PDEG
# E[z | y1,...,yn, X1, X2, X3] 
  # = -0.02500787 CHILD, if PDEG = 0;
  # = -0.02500787 CHILD + 0.81613992 * (PDEG = 1) + 0.07932693 CHILD * (PDEG = 1)
    # =  0.816 + 0.054 CHILD, if PDEG = 1



#-------
# Plots
#-------

str(Z) # num [1:1000, 1:959]
str(X) # num [1:959, 1:3]
X[, 2] # 0, 1
# pch = 15 + X[, 2]: pch = 15(PDEG = 0), 16(PDEG = 1)


pdf("Ordinal_Regression_and_Posterior_of_Beta3. pdf", 
    family = "Times")
par(mfrow = c(1, 2), mar = c(3, 3, 1, 1), mgp = c(1.75, 0.75, 0))

plot(X[,1] + 0.25*X[, 2], Z[1000, ],
     pch = 15 + X[, 2], col=c("gray","black")[X[,2]+1],
     xlim = c(0, 8), ylim = round(range(Z[1000, ])),
     xlab = "number of children", 
     ylab = "z")


z_PM <- apply(Z, 2, mean)
beta_PM <- apply(BETA, 2, mean)

abline(a = 0, b = beta_PM[1], col = "gray", lwd = 2 )
abline(a = beta_PM[2], b = beta_PM[1] + beta_PM[3], col = "black", lwd = 2)

legend(c(5, 4), legend = c("PDEG = 0", "PDEG = 1"), 
       pch = c(15, 16), col = c("gray", "black"))



plot(density(BETA[, 3], adj = 2), 
     xlab = expression(beta)[3], 
     ylab = "density",
     xlim = c(-0.5, 0.5), 
     main = "")


## prior on beta ~ N(mu0, Sigma_0)
  # mu_0 = 0
  # Sigma_0 via unit infomation
    # inv_Sigma_0 = t(X) %*% X / n
    # Sigma_0 = n * t(X) %*% X

beta3_sd_prior <- sqrt(n * (t(X) %*% X) [3, 3]) # [1] 732.8301

x <- seq(-0.7, 0.7, length = 100)
lines(x, dnorm(x, mean = 0, sd = beta3_sd_prior), 
      col = "gray", lwd = 2)

#range(dnorm(x, mean = 0, sd = beta3_sd_prior)) 
# [1] 0.0005443855 0.0005443858
# a flat prior

legend(c(-0.3, 5.5), legend = c("prior", "posterior"),
      lwd = c(2, 2), col = c("gray", "black"), bty = "n")

dev.off()





