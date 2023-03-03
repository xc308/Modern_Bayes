#=================
# MVN Normal Model
#=================

# This is function samples unknown mean vector theta (2*1) 
# and unknown covariance matrix (2*2) in a MVN model

# likelihood: 
  # y1,..., yn | theta, Sigma ~ MVN(theta, Sigma)

# priors:
  # theta ~ MVN(mu_0, Lamb_0)
  # Sigma ~ InvWishart(nu_0, S_0^{-1})  # generalization of InvGamma

# joint posterior:
  # theta, Sigma | y1,...,yn ~ no closed form

# but full conditionals of either theta or Sigma are known

# full conditionals of theta:
  # theta | y1,...,yn, Sigma ~ MVN(mu_n, Lamb_n)
    # mu_n : 2 * 1 vector 
    # Lamb_n : 2*2 matrix 

  # Lamb_n = solve(solve(Lamb0) + n*solve(Sigma))  # (.+.)^{inv}
  # mu_n = Lamb_n %*% (solve(Lamb0) %*% mu_0 + n*solve(Sigma) %*% ybar)

      # ybar = apply(Data, 2, mean) # Data is n*2


# full conditionals of Sigma:
  # Sigma | y1,...,yn, theta ~ InvWishart(nu_n, S_n^{-1})
      # nu_n: nu_0 + n
      # S_n^{-1}: S_0 + S_theta

        # S_theta = (t(Data) - theta) %*% t(t(Data) - theta)


#---------------------------------------
# Parameters in the prior -- hyperprior
#---------------------------------------
# theta ~ MVN(mu_0, Lamb_0)

# mu_0: 2*1 vector
mu_0 <- c(50, 50) # create a Vector in R using c() 

# Lambda_0: 2*2 matrix
  # Lamb_011, Lamb_022: by score must in [0, 100], 
    # by normality 95% of normal data will be in mu +/- 2 std
    # so set sqrt(Lamb_011) = 25 = 50/2, 
      # then 95% of data will be in [0, 100]
      # Pr(mu_{0j} !in [0, 100]) = 0.05, j = 1, 2

  # Lamb_012, Lamb_021 : two tests are similar, so correlation = 0.5
    # Lamb_012/Lamb_011 = 0.5
    # Lamb_012 = Lamb_021 = 312.5

Lamb_0 <- matrix(c(625, 312.5, 312.5, 625), nrow = 2)


## Sigma ~ InvWishart (nu_0, S_0^{-1})
  # S_0 = Lamb_0
  # want the Sigma loosely centered around S_0, i.e., Lamb_0
  # by E[Sigma] = 1/(nu_0 - p - 1) * S_0, 
    # nu_0 = p + 2 = 4

S_0 <- Lamb_0
nu_0 <- 4


#-------
# Data 
#-------
Y <- structure(c(59, 43, 34, 32, 42, 38, 55, 67, 64,
                 45, 49, 72, 34, 70, 34, 50, 41, 52,
                 60, 34, 28, 35, 77, 39, 46, 26, 38,
                 43, 68, 86, 77, 60, 50, 59, 38, 48,
                 55, 58, 54, 60, 75, 47, 48, 33),
               .Dim = c(22L, 2L), .Dimnames = list(NULL,
                                                   c("pretest", "posttest")))



head(Y)
ybar <- apply(Y, 2, mean)
str(ybar)
# Named num [1:2] 

Lamb_0 %*% ybar
#          [,1]
#[1,] 46321.02
#[2,] 48409.09



#----------------
# Gibbs Sampler
#----------------
install.packages("mvtnorm")
library(mvtnorm)

install.packages("stats")
library(stats) # for rWishart


MVN_Gibbs_sampler <- function(theta_start, Sigma_start, n.sim, Data) {
  
  THETA <- SIGMA <- NULL
  
  Sigma <- Sigma_start
  ybar <- apply(Data, 2, mean)
  
  for (i in 1:n.sim) {
    
    n = nrow(Data)
    
    # sample theta | Data, Sigma
    Lamb_n = solve(solve(Lamb_0) + n * solve(Sigma)) # (.)^{-1}
    mu_n = Lamb_n %*% (solve(Lamb_0) %*% mu_0 + (n * solve(Sigma)) %*% ybar)
    
    theta <- rmvnorm(1, mean = mu_n, sigma = Lamb_n) # a matrix [1, 1:2]
    theta <- c(theta) # vectorize into 2*1 for below minus calculation
    
    
    # sample Sigma | Data, theta ~ InvWishart(nu_n, S_n^{-1})
    nu_n = nu_0 + n 
    S_theta = (t(Data) - theta) %*% t(t(Data) - theta)  # 2*2
    S_n = S_0 + S_theta
    S_n_inv = solve(S_n)
    
    Sigma <- solve(rWishart(1, df = nu_n, Sigma = S_n_inv)[, , 1])
    
    
    # collect post-samplers of theta and Sigma for this run i
    THETA <- rbind(THETA, theta)
    SIGMA <- rbind(SIGMA, c(Sigma))
  }
  return(list(theta = THETA, Sigma = SIGMA))
}


#-------------------------
# starting value setting
#-------------------------

mean(Y)
# [1] 50.52273
theta_start = c(50, 50)

cov(Y)
#           pretest posttest
# pretest  182.1558 148.4069
# posttest 148.4069 243.6472

Sigma_start <- matrix(c(150, 100, 100, 200), nrow = 2)
solve(Sigma_start)


#----------
# Test on data
#----------

Res <- MVN_Gibbs_sampler(theta_start = theta_start, Sigma_start = Sigma_start, 
                  n.sim = 5000, Data = Y)

head(Res$theta)
#           [,1]     [,2]
#theta 48.28698 58.54779
#theta 41.72282 49.96556
#theta 49.98020 54.64395

Res_theta <- Res$theta


head(Res$Sigma)
#         [,1]     [,2]     [,3]     [,4]
# [1,] 253.6271 204.6613 204.6613 255.0841
# [2,] 270.7554 232.0677 232.0677 412.7416
# [3,] 127.0593 116.7394 116.7394 222.4469
# [4,] 211.6924 165.1282 165.1282 246.6646
# [5,] 132.5914  73.0294  73.0294 126.4443
# [6,] 158.1921 127.3283 127.3283 229.6900

Res_Sigma <- Res$Sigma




#---------------
# Trace plots 
#---------------

## Trace plots of theta 
plot(1:5000, Res_theta[, 1], pch = 16, cex = 0.35, 
     xlab = "Iterations", ylab = expression(theta[1]), 
     main = expression(paste("Trace plot ", theta[1])))


par(mfrow = c(1, 2))
for(idx in 1:2) {
  plot(1:5000, Res_theta[, idx], pch = 16, cex = 0.35, 
       xlab = "Iterations", ylab = bquote(theta[.(idx)]), 
       main = bquote("Trace plot "~ theta[.(idx)]))
  
}


## density plots of theta
pdf("density plt_theta1&2.pdf")
for(idx in 1:2) {
  plot(density(Res_theta[, idx]), 
       main = bquote("Density plot of "~ theta[.(idx)]))
  abline(v = mean(Res_theta[, idx]), col = "dark red")
  abline(v = quantile(Res_theta[, idx], c(0.025, 0.975)), 
         col = "blue", lty = 2)
}
dev.off()


# running average of theta
n.sim <- 5000
total <- 0
Run_Avg <- c()
for (i in 1:n.sim) {
  total = Res_theta[, 1][i] + total
  Run_Avg[i] = total / i 
}

plot(1:n.sim, Run_Avg, type = "l", 
     main = expression(paste("Running Average of ", theta[1])))





source("Fun_Running_Avg.R")
#Running_Avg(n.sim = 5000, total = 0, pst_smpl = Res_theta[, 1])

## Theta

pdf("Running_Avg_theta1_2.pdf")

par(mfrow = c(1, 2))
n.sim = 5000
Run_Avg <- matrix(NA, nrow = n.sim, ncol = 2)
for (idx in 1:2) {
  # calculate Running avg
  Run_Avg[, idx] <- Running_Avg(n.sim = 5000, total = 0, 
                                pst_smpl = Res_theta[, idx])
  
  
  # plot
  plot(1:n.sim, Run_Avg[, idx], type = "l", 
       xlab = "Iterations",
       ylab = bquote("Running Avg of " ~ theta[.(idx)]),
       main = bquote("Running Average of " ~ theta[.(idx)]))
  
}
dev.off()




#---------------------------
# Summary Statistics 
#---------------------------

## theta ##

## posterior mean 
pst_mean <- apply(Res_theta, 2, mean)
# [1] 47.08409 53.73396


## posterior Credible Intervals for two theta 
CI_theta <- apply(Res_theta, 2, function(x) quantile(x, c(0.025, 0.975)))
#        [,1]     [,2]
#2.5%  41.18628 47.21365
#97.5% 52.96473 60.65773

CI_theta <- t(CI_theta)
#         2.5%    97.5%
# [1,] 41.18628 52.96473
# [2,] 47.21365 60.65773

pst_summary_theta <- cbind(pst_mean, CI_theta)
rownames(pst_summary_theta) <- c("theta1", "theta2")



## Sigma ##

head(Res$Sigma)
#         [,1]     [,2]     [,3]     [,4]
# [1,] 253.6271 204.6613 204.6613 255.0841
# [2,] 270.7554 232.0677 232.0677 412.7416
# [3,] 127.0593 116.7394 116.7394 222.4469
# [4,] 211.6924 165.1282 165.1282 246.6646
# [5,] 132.5914  73.0294  73.0294 126.4443
# [6,] 158.1921 127.3283 127.3283 229.6900

Res_Sigma <- Res$Sigma

## posterior mean 
pst_mean_Sigma <- apply(Res_Sigma, 2, mean)
# [1] 201.0628 154.9209 154.9209 259.1341
matrix(pst_mean_Sigma, nrow = 2)
#         [,1]     [,2]
# [1,] 201.0628 154.9209
# [2,] 154.9209 259.1341


CI_Sigma <- apply(Res_Sigma, 2, function(x) quantile(x, c(0.025, 0.975)))
#          [,1]      [,2]      [,3]     [,4]
#2.5%  114.0343  67.34827  67.34827 146.6702
#97.5% 357.1705 304.77362 304.77362 456.7223

CI_Sigma <- t(CI_Sigma)
#          2.5%    97.5%
#[1,] 114.03430 357.1705
#[2,]  67.34827 304.7736
#[3,]  67.34827 304.7736
#[4,] 146.67018 456.7223

Summary_stats_Sigma <- cbind(pst_mean_Sigma, CI_Sigma)

rownames(Summary_stats_Sigma) <- c("Sig_11", 
                                   "Sig_12", 
                                   "Sig_21", 
                                   "Sig_22")


#--------------------
# Posterior Inference
#--------------------

# approxiamte posterior probabilities and confidence regions
# 
# if give this test to a larger population, 
# want to know if the average score of the 2nd one 
  # will be higher than the 1st one

# Pr(theta2 > theta1 | Data) = ?

mean(Res_theta[, 2] > Res_theta[, 1])
# [1] 0.9944























