#===================
# NormalGamma Model
#===================

# Data X1, X2, ..., Xn ~ N(mu, lambda^{-1})
  # mu, lambda unknown

# NormalGamma Model can be 
  # conjugate case: 
      # prior of mu and lambda is NormalGamma, 
        # (mu, lambda) ~  NormalGamma(m, c, a, b)

      # posterior mu and lambda is also Normal Gamma  
        # (mu, lambda | x1, ..., xn) ~ NormalGamma(M, C, A, B)
          # M = m * (c/c+n) + xbar * (n/c+n)  
          # C = c + n      //sample size used for estimating posterior mean mu|data
          # A = a + n/2   //posterior rate for lambda|data
          # B = b + 1/2 * sum(data_i - xbar)^2 + 1/2 (nc)/(n+c) (xbar - m)^2
            # = b + 1/2 * n * var(data) + 1/2 * (nc)/(n+c) (xbar - m)^2

  # semi-conjugate case:
    # pr


#----------------
# Conjugate case
#----------------
# given priors and data, 
  # return posteriors of parameters
  # by closed known form

priors = data.frame(m = 0, c = 1, a = 1/2, b = 50)

findPostPriors <- function(prior, data) {
  m = prior$m
  c = prior$c
  a = prior$a
  b = prior$b

  n = length(data)
  xbar = mean(data)
  s2 = var(data)
  
  #PostPar = NULL
  PostPar = data.frame(
    M = m * (c/(c + n)) + xbar * (n/(n + c)),
    C = c + n,
    A = a + n/2,
    B = b + 1/2 * (n - 1) * s2 + 1/2 * (n * c) / (n + c) * (xbar - m)
  )
  return(PostPar)
}

x <- c(18, 40, 15, 17, 20, 44, 38)
PostPars_post <- findPostPriors(prior = priors, data = x)
n.sims <- 1e4

# sample from posterior distribution of (mu, lambda | data)
  # 1. sample lambda ~ gamma(A, B)
  # 2. given lambda, sample mu|lambda ~ N(M, (C*lambda)^{-1})

lambda_post = rgamma(n.sims, shape = PostPars_post$A, rate = PostPars_post$B)
mu_post = sapply(sqrt((PostPars_post$C*lambda_post)^{-1}), rnorm, n = 1, mean = PostPars_post$M)

plot(mu_post, lambda_post)


mean(mu_post) # [1] 23.97186
mean(lambda_post) # [1] 0.007478751


#--------------------
# Semi-conjugate case
#--------------------
x <- c(18, 40, 15, 17, 20, 44, 38)
priors <- data.frame(mu_0 = 0, lam_0 = 1, a = 1/2, b = 50)
Semi_conj_Gibbs_post <- function(start.mu, start.lam, n.sim, priors, data) {
  
  print("This generates posterior of mu and lambda using Gibbs
        under semi-conjugate situation")
  
  Pars_post <- matrix(NA, nrow = n.sim, ncol = 2)
  Pars_post[1, ] <- c(start.mu, start.lam)
  
  mu_0 <- priors$mu_0
  lam_0 <- priors$lam_0
  a <- priors$a
  b <- priors$b
  
  xbar <- mean(data)
  S2 <- var(data)
  n <- length(data)
  
  for(i in 2:n.sim) {
    # sample from full conditional of mu
    mu_n <- mu_0 * lam_0 / (lam_0 + n * Pars_post[i-1, 2]) + 
      xbar * n * Pars_post[i-1, 2] / (n * Pars_post[i-1, 2] + lam_0)
    
    sd_n <- sqrt((lam_0 + n * Pars_post[i-1, 2])^{-1})
    Pars_post[i, 1] <- rnorm(1, mean = mu_n, sd = sd_n)
    
    # sample from full conditional of lam
    A = a + n/2
    B = 1/2 * ((n-1) * S2 + n * (xbar - Pars_post[i, 1])^2 + 2*b)
    Pars_post[i, 2] <- rgamma(1, shape = A, rate = B)
  }
  return(Pars_post)
}

semi_post <- Semi_conj_Gibbs_post(start.mu = 1, start.lam = 1, 
                                  priors = priors, 
                                  n.sim = 1e4, data = x)

head(semi_post)
mu_post <- semi_post[, 1]
lam_post <- semi_post[, 2]

mean(mu_post) # [1] 0.2356989
mean(lam_post) # [1] 0.001388827


priors2 <- data.frame(mu_0 = 0, lam_0 = 1, a = 1/2, b = 10)
semi_post2 <- Semi_conj_Gibbs_post(start.mu = 1, start.lam = 1, priors = priors2,
                     n.sim = 1e4, data = x)

mean(semi_post2[, 1]) # [1] 0.2598341
mean(semi_post2[, 2]) # [1] 0.001402335


priors3 <- data.frame(mu_0 = 0, lam_0 = 1, a = 1/2, b = 1)
semi_post3 <- Semi_conj_Gibbs_post(start.mu = 1, start.lam = 1, priors = priors3,
                                   n.sim = 1e4, data = x)

mean(semi_post3[, 1]) # [1] 0.2502062
mean(semi_post3[, 2]) # [1] 0.0014111



## data in Hoff P94
y <- c( 1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
mean(y) # [1] 1.804444
var(y) # [1] 0.01687778 sigma^2
sqrt(var(y)) # [1] 0.1299145
1/var(y) # [1] 59.24951  lambda = 1/sigma^2
# for prior parameter of lambda
  # E[lamda] = a/b ~= 60
  # so if a set to 1, b ~= 1/60
# 1/60 = [1] 0.01666667

priors <- data.frame(mu_0 = 1, lam_0 = 1, a = 1, b = 1)
mu_lam_pst <- Semi_conj_Gibbs_post(start.mu = 1, start.lam = 1, n.sim = 1e4, 
                     priors = priors, data = y)

head(mu_lam_pst)
mean(mu_lam_pst[, 1]) # [1] 1.760803
mean(mu_lam_pst[, 2]) # [1] 2.422582 precision on average
1/mean(mu_lam_pst[, 2]) # [1] 0.4123006 var on average


priors2 <- data.frame(mu_0 = 1, lam_0 = 1, a = 1, b = 0.01)
mu_lam_pst2 <- Semi_conj_Gibbs_post(start.mu = mean(y), start.lam = 1/var(y), n.sim = 1e4, 
                                   priors = priors2, data = y)

mean(mu_lam_pst2[, 1]) # [1] 1.801969
mean(mu_lam_pst2[, 2]) # [1] 64.47763


#--------------------------
# Compare b = 1 vs b = 0.01
#--------------------------

par(mfrow = c(1, 2))
plot(mu_lam_pst[, 1], mu_lam_pst[, 2], main = "b = 1")
plot(mu_lam_pst2[, 1], mu_lam_pst2[, 2], main = "b = 0.01")

# conclusion: how to set rate parameter in gamma distri for precison lambda
  # y <- c( 1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
  # var(y) = sigma^2 = 1/lambda =  [1] 0.01687778
  # lambda on average shall be 1/0.01687778 ~= 60
  # by formula, E[lambda] = a/b
    # so if set a = 1, b ~= 1/60 ~= 0.016


#----------
# plots
#----------
plot(mu_lam_pst[, 1], mu_lam_pst[, 2])

# trace plots
plot(mu_lam_pst[, 1], type = "l", main = "Trace plot of mu_post")
plot(mu_lam_pst[, 2], type = "l", main = "Trace plot of lambda_post")


hist(mu_lam_pst[, 1], probability = T)
lines(density(mu_lam_pst[, 1]), col = "red")

hist(mu_lam_pst[, 2], probability = T, 
     main = "Histogram of precision lambda")
lines(density(mu_lam_pst[, 2]), col = "red")


#--------------------------------------
# empirical quantiles of Gibbs samples
#--------------------------------------

## b = 1
# quantile for mu_post
quantile(mu_lam_pst[, 1], c(.025, .5, .975))
#     2.5%      50%    97.5% 
# 1.286193 1.764634 2.218488 

# quantile for lambda_post
quantile(mu_lam_pst[, 2], c(.025, .5, .975))
#      2.5%       50%     97.5% 
#   1.516578 4.338744 9.521027 

# CI for population standard deviation
quantile(sqrt(mu_lam_pst[, 2]^{-1}), c(.025, .5, .975))
#   2.5%       50%     97.5% 
# 0.3240844 0.4800849 0.8120215 


## b = 0.01
quantile(mu_lam_pst2[, 1], c(.025, .5, .975))
#    2.5%      50%    97.5% 
# 1.707953 1.803137 1.895942

quantile(mu_lam_pst2[, 2], c(.025, .5, .975))
#      2.5%       50%     97.5% 
# 20.45162  60.17795 133.14719 

quantile(sqrt(mu_lam_pst2[, 2]^{-1}), c(.025, .5, .975))
#      2.5%        50%      97.5% 
# 0.08666306 0.12890842 0.22112416 
# sqrt(var(y)) # [1] 0.1299145


