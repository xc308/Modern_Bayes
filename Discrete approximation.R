#============================
# 6.2 Discrete approximation
#============================

# evaluates posterior distribution
  # p(theta, lambda| y1,..., yn) on 100*100 grid of  
    # evenly spaced parameter values

# Model;
  # Y1, ..., Yn | theta, sigma2 ~ N(theta, sigma2)
  # theta|sigma2 ~ N(mu0, sigma2/k0)
  # 1/sigma2 = lambda ~ Gamma(a, b)


# if prior joint distribution is discrete on grid
  # then p(theta_k, lamb_l | y1,..., yn) / 
            # Sum_{g=1}^G Sum_{h = 1}^H p(theta_g, lamb_h | y1, ..., yn)
  # is the actual posterior dist of {theta_k, lamb_l}


y <- c( 1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
mean(y) # [1] 1.804444
var(y) # [1] 0.01687778
1/var(y) # [1] 59.24951
# E[lamb] = a/b = 60; a = 1
# b = 1/60 = 0.0168


# set prior:
# mu0 = 1.9; sigma2 = 1; 
# a = 1, b = 0.01

mu0 <- 1.9
sigma2 <- 1
a = 1
b = 0.01

G = 100; H = 100

mean.grid <- seq(1.5, 2, length = G)
pres.grid <- seq(17.5, 175, length = H)

#quantile(mu_lam_pst2[, 1], c(.025, .5, .975))
#    2.5%      50%    97.5% 
# 1.707953 1.803137 1.895942

#quantile(mu_lam_pst2[, 2], c(.025, .5, .975))
#      2.5%       50%     97.5% 
# 20.45162  60.17795 133.14719 

post.grid <- matrix(NA, nrow = G, ncol = H)
for (g in 1:G) {
  for (h in 1:H) {
    post.grid[g, h] <- 
      dnorm(mean.grid[g], mu0, sqrt(sigma2)) * 
      dgamma(pres.grid[h], a, b) *
      prod(dnorm(y, mean.grid[g], sqrt(1/pres.grid[h])))
  }
}
post.grid <- post.grid / sum(post.grid)


theta.marg <- rowSums(post.grid)
pres.marg <- colSums(post.grid)


hist(rgamma(100, a, b), probability = T)
lines(density(rgamma(100, a, b)), col = "red")


## Remark;
  # evaluate two-parameter posterior distribution on 100 grid of each parameter
    # require grid size 100*100 = 100^2
  # in general, to evaluate p-parameter posterior dist on 100 grid 
    # of each of the parameters require grid size 100^p
# so discrete evaluation will only be feasible for a small number
  # of parameters 










    



















