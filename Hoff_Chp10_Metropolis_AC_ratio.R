#========================
# Metropolis algo in GLM
#========================

# Using know Normal model and normal prior with known variance

# y1,...,yn | theta ~ iid N(theta, sigma2); sigma2 known
  # theta ~ N(mu_0, tau_0), 
    # theta | y1,...,yn ~ N(mu_n, tau_n)
      
      # tau_n = (1/tau_0 + n/sigma2)^{-1}  # variance
      # mu_n = (mu_0 * 1/tau_0 + ybar * n/sigma2) * tau_n


# so the posterior theta | y1,...,yn can be computable
  # p(theta | y1,...,yn) = dnorm(mean = mu_n, sd = sqrt(tau_n))


y <- c(9.37, 10.18, 9.16, 11.60, 10.33)
ybar <- mean(y)
n <- length(y)
#mean(y) # [1] 10.128
#var(y) # [1] 0.93047

sigma2 <- 1
tau_0 <- 100
mu_0 <- 5

tau_n <- (1/tau_0 + n/sigma2)^{-1}
sqrt(tau_n)
mu_n <- (mu_0 * 1/tau_0 + ybar * n/sigma2) * tau_n

# posterior density of theta | y1, .., yn
p(theta|y1, ..., yn) = dnorm(x, mean = mu_n, sd = sqrt(tau_n))


## Now assume we don't know the form of the posterior p(theta|y1, ..., yn)
  # so has to resort to the Metropolis algo
delta2 <- 2

Metropolis_NN <- function(theta_0, n.sim, delta2, Data) {
  set.seed(1)
  
  THETA <- NULL
  theta <- theta_0
  count <- 0
  for (i in 1:n.sim) {
    
    # sample theta_star ~ N(theta^(s), delta2)
    theta_star <- rnorm(1, mean = theta, sd = sqrt(delta2))
    
    
    # compute logr 
    logr = sum(dnorm(Data, theta_star, sqrt(sigma2), log = T) - 
      dnorm(Data, theta, sqrt(sigma2), log = T)) +
      dnorm(theta_star, mu_0, sqrt(tau_0), log = T) - 
      dnorm(theta, mu_0, sqrt(tau_0), log = T)
    
    # u ~ unif(0, 1)
    u <- runif(1)
    if(log(u) < logr) {
      # update theta using theta_star
      theta <- theta_star
      count <- count + 1
    }
  
    # update
    THETA <- rbind(THETA, theta)
  }
  AC_ratio <- count / n.sim
  return(list(theta = THETA, AC_ratio = AC_ratio))
}

Theta_pst <- Metropolis_NN(theta_0 = 0, n.sim = 1e4, Data = y)

par(mfrow = c(1,1))
plot(1:1e3, Theta_pst, type = "l")

hist(Theta_pst[10:1e4], freq = F, breaks = 1e2)
lines(density(Theta_pst), col = "red")


#------------------
# Selecting delta2: ideal AC_ratio : 20%-50%
#------------------

# delta2 = 1/32, 1/2, 2, 32, 64

delta <- c(1/32, 1/2, 2, 32, 64)
sapply(delta, function(x) Metropolis_NN(theta_0 = 0, n.sim = 1000, 
                                    delta2 = x, Data = y)$AC_ratio)


# [1] 0.816 0.579 0.349 0.098 0.068
# want AC_ratio = 20% - 50% 
  # when delta2 = 32, 64, AC_ratio is too low with more than 90% RJ
  # when delta2 = 1/32, 1/2 AC_ratio is too high with most of them AC
  # both of these two scenarios end up with highly correlated MC values
    # either correlated with theta_start or theta_{s-1}

  # so choose delta2 = 2 with AC_ratio = 35% 











