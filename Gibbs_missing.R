#=================
# Truncated Gamma
#=================
pgamma()
qgamma()

# This function samples from a truncated gamma
# Truncation (c, infinity)
# input: c, shape = a, rate = b
# output: a sample from truncated Gamma I(Zi > ci)Gamma(Zi|a, b)

SampleTruncGamma <- function(c, a, b) {
  # feed in truncation quantile to get the corresponding 
    # lower bound probability
  p0 <- pgamma(c, shape = a, rate = b)
  
  ## iverse CDF:
  
  # 1. generate a random unif 
  x <- runif(1, min = p0, max = 1)
  
  # 2. feed in the probability corresponds to truncated quantile
  zi <- qgamma(x, shape = a, rate = b)
  
  return(zi)
}


#=========================================
# Gibbs sampler for missing (censored) data
#=========================================

# Arg:
  # c: a vector of values for censored data (missing)
  # z: a vector of values for observed data

  # theta.start: starting value of desired theta
  # zmiss.start: starting value of desired missing data

  # n.sim: # of simulations
  # burnin : # of iterations of burnin

  # r: shape of Gamma, known
  # a: shape of Gamma, known
  # b: rate of Gamma, known


Missing_Gibbs <- function(c, z, theta.start, zmiss.start, n.sim, burnin, r, a, b) {
  m = length(c) # number of missing / censored
  n = m + length(z) # total length of missing and obs
  
  res <- matrix(NA, nrow = n.sim, ncol = 1 + m) # 1 for theta, m for missing
  sumz = sum(z)
  zmiss <- zmiss.start
  for (i in 1:n.sim) {
    # sample theta from Gamma(theta| nr + a, b + sum(miss+obs))
    
    sum_all <- sumz + sum(zmiss) # 1st run use start, the following runs use sampled zmiss
    theta <- rgamma(1, shape = n*r + a, rate = b + sum_all)

    # sample desired missing 
    zmiss <- vapply(c, function(x) SampleTruncGamma(x, a = r, b = theta), numeric(1))
            # each apply generate 1 zmiss,
            # a vector of c will generate m zmissings
    
    res[i, ] <- c(theta, zmiss) # ncol = 1+m
    
  }
  return(res[burnin:n.sim, ])
}



#--------------
# Test on data
#--------------
z <- c(3.4,2.9,1.4,3.2,1.8,4.6,2.8)
c <- c(1.2,1.7,2.0,1.4,0.6)

r <- 12
a <- 1
b <- 1
theta.start <- 1
zmiss.start <- rgamma(length(c), shape = r, rate = theta.start)

Res <- Missing_Gibbs(c, z, theta.start = theta.start, zmiss.start = zmiss.start,
                     n.sim = 100, burnin = 1, r, a, b)



#-------------------
# Diagnostics Plots
#-------------------





