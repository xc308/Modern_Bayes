#==========================
# Three-components Mixture
#==========================

install.packages("MCMCpack")
library(MCMCpack)


#-----------
# Data set
#-----------
Mix_Data <- read.csv("Lab8Mixture.csv", header = F)
head(Mix_Data)
Mix_Data <- Mix_Data$V1

str(Mix_Data)
# 'data.frame':	500 obs. of  1 variable:

mean(Mix_Data) # [1] 2.781925
var(Dat_Mix)
# [1] 4.637941



#--------------------------------------
# Parameters for prior and hyperpriors
#--------------------------------------

## for prior epsilon^2 ~ InverseGamma(2, 2)
e.a = 2 
e.b = 2

## for hyper-priors mu0 ~ N(0, 3)
mu0.m = 0
mu0.v = 3

## for hyper-prior sigam0^2 ~ InverseGamma(2, 2)
s0.a = 2
s0.b = 2


#----------------------------------------
# Pre-Functions for posterior calculation
#----------------------------------------

NN_update <- function(prior.mu, prior.v, like.v, data) {
  # computes the posterior parameters for Normal-Normal model
  # prior.mu: prior mean
  # prior.v: prior variance
  # like.v : variance of likelihood 
  # data: observations
  
  n = length(Data)
  x = sum(Data)
  
  var = 1 / prior.v + n / like.v
  mean = (prior.mu / prior.v + x / like.v) / var
  
  return(c(mean, var))
}



NIG_update <- function(a, b, mu, Data) {
  # compute the parameters for posterior of IG(shape, rate)
  # a: prior shape
  # b: prior rate
  # mu: mean of the likelihood
  # Data: obs
  
  n = length(Data)
  a.post <- a + n/2 
  b.post <- b + 1/2 * sum((Data - mu)^2)
  
  return(c(a.post, b.post))
}



# counter for counts in each category
counter <- function(zs, k) {
  # counts the number of zi = j, j = 1, 2, ..., k
  # k: length/ number of categories
  
  counts <- vector(mode = "numeric", length = k)
  for (j in 1:k) {
    counts[j] <- sum(zs == j)
  }
  return(counts)
}



## probs for zi ~ Cat(probs)
CatProbs <- function(w, yi, mu, epsilon) {
  # calculate probability of Zi that defines Catogorical distrubtion
  # Pr(Zi == j)
  # w = (w1, w2, w3)
  numerator <- w * dnorm(yi, mean = mu, sd = sqrt(epsilon)) # W is a vector of 3
  probs <- numerator / sum(numerator)
  
  return(probs)
}


## sample from categorical distribution
CateSampler <- function(probs) {
  
  cdf <- cumsum(probs)
  x <- runif(1)
  samp <- which(x <= cdf)[1] 
  # if more than one indice corresponds to condition, just use the 1st output
  # the first indices that satisfy the condition
  return(samp)
}


#---------------------------------------
# Gibbs sampler for 3-Components Mixture
#---------------------------------------

Gibbs_MultiCmpts <- function(ys, zs, mus, epsilon, mu0, s0, w, nsim, burnin) {
  # ys : y_{1:n} data
  # zs:z_{1:n} initial value of z
  # mus: (mu_1, mu_2, mu_3) initial values of mu of 3 category
  # epsilon: like variance initial value
  # mu0: initial value for prior mean on mu
  # s0: sigma_0^2 initial value for prior var on mu
  # w = (w1, w2, w3), initial weights of mixture 
  
  k = length(mus)
  n = length(ys)
  m = length(mus) + length(w) + 3 
  
  res <- matrix(NA, nrow = n.sim, ncol = m)
  for (i in 1:n.sim) {
    # sample mu0 | ...
    mu0_pst_par <- NN_update(prior.mu = mu0.m, prior.v = mu0.v, like.v = s0, 
                             data = mus)
    mu0 <- rnorm(1, mean = mu0_pst_par[1], sd = sqrt(mu0_pst_par[2]))
    
    # sample s0|...
    s0_pst_par <- NIG_update(a = s0.a, b = s0.b, mu = mu0, Data = mus)
    s0 <- 1/rgamma(1, shape = s0_pst_par[1], rate = s0_pst_par[2])
    
    # sample epsilon |...
    eps_pst_par <- NIG_update(a = e.a, b = e.b, mu = mus[zs],Data = ys)
    epsilon <- 1/rgamma(1, shape = eps_pst_par[1], rate = eps_pst_par[2])
    
    # sample w = (w1, w2, w3)|... ~ Dir(N1+1, N2+1, N3+1)
    # get the counts of each category first
    counts <- counter(zs, k)
    w <- rdirichlet(n = 1, alpha = counts + 1)
    
    # sample mu[j]|...
    for (j in 1:k) {
      data.j <- ys[zs == j]  # zs is index of each category 1, 2, 3
      mj_pst_par <- NN_update(prior.mu = mu0, prior.v = s0, like.v = epsilon, 
                              data = data.j)
      
      mus[j] <- rnorm(1, mean = mj_pst_par[1], sd = sqrt(mj_pst_par[2]))
    }
    
    
    # sample zi|....
    # calculate the probs determine Categorical(probs)
    z.probs <- c()
    for (i in 1:n) {
      z.probs[i] <- CatProbs(w = w, yi = ys[i], mu = mus, epsilon = epsilon)
    }
    
    #probs <- vapply(1:n, function(i) {CatProbs(w = w, yi = ys[i], 
                                      #mu = mus, epsilon = epsilon)}, numeric(1))
    
    # sample from Categorical(probs)
    zs <- apply(z.probs, 2, function(x) {CateSampler(x)})
    
    
    # collect posterior samples at the ith run
    res[i, ] <- c(mu0, s0, epsilon, w, mus)
  }
  return(res[burnin:n.sim, ])
}



#----------------
# Initial values
#----------------

# zs, 
# mus, 
# epsilon, 
# mu0, 
# s0, 
# w,

w_start = rdirichlet(1, alpha = c(1, 1, 1))
dim(w_start) #[1] 1 3
zs_start = replicate(length(Mix_Data), CateSampler(w_start))

zs_start = replicate(9, CateSampler(w_start))
mus_start = rnorm(3)
mus_start[zs_start]


mu0_start = 1  # mean para of prior on mu
s0_start = 1  # var para of prior on mu
mean(Mix_Data) # [1] 2.781925
# not exactly the same as data mean, but not hugely different
  # allow better variation and mixing and wont' take very long 
    # to converge

epsilon_start = 5 # from var(Mix_Data) # [1] 4.637941


#------------------------
# starting value setting
#------------------------

w_start = rdirichlet(1, c(1, 1, 1))
zs_start = replicate(length(Mix_Data), CateSampler(probs = w_start))
# 500 starting value for z_i, one for each yi
#mus = (mu1, mu2, mu3)
mus_start = rnorm(3)
mu0_start = 1
s0_start = 1
epsilon_start = 5
#var(Mix_Data) # 4.63 


#-----------------------
# To understand the code
#-----------------------

Mix_Data[1]
P <- matrix(NA, nrow = 5, ncol = 3)
for(i in 1:5) {
  P[i, ] <- CatProbs(w = w_start, yi = Mix_Data[i], mus = mus_start,eps = 5)
}
P #5*3
#          [,1]      [,2]      [,3]
# [1,] 0.1808866 0.4063159 0.4127975

apply(P, 2, CateSampler)
# [1] 1 1 2     underestimate the zs numbers, 
  # affect data.j <- ys[zs == j]
  # less ys will be remained/selected for each category
apply(P, 1, CateSampler)
# [1] 1 3 2 3 2

Mix_Data[1:5][c(T, F)] 
# [1] 2.453943 1.965925 4.875129
# logical assertion only 2 outcomes, data has 1:5
# recycle, T F T F T








Mix_Data <- read.csv("Lab8Mixture.csv", header = F)
head(Mix_Data)
Mix_Data <- Mix_Data$V1

str(Mix_Data)
# 'data.frame':	500 obs. of  1 variable:

mean(Mix_Data) # [1] 2.781925
var(Mix_Data)
# [1] 4.637941









CateSampler
cumsum(w_start)
x <- runif(1)

which(x <= cdf)[1]
replicate(5, CateSampler(w_start))


y <- rnorm(15)
str(y)
y[zs_start]


vector <- c(3, 5, 1, 6, 12, 4)
which(vector > 5)[1]

































