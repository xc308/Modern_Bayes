#==================
# Mixture 3-compts
#==================
# This is a function of Gibbs sampler for 3-cmpts mixture posteriors

# Arg: 
  # ys: data, y_{1:n}
  # zs: latent tag z_{1:n}, initial values 
  # mus: (mu1, mu2, mu3), initial values 
  # eps: var of likelihood, initial values
  # mu0: mean parameter of prior on mu, initial values
  # s0: var parameters of prior on mu, initial values

# return:
  # posteriors of unknowns parameters:
    # mu0, s0, eps, w, mus[j], zs[i] | Data


#-----------------
# Prior parameters
#-----------------

mu0.m = 0
mu0.v = 3

s0.a = 2
s0.b = 2

e.a = 2
e.b = 2


#-------------
# Functionals
#-------------

## Normal-Normal Model
NN_model <- function(prior.m, prior.v, like.v, Data) {
  n = length(Data)
  x = sum(Data)
  
  var = 1/ (1/prior.v + n / like.v)  # L^{-1} is Variance
  mean = (prior.m / prior.v + x / like.v) * var
  
  return(c(mean, var))
}


## Normal- InverseGamma model
NIG_model <- function(a, b, mu, Data) {
  # compute the posterior pars for IG(shape, rate)
  # a = shape parameter of prior
  # b = rate parameter of prior
  # mu: mean of the likelihood
  # Data: obs
  
  n = length(Data)
  
  a_pst = a + n / 2 
  b_pst = b + 1/2 * sum((Data - mu)^2)

  return(c(a_pst, b_pst))
}


## counts I(zi == j)
counter <- function(zs, k) {
  # counts the number when zi == j, j = 1...,k categories
  
  counts <- vector(mode = "numeric", length = k)
  for (j in 1:k) {
    # zs == 1
    # [1]  TRUE FALSE FALSE FALSE  TRUE
    counts[j] <- sum(zs == j)
  }
  return(counts)
}


## Calculate the probability of P(zi == j)
CatProbs <- function(w, yi, mus, eps) {
  num <- w * dnorm(yi, mean = mus, sd = sqrt(eps))
  z.probs <- num / sum(num)
  
  return(z.probs)
}


## sample from a categorical function
CatSampler <- function(probs) {
  cdf <- cumsum(probs)
  x <- runif(1)
  samp <- which(x <= cdf)[1]
  
  return(samp)
}



#-------------
# Gibbs Sampler
#-------------
Mix_3cmpts_Gibbs <- function(ys, zs, mus, eps, mu0, s0, w, n.sim, burnin) {
  
  n = length(ys)
  k = length(mus) # number of categories
  m = length(mus) + length(w) + 3
  
  res <- matrix(NA, nrow = n.sim, ncol = m)
  for (i in 1:n.sim) {
    # sample mu0|...., NN model
    mu0_pst_par <- NN_model(prior.m = mu0.m, prior.v = mu0.v, like.v = s0, Data = mus)
    mu0 <- rnorm(1, mean = mu0_pst_par[1], sd = sqrt(mu0_pst_par[2]))
    
    # sample s0|... NInvG model
    s0_pst_par <- NIG_model(a = s0.a, b = s0.b, mu = mu0, Data = mus)
    s0 <- 1/rgamma(1, shape = s0_pst_par[1], rate = s0_pst_par[2])
    
    # sample epsilon | ..., NInvG model
    eps_pst_par <- NIG_model(a = e.a, b = e.b, mu = mus[zs], Data = ys)
    eps <- 1/rgamma(1, shape = eps_pst_par[1], rate = eps_pst_par[2])
    
    # sample w | ...
    # get the counts N1, N2, ..., Nk
    Ns <- counter(zs, k) # a vector of k counts, one for each category
    w <- rdirichlet(1, alpha = (Ns + 1))
    
    # sample muj | ... Normal-Normal model
    for (j in 1:k) {
      data.j <- ys[zs == j] # zs == j is logical
      muj_pst_par <- NN_model(mu0, s0, like.v = eps, Data = data.j)
      mus[j] <- rnorm(1, mean = muj_pst_par[1], sd = sqrt(muj_pst_par[2]))
    }
    
    
    # sample zi: i=1:n
    # calculate the categorical probability 
    # n = length(ys)
    CatsProbs <- sapply(1:n, function(i) {CatProbs(w = w, yi = ys[i], 
                                                   mus = mus, eps = eps)}) 
    # 3 * n !! see below test on 1:5
    # Alternative 
    #CatsProbs <- vapply(1:n, function(i) {CatProbs(w = w, yi = ys[i], 
     #                                              mus = mus, eps = eps)},
      #                  numeric(3)) 
    
    
    
    # sample from categorical(CatsProbs)
    zs <- apply(CatsProbs, 2, CateSampler) 
        # n categories, one for each yi
        # margin=2 by col as CatsProbs is 3 by n
    
    
    # collect obs for this obs
    res[i, ] <- c(mu0, s0, eps, w, mus)
  }
  return(res[burnin:n.sim, ])
}


#----------------------------------------------
# Starting values settings for desired unknowns
#----------------------------------------------
mean(Mix_Data) # [1] 2.781925
mu0_start = 1
s0_start = 1

var(Mix_Data) # [1] 4.637941
eps_start = 5 # variance of data, obtained from data 

w_start = rdirichlet(1, c(1, 1, 1))
# sum(w_start) = 1 
# 1*3

mus_start = rnorm(3) # if air pollution, need truncated normal

## zs_starting
CatSampler(w_start)  # [1] 2 
# generate one category follow categorical distribution

zs_start = replicate(length(Mix_Data), CatSampler(w_start))
str(zs_start)
# int [1:500] 1 1 1 3 1 1


#-------------------------------------
# Test the Gibbs sampler on Mix_Data
#-------------------------------------

Res <- Mix_3cmpts_Gibbs(ys = Mix_Data, zs = zs_start, mus = mus_start, 
                 eps = eps_start, mu0 = mu0_start, s0 = s0_start, 
                 w = w_start, n.sim = 100, burnin = 1)


#-------
# Plots
#-------
# c(mu0, s0, eps, w, mus)
n.sim = 100
plot(1:n.sim, Res[, 1], pch = 16, cex = 0.35, 
     xlab = "Iterations", ylab = expression(mu[0]), 
     main = expression(paste("Traceplot of ", mu[0])))


plot(1:n.sim, Res[, 2], pch = 16, cex = 0.35, 
     xlab = "Iterations", ylab = expression(sigma[0]^2),
     main = expression(paste("Tracplot of ", sigma[0]^2)))


plot(1:n.sim, Res[, 3], pch = 16, cex = 0.35, 
     xlab = "Iterations", ylab = expression(epsilon^2),
     main = expression(paste("Traceplot of ", epsilon^2)))


par(mfrow = c(1, 3))
for (idx in 1:3) {
  plot(1:n.sim, Res[, 3 + idx], pch = 16, cex = 0.35, 
       xlab = "Iterations", ylab = bquote(omega[.(idx)]), 
       main = bquote(paste("Traceplot of ", omega[.(idx)])))
}

for (idx in 1:3) {
  plot(1:n.sim, Res[, 6 + idx], pch = 16, cex = 0.35, 
       xlab = "Interactions", ylab = bquote(mu[.(idx)]),
       main = bquote(paste("Traceplot of ", mu[.(idx)])))
  
}


#------------------------------------------------------
# calculate summary statistics and display using xtable
#------------------------------------------------------

## quantile for each posterior output
# Res n*9
Cred_intervals <- apply(Res, 2, function(x) quantile(x, c(0.025, 0.975)))
# [1:2, 1:9]

## Posterior mean
Pst_mean <- apply(Res, 2, function(x) mean(x))
# [1:9]

## pst mean together with credible intervals
summary_stats <- cbind(Pst_mean, t(Cred_intervals))


# c(mu0, s0, eps, w, mus)
par.names <- c("$\\mu_0$", "$\\sigma_0^2$", "$\\epsilon$",
               "$\\w_1$", "$\\w_2$", "$\\w_3$", 
               "$\\mu_1$", "$\\mu_2$", "$\\mu_3$")

rownames(summary_stats) <- par.names

head(summary_stats)

library(xtable)
xtable(summary_stats, sanitize.rownames.function = identity)
# ref: https://tex.stackexchange.com/questions/358462/including-greek-letters-when-saving-matrix-from-r-to-tex-table


#---------
# Understand structure of output of CatProbs for more than 1 data point
# for one data point, it produces a vector of 3, 
# for n data points, it produces 3 * n, matrix is stacked by col in defalut
sapply(1:5, function(i) {CatProbs(w = w_start, yi = Mix_Data[i], 
                                   mus = mus_start, eps = eps_start)}) # 

vapply(1:5, function(i) {CatProbs(w = w_start, yi = Mix_Data[i], 
                                  mus = mus_start, eps = eps_start)},
       numeric(3))

#0.02978728 + 0.87352602 + 0.09668670 = 1














