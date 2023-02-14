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
PostPars <- findPostPriors(prior = priors, data = x)
n.sims <- 1e4

# sample from posterior distribution of (mu, lambda | data)
  # 1. sample lambda ~ gamma(A, B)
  # 2. given lambda, sample mu|lambda ~ N(M, (C*lambda)^{-1})

lambda_post = rgamma(n.sims, shape = PostPars$A, rate = PostPars$B)
mu_post = sapply(sqrt((PostPars$C*lambda_post)^{-1}), rnorm, n = 1, mean = PostPars$M)

plot(mu_post, lambda_post)









