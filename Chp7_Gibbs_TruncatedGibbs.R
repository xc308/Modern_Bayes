#===============
# Gibbs sampler
#===============
# Function of Gibbs sampler

# Args
  # start.a: initial value for a in posterior distribution;
  # start.b: initial value for b in posterior distribution;
  # n.sims: number of simulations for Gibbs sampler;
  # data: observed data, in data frame with one column variable;

# return
  # a matrix with two colums: 
    # one col is for samples of a, one for samples of b;
    # rows of the matrix is the n.sims;

Gibbs_smp <- function(start.a, start.b, n.sims, data) {
  sum.x <- sum(data)
  n <- nrow(data)
  
  # initialise a null matrix for results, allocate memory
  res <- matrix(NA, nrow = n.sims, ncol = 2)
  res[1, ] <- c(start.a, start.b)
  
  for (i in 2:n.sims) {
    res[i, 1] <- rgamma(1, shape = n + 1, rate = res[i-1, 2] * sum.x + 1)
    res[i, 2] <- rgamma(1, shape = n + 1, rate = res[i, 1] * sum.x + 1)
  }
  return(res)
}

# test on data
data <- read.csv("data-exponential.csv", header = F)

library(MASS)
Data <- read.csv("data-exponential.csv", header = FALSE)

n.sims <- 1e4
start.a <- 1
start.b <- 1
res <- Gibbs_smp(start.a = .25, start.b = .25, n.sims = 1e4, data = Data)

head(res)

# plot
hist(res[, 1], probability = T, main = "Histogram of posterior of a")
lines(density(res[, 1]), col = "red")

hist(res[, 2], probability = T, main = "Histogram of posterior of b")
lines(density(res[, 2]), col = "red")

# trace plot
#plot(1:n.sims, res[, 1], type = "l", main = "Trace plot of a")
#plot(1:n.sims, res[, 2], type = "l", main = "Trace plot of b")

# thining 5000
plot(5001:n.sims, res[5001:n.sims, 1], type = "l", main = "Trace plot of a")
plot(5001:n.sims, res[5001:n.sims, 2], type = "l", main = "Trace plot of b")


## sample posterior mean and var for a and b
a.post.mean <- mean(res[, 1]) # [1] 1.020001
a.post.var <- var(res[, 1]) # [1] 0.7447557

b.post.mean <- mean(res[, 2]) # [1] 0.8628327
b.post.var <- var(res[, 2]) # [1] 0.4692718

# theoretical posterior mean for a and b
# E[a|x_1:n] = shape / rate = (n + 1) / (b * sum.x + 1)
# E[b|x_1:n] = shape / rate = (n + 1) / (a * sum.x + 1)



#-------------------------
# Truncated Gibbs Sampler 
#-------------------------
# A toy: 
  # pure sample from conditional distribution 
  # not posterior distribution that will involve data


# p(x, y) propto exp^{-xy} I(x, y \in (0, c))

# want to sample from a truncated full conditional distri

# p(x|y) propto_x p(x, y) 
          # propto exp(-xy) I(0< x < c)
          # propto exp(x|y) I(x < c)

# p(x|y) is a truncated version of Exp(y) distribution
# is equivalent to 
  # x ~ Exp(y)
  # then conditioning on x, x<c

# refer to this truncated Exp(y) as TExp(y, (0, c))

# to sample from this truncated TExp(y, (0, c))
  # 1. sample u ~ Unif(0, F(c|theta = y))
      # F(c|theta = y) = 1 - exp^{-theta*c}

  # 2. sample z = F_inv(u|theta = y)
              # = - 1/theta * log(1-u)


## Gibbs sampler for truncated distribution

  # initialise x0, y0
  # sample x1 ~ TExp(y0, (0, c)), i.e., 
    # 1. sample u ~ Unif(0, 1 - exp^{-theta = y0*c})
    # 2. calculate z = - 1/theta = y0 * log(1-u)

  # sample y1 ~ TExp(x1, (0, c)), i.e.,
    # 1. sample u ~ Unif(0, 1 - exp^{-theta = x1*c})
    # 2. calculate z = - 1/theta = x1 * log(1-u)   


## args:
  # x0, y0 : initial values
  # c : turncate limit 
  # n.sims


## return:
  # a matrix: 2 colums, n.sims rows
    # the 1st col: samples for x in p(x|y)
    # the 2nd col: samples for y in p(x|y)

Truncat_Gibbs_smp <- function(start.x, start.y, C, n.sims) {
  Res <- matrix(NA, nrow = n.sims, ncol = 2)
  Res[1, ] <- c(start.x, start.y)
  
  for (i in 2:n.sims) {
    # sample xi from TExp(y_i-1, (0, c))
    CDF = 1 - exp(- C * Res[i-1, 2])
    u = runif(1, min = 0, max = CDF)
    z = - 1/Res[i-1, 2] * log(1 - u)
    Res[i, 1] = z
    
    # sample yi from TExp(xi, (0, c))
    CDF = 1 - exp(- C * Res[i, 1])
    u = runif(1, 0, CDF)
    z = - 1/Res[i, 1] * log(1 - u)
    Res[i, 2] = z
  }
  return(Res)
}

Trunc_toy <- Truncat_Gibbs_smp(start.x = 1, start.y = 1, C = 2, n.sims = 1e4)

par(mfrow = c(1, 1))
plot(Trunc_toy, type = "p", cex = 0.025, col = "blue")



      
















