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

a.post.mean.thry <- (n + 1) / (b * sum.x + 1)




