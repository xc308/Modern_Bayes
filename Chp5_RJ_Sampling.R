#==========================
# chp 5 Monte Carlo Method
#==========================

# Q1:
# f(x) propto sin(pi * x), x in [0, 1]
# Plot the densities of f (x) and the Unif(0,1) 
  # on the same plot


# density function of f(x) =  sin(pi * x))^2
f <- function(x) (sin(pi * x))^2

# density function of Unif(0, 1)
Unif <- Vectorize(function(x) 1/(1-0))

x <- seq(0, 1, length.out = 200)
plot(x, f(x), type = "l", lwd = 2, col = 1, ylab = "")
lines(x, Unif(x), col = 2, lwd = 2)

legend("topleft", legend = c("f", "unif"),
       col = c(1, 2), lwd = 2, lty = 1, cex = 0.7)



# Q2: 
# According to the rejection sampling approach sample from f(x)
# using the Unif(0,1) pdf as an enveloping function.

# 1. implement rejection sampling for a single data point

# Initialize a vector Samp for future sample collection
Samp <- NULL
while(is.null(Samp)) {
  # as long as nothing in Samp vector, start working
  
  # step 1, sample from Unif(0,1) pdf
  x <- runif(1, min = 0, max = 1)
  y <- runif(1, min = 0, max = Unif(x))
  
  # step 2, check condition and decide whether to accept x or RJ
    # f(x) = (sin(pi * 2))^2
  if (y < f(x)) 
    Samp = x
}
x
# [1] 0.6904894


# 2. a more general RJ sampling
RJ_sampling <- function(f, rq, q) {
  while(T) { # infinite loop
    # step 1 sample from proposal
    x <- rq(1)
    y <- runif(1, min = 0, max = q(x))
    
    # step 2 check acceptance condition
    if (y < f(x)) 
      Samp <- x
    
    return(Samp)
  }
}

RJ_sampling(f, rq = runif, q = Unif)
# [1] 0.6904894



## Q3:
# Plot a histogram of the points that fall in the acceptance region.
# Do this for a simulation size of 10^2 and 10^5
# and report your acceptance ratio
# Compare the ratios and histograms

# when do points can fall into an acceptance region? 
k <- 10^2

x <- runif(k)
y <- runif(k)

samples <- x[y < f(x)] # points can fall into an acceptance region 

hist(samples)
mean(y < f(x))
# [1] 0.53 # acceptance ratio


samples_acratio <- function(k) {
  x <- runif(k)
  #y <- runif(k)
  y <- runif(k, min = 0, max = Unif(x))
  
  samples <- x[y < f(x)]
  ratio <- mean(y < f(x)) # ac ratio
  
  output <- list(samples = samples, AC_ratio = ratio)
  
}

k <- 100
outputs <- samples_acratio(k)
outputs$AC_ratio
# [1] 0.49
hist(outputs$samples)


k2 <- 1e5
outputs2 <- samples_acratio(k2)
outputs2$AC_ratio
# [1] 0.50217
hist(outputs2$samples)


k3 <- 1e6
outputs3 <- samples_acratio(k3)
outputs3$AC_ratio
# [1] 0.500278
hist(outputs3$samples)



#-----------------------------------------
# Repeat above for envelope function Beta(2,2)
#-----------------------------------------

# Q1:
# Plot the density of f(x) = (sin(pi * x))^2 
  # and Beta(2,2) in the same graph

# Q2:
# use RJ sampling to sample from f(x) using Beta(2, 2) envelope

# Q3:
# compare samples size 100 and 1e5, 
# plot histgram of points that fall into the acc region
# report acceptance ratio


## Answer for Q1
## density of new proposal Beta(2, 2)
a <- 2
b <- 2

Beta_dens <- function(x, a, b) {
  num <- x^{a - 1} * (1 - x)^{b - 1}
  B <- (gamma(a) * gamma(b)) / gamma(a + b)
  
  return(num / B)
  
}

x <- seq(0, 1, length.out = 200)
plot(x, f(x), type = "l", lwd = 2, col = 1, ylim = c(0,2))
lines(x, Beta_dens(x, a, b), lwd = 2, col = 2)

# max(Beta_dens(x, a, b))
# [1] 1.499962 for ylimit

legend("topright", legend = c(expression((sin(pi*x))^2), "Beta(2,2)"), 
       col = c(1, 2), cex = 0.7, lwd = 2, lty = 1 )


## Answer for Q2:
# for 1 sample point
# initialise a vector to collect Sample
Sample <- NULL
while(T) {
  x <- rbeta(1, shape1 = a, shape2 = b)
  y <- runif(1, min = 0, max = Beta_dens(x, a, b))
  
  if (y < f(x))
    Sample <- x
  
  return(Sample)
}
Sample
# [1] 0.6209267


# or using general RJ_sampling function, but need modification

Sample <- NULL
RJ_Samp <- function(f, rq, q, ...) {
  x <- rq(1, ...)
  y <- runif(1, min = 0, max = q(x, ...))
  
  if (y < f(x))
    Sample <- x
  
  return(Sample)
}

RJ_Samp(f, rq = rbeta, Beta_dens, 2, 2)
# [1] 0.591476
# [1] 0.6805201
# [1] 0.4321957


## Answer for Q3: 
# modify above function to involve simulation size k
Samples <- NULL
Smp_Ratio <- function(k, f, rq, q, ...) {
  x <- rq(k, ...)
  y <- runif(k, min = 0, max = q(x, ...))
  
  #if(y < f(x))
     #Samples <- x
  Samples <- x[y < f(x)] # accepted x are samples
  ratio <- mean(y < f(x)) # how many times condition is ture
  
  list(Samples = Samples, ratio = ratio)
  
}

k1 <- 1e2
outputs1_Beta <- Smp_Ratio(k = k1, f, rq = rbeta, q = Beta_dens, 2, 2)
outputs1_Beta$ratio
# [1] 0.53 [1] 0.49
hist(outputs1_Beta$Samples)

k2 <- 1e6
outputs2_Beta <- Smp_Ratio(k = k2, f, rq = rbeta, q = Beta_dens, 2, 2)
outputs2_Beta$ratio
# [1] 0.500578
hist(outputs2_Beta$Samples)

