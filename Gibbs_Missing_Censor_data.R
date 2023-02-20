#============================
# Truncated Gamma Inverse CDF
#============================
library(coda)


# This is a function that generate samples from tuncated Gamma

# Arg:
  # c: trunction limit, lower bound (c, infty)
  # a, b: parameters for Gamma

SampleTruncGamma <- function(c, a, b) {
  # from quantile c, find the corresponding lower bound of probability
  p0 <- pgamma(c, shape = a, rate = b)
  
  ## inverse CDF:
  # 1. sample u ~ unif
  u <- runif(1, min = p0, max = 1)
  
  # 2. inverse CDF
  z <- qgamma(u, shape = a, rate = b)
  
  return(z)
}


#==========================================
# Gibbs sampler for missing / censored data
#==========================================

# This function samples the desired theta parameter and missing/censored data 
# using Gibbs

# Arg:
  # c: a vector of values for censored data
  # z: a vector of values for obs data

  # theta.start: starting value for theta
  # zmissing.start: starting value for missing/censored data

  # n.sim : number of iteration for Gibbs sampler
  # burnin: number of samples throw away

  # r, a, b: known parameter for Gamma

Missing_Gibbs <- function(c, z, theta.start, zmiss.start, n.sim, burnin, r, a, b) {
  m = length(c) # length of missing values
  n = m + length(z) # total lenght of data: missing + obs
  
  res <- matrix(NA, nrow = n.sim, ncol = 1 + m)
  sumz <- sum(z)
  zmiss <- zmiss.start
  for (i in 1:n.sim) {
    # sample theta from Gamma(theta| nr+a, b+sum(obs+missing))
    sum_all <- sumz + sum(zmiss)
    theta <- rgamma(1, shape = n*r + a, rate = b + sum_all)
    
    # sample zmissing
    zmiss <- vapply(c, function(x) SampleTruncGamma(x, r, theta), numeric(1))
    
    # collect ith run samples, cols = 1+m
    res[i, ] <- c(theta, zmiss)
  }
  return(res[burnin:n.sim, ])
}




Missing_Gibbs2 <- function(c, z, theta.start, zmiss.start, n.sim, burnin, r, a, b,
                           Thining = TRUE, Thin = 5) {
  m = length(c)
  n = m + length(z)
  
  sumz = sum(z)
  
  res <- matrix(NA, nrow = n.sim, ncol = 1+m)
  res[1, ] <- c(theta.start, zmiss.start) 
  
  for (i in 2:n.sim) {
    # sample theta
    sum_all = sumz + sum(res[i-1, -1])
    res[i, 1] <- rgamma(1, shape = n*r+a, rate = b + sum_all)
    
    # sample zmiss
    res[i, -1] <- vapply(c, function(x) SampleTruncGamma(x, r, res[i, 1]), numeric(1))
  }
  
  if (Thining) {
    res <- res[burnin:n.sim, ]
    return(res[seq(1, nrow(res), Thin), ])
  } else {
    return(res[burnin:n.sim, ])
  }
}

# ref: thinning 
# https://stackoverflow.com/questions/23279550/select-every-nth-row-from-dataframe
# df.new = df[seq(1, nrow(df), 5), ]

#-------------
# Test on data
#-------------
#total data:
#  3.4, 2.9, 1.2+, 1.4, 3.2, 1.8, 4.6, 1.7+, 2.0+, 1.4+, 2.8, 0.6+
# z3, z8, z9, z10, z12 missing/censored

c <- c(1.2,1.7,2.0,1.4,0.6)
z <- c(3.4,2.9,1.4,3.2,1.8,4.6,2.8)

r <- 12 # twelve in total
a <- b <- 1

theta.start <- 1
zmiss.start <- rgamma(length(c), r, theta.start)
# z1, ..., zn ~ Gamma(r, theta)
# so missing z must also from this distribution

n.sim <- 1e4
burn.in <- 1000
Thin <- 5

Res2 <- Missing_Gibbs2(c, z, theta.start = theta.start, zmiss.start = zmiss.start,
               n.sim = n.sim, burnin = burn.in, r, a, b)


# with thinning 
Res3 <- Missing_Gibbs2(c, z, theta.start = theta.start, zmiss.start = zmiss.start,
                       n.sim = n.sim, burnin = burn.in, r, a, b, Thining = T, Thin = 5)


#==================
# Diagnostic plots
#==================

# trace plot and running averages for non-convergence
# lag-1 acf and N_eff for independence of samples
# 

#------------------------------------
# Trace plot for theta w/o thinning
#------------------------------------
plot(burn.in:n.sim, Res2[, 1], type = "l", 
     main = bquote("Tracplot of " ~ theta), 
     xlab = "Iterations", 
     ylab = expression(theta ~"|"~ "Oberserved Data"))

lag1.acf <- acf(Res2[, 1])$acf[2]
# [1] 0.3665532

effectiveSize(Res2[, 1])
# var1 
# 4547.507 


#------------------------------------
# Trace plot for theta w thinning
#------------------------------------
par(mar = c(4.5, 4.5, 6, 3))

pdf("Traceplot_theta_thin.pdf")
plot(burn.in:(burn.in + length(Res3[, 1]) - 1), Res3[, 1], 
     type = "l",
     xlab = paste("Iterations after burn-in: ", burn.in),
     ylab = expression(theta ~"|" * "Observed data"),
     main = bquote(atop("Trace plot of "~ theta, 
                        atop("Thinning ="~ .(Thin), 
                             atop("Lag-1 acf: " ~ .(round(acf_theta_thin, 2)),
                             "N_eff: " ~ .(effectiveSize(Res3[, 1]))))))
     )

dev.off()

acf_theta_thin <- acf(Res3[, 1])$acf[2]
# [1] 0.0146614
effectiveSize(Res3[, 1])
# 1801


## 
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 6, 3))

#pdf("Traceplot_Z3.pdf")
plot(burn.in:(burn.in + length(Res3[, 1]) - 1), Res3[, 2], 
     type = "l",
     xlab = paste("Iterations after burn-in: ", burn.in),
     ylab = expression(theta ~"|" * "Observed data"),
     main = bquote(atop("Trace plot of "~ Z_3, 
                        atop("Thinning ="~ .(Thin), 
                             atop("Lag-1 acf: " ~ .(round(acf_z3_thin, 2)),
                                  "N_eff: " ~ .(effectiveSize(Res3[, 2]))))))
)
#dev.off()

acf_z3_thin <- acf(Res3[, 2])$acf[2]


plot(burn.in:(burn.in + length(Res3[, 1]) - 1), Res3[, 3], 
     type = "l",
     xlab = paste("Iterations after burn-in: ", burn.in),
     ylab = expression(Z_8 ~"|" * "Observed data"),
     main = bquote(atop("Trace plot of "~ Z_8, 
                        atop("Thinning ="~ .(Thin), 
                             atop("Lag-1 acf: " ~ .(round(acf(Res3[, 3], plot = F)$acf[2], 2)),
                                  "N_eff: " ~ .(effectiveSize(Res3[, 3]))))))
)

acf_z8_thin <- acf(Res3[, 3], plot = F)$acf[2]

plot(burn.in:(burn.in + length(Res3[, 1]) - 1), Res3[, 4], 
     type = "l",
     xlab = paste("Iterations after burn-in: ", burn.in),
     ylab = expression(Z_9 ~"|" * "Observed data"),
     main = bquote(atop("Trace plot of "~ Z_9, 
                        atop("Thinning ="~ .(Thin), 
                             atop("Lag-1 acf: " ~ .(round(acf(Res3[, 4], plot = F)$acf[2], 2)),
                                  "N_eff: " ~ .(effectiveSize(Res3[, 4]))))))
)



#----------------
# Running average
#----------------

Running_Avg <- function(Data) {
  n <- length(Data)
  
  avg <- c()
  total <- 0
  for (i in 1:n) {
    total = total + Data[i]
    avg[i] <- total / i
  }
  plot(1:n, avg, type = "l", main = "Running Average")
}



## for theta
Res3[, 1]

par(mfrow = c(1, 1))
Running_Avg(Res3[, 1])
Running_Avg(Res[, 1])


#-----------
# Density
#----------

lines(density(Res[, 1]))

plot(density(Res[, 1]), type = "l")


#------------
# HPD interval
#------------
# 90% Credible interval
library(coda)

HPD_int <- HPDinterval(as.mcmc(Res[, 1]), prob = 0.9)
#         lower    upper
#  var1 3.450658 4.808671
# attr(,"Probability")
# [1] 0.9

abline(HPD_int[1], col = "red")

abline(v = c(HPD_int[1], HPD_int[2]), col = "red")

data.frame(HPD_int)

library(xtable)
xtable(data.frame(HPD_int))


#-------------
# HPD intervals for all parameters
#-------------
HPD_int <- matrix(NA, nrow = 6, ncol = 2)
for (i in 1:6) {
  HPD_int[i, ] <- HPDinterval(as.mcmc(Res[, 1]), prob = 0.9)
  
}

HPD_df <- data.frame(HPD_int)
rownames(HPD_df) <- c("theta", "Z3", "Z8", "Z9", "Z10", "Z12")
colnames(HPD_df) <- c("Lower", "Upper")

xtable(HPD_df)


#-------------------
# Change hyper-prior
#-------------------

# a = 1, b = 100
Res_b100 <- Missing_Gibbs2(c, z, theta.start = theta.start, 
              zmiss.start = zmiss.start, n.sim = n.sim, 
              burnin = 1, r = 12, a, 1, b = 100, Thining = F)

Res_b100[, 1]

par(mfrow = c(1, 2))
plot(density(Res[, 1]), type = "l")
abline(v = c(HPD_int[1], HPD_int[2]), col = "red")
abline(v = c(HPD_int[1], HPD_int[2]), col = "red")
abline(v = HPD_int[2], col = "red")


plot(density(Res_b100[, 1]), type = "l")
HPD_int_b100 <- HPDinterval(as.mcmc(Res_b100[, 1]), prob = 0.9)
abline(v = c(HPD_int_b100[1], HPD_int_b100[2]), col = "red")


# a = 1, b = 10000
Res_b1e4 <- Missing_Gibbs2(c, z, theta.start = theta.start, 
                           zmiss.start = zmiss.start, n.sim = n.sim, 
                           burnin = 1, r = 12, a, 1, b = 1e4, Thining = F)


plot(density(Res_b1e4[, 1]), type = "l")
HPD_int_b1e4 <- HPDinterval(as.mcmc(Res_b1e4[, 1]), prob = 0.9)
abline(v = c(HPD_int_b1e4[1], HPD_int_b1e4[2]), col = "red")


## The larger the rate b, the smaller the scale = 1/rate
  # i.e., the number of time waited betw two events get smaller
  # graph contracts


## a = 100, b = 1
Res_a100 <- Missing_Gibbs2(c, z, theta.start = theta.start, 
                           zmiss.start = zmiss.start, n.sim = n.sim, 
                           burnin = 1, r = 12, a = 100, b = 1, Thining = F)


plot(density(Res_a100[, 1]), type = "l")
# shifted mean a/b = 100 and density curve to the right
# longer waiting time for 100 events. 




