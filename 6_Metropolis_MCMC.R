#=============
# Chp 6 MCMC
#=============

# In low-dim, 
  # Important sampling and RJ sampling works pretty well

# but in high-dim
  # proposal worked in 2-D might not work in any dimension

# It's hard to capture high-dim space! 
  # turn to MCMC. 

# Markov Chain 
  # – where we go next only depends on our last state 
    # (the Markov property).
# Monte Carlo – just simulating data.

# Why MCMC?
  # (a) the region of high probability tends to be “connected”
    # get from one point to another without going through a low-probability region
  
  # (b) we tend to be interested in the expectations of functions 
    # that are relatively smooth and have lots of “symmetries”
    # That is, one only needs to evaluate them at a small number 
    # of representative points in order to get the general picture.
  

#---------------------------------  
# Advantages/Disadvantages of MCMC:
#---------------------------------  

## Advantanges:
  # applicable even when we can’t directly draw samples
  # works for complicated distributions in high-dimensional spaces,
    #even when we don’t know where the regions of high probability are
  # relatively easy to implement
  # fairly reliable


## Disadvantages:
  # slower than simple Monte Carlo or importance sampling 
    # (i.e.,requires more samples for the same level of accuracy

  # can be very difficult to assess accuracy and evaluate convergence,
    # even empirically


#---------------
# Ergodic Thrm
#---------------

# The ergodic theorem says: “if we start at a point xo and we
# keeping moving around in our high dimensional space, then we
# are guaranteed to eventually reach all points in that space with
# probability 1."


##-----------
# Metropolis Algorithm
##-----------

# Setup: Assume p.m.f.(mass) pi on Chi (countable)
          # i.e. measure sapce (Chi, sigma-algebra, pi).
        # f : Chi -> R

# Goal:
  # sample x from pi
  # apporximate E_pi[f(x)]

# supporting theory:
  # “If we take samples X = (X0, X1, . . . , ) 
    # then by the ergodic theorem,
      # they will eventually reach pi, 
        # which is known as the stationary distribution 
          # (the true pmf)."


# Approach : ergodic thrm
  # 1. If we run the Markov chain long enough, 
      # then the last state is approximately from pi.
  # 2. Under some regularity conditions,
      # 1/n * Sum_i f(x_i) --> E_pi[f(x)]



#------------------------
# Metropolis Toy Example
#------------------------

## data 
x <- c(9.37, 10.18, 9.16, 11.60, 10.33)
n <- length(x) 

## parameters for data likelihood 
Sig2 <- 1

## parameters for prior
mu0 <- 5
tau2 <- 10

## parameters for initial proposal
theta <- 0  # starting theta, to be updated
delta2 <- 1


## Metropolis sampling
T <- 500
THETA <- c()
for (t in 1:T) {
  #set.seed(1)
  
  ## sample a candidate theta_star ~ N(theta, delta2)
    # !! the 2nd run, theta needs to be updated to theta_start!!
  theta_star <- rnorm(1, theta, sqrt(delta2))

  ## calculate alpha
  log_alpha <- sum(dnorm(x, theta_star, sqrt(Sig2), log = T)) + 
    dnorm(theta_star, mu0, sqrt(tau2), log = T) -
    sum(dnorm(x, theta, sqrt(Sig2), log = T)) - 
    dnorm(theta, mu0, sqrt(tau2), log = T)
  
  ## compare with log(u), u ~ U[0, 1]
  log_u <- log(runif(1, 0, 1))
  
  if (log_u < log_alpha) { # AC
    ## !! update theta for next proposal N(theta, delta2)
                                   # = N(theta_star, delta2)
    theta <- theta_star  
  } 
  #THETA <- c(THETA, theta)
  THETA[t] <- theta
}

length(THETA)
head(THETA)
tail(THETA)

hist(THETA)
plot(1:T, THETA, type = "l")


## posterior mean from samples from Metropolis
mean(THETA) # [1] 10.06747
## posterior variance from samples from Metropolis
var(THETA) # [1] 0.2118706

## compare to the formula results



#--------------
# Toy example
#--------------
## data 
x <- c(9.37, 10.18, 9.16, 11.60, 10.33)
n <- length(x) 
mean(x)

## likelihood parameter
sig2 <- 1

## prior parameters
mu0 <- 5
tau2 <- 10  # larger 

## proposal parameter
delta2 <- 1

## initialise 
theta <- 0
THETA <- c()

## sample a candidate theta_star from symetric proposal N(0, delta2)
theta_star <- rnorm(1, 0, sqrt(delta2))

## calculate alpha
log_alpha <- sum(dnorm(x, theta_star, sqrt(sig2), log = T)) +
  dnorm(theta_star, mu0, sqrt(tau2), log = T) - 
  sum(dnorm(x, theta, sqrt(sig2), log = T)) - 
  dnorm(theta, mu0, sqrt(tau2), log = T)
  
## sample a u ~ U[0, 1]
u <- runif(1, 0, 1)

## AC condition
if (log(u) < log_alpha) {
  theta <- theta_star  # update theta for next candidate
}

## collect the theta into a vector
THETA[1] <- theta

## Remark:
prod(dnorm(x, theta, sqrt(sig2))) # [1] 6.69053e-115
# easy for numeric overflow
# so take log


## repeat above process for S times
S <- 5000
theta <- 0 ## initial
THETA <- c()
for (s in 1:S) {
  # sample a candidate theta_star ~ N(theta, sqrt(sig2))
  theta_star <- rnorm(1, theta, sqrt(delta2))
  
  # calculate alpha
  log_alpha <- sum(dnorm(x, theta_star, sqrt(sig2), log = T)) + 
    dnorm(theta_star, mu0, sqrt(tau2), log = T) - 
    sum(dnorm(x, theta, sqrt(sig2), log = T)) - 
    dnorm(theta, mu0, sqrt(tau2), log = T)

  
  # u ~ U(0, 1)
  u <- runif(1, 0, 1)
  
  # AC condition
  if (log(u) <  log_alpha) {
    theta <- theta_star ## update theta for next proposal distribution
  }
  
  AC_ratio <- mean(log(u) <  log_alpha) / (S * 0.8) # 20% thining
  
  #THETA <- c(THETA, theta)
  THETA[s] <- theta
  
  list(THETA = THETA, AC_rato = AC_ratio)
}



## plots
par(mfrow = c(1, 2))
# trace plot
plot(1:S, THETA, type = "l", xlab = "iteration s", 
     ylab = expression(theta^{(s)}))

hist(THETA, probability = T)
lines(density(THETA), col = "red")
mean(THETA)
# [1] 10.00791
var(THETA)
# [1] 0.190677

## theoretical posterior mean of theta
mu_post <- mu0 * (1/tau2) / (1/tau2 + n / sig2) + 
  mean(x) * (n / sig2) / (n/sig2 + 1/tau2)
# [1] 10.02745

## theoretical posterior variance of theta
tau2_post <- (1/tau2 + n/sig2)^{-1}
# [1] 0.1960784


#-----------------------------------------
# Sensitivity analysis of prior parameters
#-----------------------------------------

S <- 5000
Metropolis_smp <- function(mu0 = 5, tau2 = 10, delta2 = 1, plot = T) {
  set.seed(1)
  
  THETA <- c()
  AC_indx <- NULL
  for (s in 1:S) {
    # sample a candidate from proposal 1-D gaussian
    theta_star <- rnorm(1, theta, sqrt(delta2))
    
    # calculate alpha
    log_alpha <- sum(dnorm(x, theta_star, sqrt(sig2), log = T)) +
      dnorm(theta_star, mu0, sqrt(tau2), log = T) - 
      sum(dnorm(x, theta, sqrt(sig2), log = T)) - 
      dnorm(theta, mu0, sqrt(tau2), log = T)
    
    # u ~ U[0, 1]
    u <- runif(1, 0, 1)
    
    # AC condition 
    if (log(u) < log_alpha) {
      theta <- theta_star # update theta for next proposal
      
      AC_indx <- c(AC_indx, sum(log(u) < log_alpha)) # count AC times
    }
    THETA[s] <- theta
  }
  
  # AC ratio, after 15% thinning
  AC_ratio <- sum(AC_indx) / (S * 0.85) 
  
  # plot
  if (plot) {
    par(mfrow = c(1, 2))
    par(mar=c(3,3,4,3))
    # trace plot
    plot(1:S, THETA, type = "l", xlab = "iteration s",
         ylab = expression(theta^{(s)}), 
         main = 
           #expression(atop("Trace plot of "*theta)), 
         paste0("Trace plot of theta",  
               "\n Acceptance ratio = ", round(AC_ratio, 2),
               "\n mu0 = ", mu0, ", tau2 = ", tau2,
               "\n delta2 in Q(.) = ", delta2),
         cex = 0.55)
    
    # hist
    hist(THETA, probability = T, xlab = expression(theta),
         main = expression('Histogram of' ~ theta))
    lines(density(THETA), col = "red")
  }
  
  list(THETA = THETA, AC_ratio = AC_ratio)
}

try1 <- Metropolis_smp(mu0 = 5, tau2 = 10, delta2 = 1, plot = T)
try1$AC_ratio
# [1] 0.5289412

try2 <- Metropolis_smp(mu0 = 5, tau2 = 5, delta2 = 1, plot = T)
try3 <- Metropolis_smp(mu0 = 5, tau2 = 1, delta2 = 1, plot = T)
try3 <- Metropolis_smp(mu0 = 1, tau2 = 1, delta2 = 1, plot = T)

## conclusion:
  # it's irrelavent to mu0 and tau2


## now change delta2 in proposal
try4 <- Metropolis_smp(mu0 = 1, tau2 = 1, delta2 = 10, plot = T)
# detla2 = 10, AC ratio = 0.19 
  # NOT ideal!! in Normal-Normal, AC ratio = 50%

try5 <- Metropolis_smp(mu0 = 1, tau2 = 1, delta2 = 50, plot = T)
# detla2 = 50, AC ratio = 0.09 

try6 <- Metropolis_smp(mu0 = 1, tau2 = 1, delta2 = 0.5, plot = T)
# detla2 = 0.5, AC ratio = 0.64

try7 <- Metropolis_smp(mu0 = 1, tau2 = 1, delta2 = 1, plot = T)
# detla2 = 0.5, AC ratio = 0.5

## Remark:
  # AC ratio is highly connected to variance in proposal distribution


#sum(c(sum(T), sum(F)))


#------
# expression for main and lables
#------
# ref: https://stackoverflow.com/questions/4973898/combining-paste-and-expression-functions-in-plot-labels
# should be able to use expression without paste. 
# If you use the tilda (~) symbol within the expression function 
# it will assume there is a space between the characters, 
# or you could use the * symbol and 
# it won't put a space between the arguments


# line break in expression
# ref: https://stackoverflow.com/questions/18237134/line-break-in-expression

#use atop, which will break at comma and center
#expression(atop("Histogram of "*hat(mu), Bootstrap~samples*','~Allianz))

#
alpha = rnorm(1e3)
hist(alpha,cex.main=2,cex.axis=1.2,cex.lab=1.2,main=NULL )
title <- list( bquote( paste( "Histogram of " , hat(mu) ) ) ,
               bquote( paste( "Bootstrap samples, Allianz" ) ) )

mtext(do.call(expression, title ),side=3, line = c(1,-1) , cex = 2 )


hist(1:10,cex.main=2,cex.axis=1.2,cex.lab=1.2,
     main=expression(atop("Histogram of "*hat(mu), 
                          Bootstrap~samples * ',' ~Allianz)))


