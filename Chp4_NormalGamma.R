#===================
# Chp 4 NormalGamma
#===================

install.packages("dplyr")
install.packages("reshape")
library(dplyr)
library(reshape)

library(ggplot2)


# case study: # (Have fun!)
  # Do teachers' expectation affect students performance?

# At the beginning of the year,
  # all students were given an IQ test

# For each class, the researchers randomly selected around
# 20% of the students, and told the teacher that these students were
# “spurters’ ’ that could be expected to perform particularly well that
# year.

# At the end of the year, 
  # all students were given another IQ test.

# The change in IQ score for the first-grade students was:

#spurters
x <- c(18, 40, 15, 17, 20, 44, 38)
length(x) 
# [1] 7
#control
y <- c(-4, 0, -19, 24, 19, 10, 5, 10,
       29, 13, -9, -8, 20, -1, 12, 21,
       -7, 14, 13, 20, 11, 16, 15, 27,
       23, 36, -33, 34, 13, 11, -19, 21,
       6, 25, 30,22, -28, 15, 26, -1, -2,
       43, 23, 22, 25, 16, 10, 29)
length(y)
# [1] 48

iqChange <- data.frame(Treatment = c(rep("spurter", length(x)), 
                         rep("control", length(y))),
           Gain = c(x, y))

head(iqChange)
tail(iqChange)

## An initial exploratory analysis
# Plot the number of students versus the change in IQ score for
# the two groups.

# How strongly does this data support the hypothesis
# that the teachers’ expectations caused the spurters to perform
# better than their classmates?

xLimits = seq(min(iqChange$Gain) - (min(iqChange$Gain) %% 5),
              max(iqChange$Gain) - (max(iqChange$Gain) %% 5),
              by = 5)

ggplot(data = iqChange, 
                aes(x = Gain, 
                       fill = Treatment, 
                       colour = I("black")))+ 
  geom_histogram(position = "dodge", alpha = 0.5,
               breaks = xLimits, closed = "left")+ 
  scale_x_continuous(breaks = xLimits,
                     expand = c(0,0))+
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 10, by = 1))+
  ggtitle("Histogram of Change in IQ Scores") + 
  labs(x = "Change in IQ Score", fill = "Group") + 
  theme(plot.title = element_text(hjust = 0.5))


#--------------------------------------------
## generate from rd normalgamma distribution
#--------------------------------------------
# ref: https://rdrr.io/rforge/nclbayes/
install.packages("nclbayes", repos="http://R-Forge.R-project.org")
library(nclbayes)

rnoga <- rnormgamma(1, 1, 2, 5, 3)
str(rnormgamma(1, 1, 2, 5, 3))
# num [1, 1:2] 0.706 1.12

rnoga[, 1]

rnoga_x <- rnormgamma(1e2, 1, 2, 5, 3)
rnoga_y <- rnormgamma(1e2, 1, 2.5, 5.5, 3)
  
cmp_xy <- rnoga_x[, 1] > rnoga_y[, 1] 
str(cmp_xy)

head(cmp_xy)

## to sum up the T count in logical vector
# ref:https://www.statology.org/r-count-true-values/#:~:text=Method%201%3A%20Use%20sum()%20sum(x%2C%20na.&text=This%20method%20will%20return%20the%20count%20of%20TRUE%20values%20in%20a%20vector.&text=This%20method%20will%20return%20the%20count%20of%20TRUE%2C%20FALSE%2C%20and,NA%20values%20in%20a%20vector.

sum_T <- sum(cmp_xy, na.rm = T)

## to calculate the probability that mu_x > mu_y
sum_T / 100
# [1] 0.6
# 60% out of 100 students 


##-----------------------------------
# sample from posterior distribution
##-----------------------------------

# Likihood is gaussian for X_s, X_c
  # X_s ~ N(mu_s, lambda_s^(-1))
  # X_c ~ N(mu_c, lambda_c^(-1))

# prior for two unknowns jointly is NormalGamma
  # (mu_s, lambda_s) ~ 
    # NormalGamma(mu_s, lambda_s | m, c, a, b)

  # (mu_c, lambda_c) ~ 
    # NormalGamma(mu_c, lambda_c | m, c, a, b)

# posterior for two unknow parameters is also NormalGamma
  # (mu_s, lambda_s | X_1:ns) ~ 
    # NormalGamma(mu_s, lambda_s | M, C, A, B)

  # (mu_c, lambda_c | X_1:nc) ~ 
    # NormalGamma(mu_c, lambda_c | M, C, A, B)

## hyperprior settings for m, c, a, b
  # m_s = m_c = 0, 
    # priorly don't know average for the mean of IQ difference

  # c = 1, don't know how big the mean of IQ difference will variate

  # a = 1/2, don't know how big the std of IQ difference will be
    # so lambda^{1/2} in NomalGamma prior will disappear

  # b  = 10^2a, the expected std (lambda) of IQ difference will be 10
    # E(lamda^{-1/2}) = {E(lamda)}^{-1/2} = a / b 
    # 10^ ^{-1/2} = a / b
    # b = 10^{1/2}a


## By derivation, 
  # parameters in posterior can be expressed in terms of hyperpriors

# Spruters x
n_x = length(x)

a = 1/2
c = 1
m = 0
b = 10^2 * a

sum_x = sum(x)
sum_x2 = sum(x^2)

A_x = a + n_x / 2
C_x = c + n_x
M_x = (c*m + sum_x) / C
B_x = b + (c*m^2 - C*M^2 + sum_x2)/2


# control y
n_y = length(y)
sum_y = sum(y)
sum_y2 = sum(y^2)

A_y = a + n_y / 2
C_y = c + n_y
M_y = (c*m + sum_y) / C_y
B_y = b + (c*m^2 - C_y*M_y^2 + sum_y2) /2


## generate 1e6 samples from NomalGamma with above
  # parameters for posteriors for two groups

# rd samples will be in the form of 
  # mu_i tau_i (lambda precision)

smpls_x = rnormgamma(1e6, b = M_x, c= C_x, g = A_x, h = B_x)
post_mu_x = smpls_x[, 1]

smpls_y = rnormgamma(1e6, b = M_y, c = C_y, g = A_y, h = B_y)
post_mu_y = smpls_y[, 1]


#------------------------------------
# calculate the P(mu_x > mu_y | data)
#------------------------------------

cmpare = post_mu_x > post_mu_y
count_T = sum(cmpare, na.rm = T)

Prob_T = count_T / 1e6
# [1] 0.970722

# conclusion: almost sure that teachers' expectations 
  # affect children's IQ change. 


##--------------------------------------------
# Alternative method for posterior derivation
##--------------------------------------------
# ref: https://github.com/resteorts/modern-bayes/blob/master/TA-lab-material/Lab%2004/problem-statement/lab-04.pdf

a_pos = a + n / 2
c_pos = c + n
m_pos = c * m + n * mean(data) / (c + n)
b_pos = b + 1/2 * sum((data - mean(data))^2) + n * c * (mean(data) - m)/(2 * (c + n))

## corresponding code

prior = data.frame(m = 0, c = 1, a = 1/2, b = 10^2*a)
findPostPars = function(prior, data) {
  PostPar = NULL
  a = prior$a
  c = prior$c
  m = prior$m
  b = prior$b
  n = length(data)
  
  PostPar = data.frame(
    m = (c * m + n * mean(data)) / (c + n),
    c = c + n,
    a = a + n / 2,
    b = b + 1/2*sum((data - mean(data))^2) + n*c*(mean(data) - m)^2 / (2*(c + n))
  )
  
  return(PostPar)
}

PostParS <- findPostPars(prior, data = x)
PostParC <- findPostPars(prior, data = y)

rbind(`Prior Pars` = prior, 
      `Spurters Post Pars` = round(PostParS, 2),
      `Control Post Pars` = round(PostParC, 2))

xtable::xtable(rbind(`Prior Pars` = prior, 
                     `Spurters Post Pars` = round(PostParS, 2),
                     `Control Post Pars` = round(PostParC, 2)),
               caption = "Prior parameters and posterior parameter")



#----------------------------------
# Sample from two posterior group
#----------------------------------
# number of posterior simulations
#sim = 1000 for plot
sim = 1e6 # for next section prob calculation

# initialize vectors to store samples
mu_s = NULL
lambda_s = NULL
mu_c = NULL
lambda_c = NULL

# first sample from gamma to get post samples for lamda_s and lamda_c
set.seed(03-02-2023)

lambda_s = rgamma(n = sim, shape = PostParS$a, rate = PostParS$b)
lambda_c = rgamma(sim, shape = PostParC$a, rate = PostParC$b)

# with the post samples for lambda
# sample posterior mean given lambda 
  # mu | lambda ~ N(mu, (c*lambda)^{-1})

mu_s = sapply(sqrt((c*lambda_s)^{-1}), rnorm, n = 1, mean = PostParS$m)
# rnorm(n, mean, sd )
mu_c = sapply(sqrt((c*lambda_c)^{-1}), rnorm, n = 1, mean = PostParC$m)


## store simulations
simDF = data.frame(lambda = c(lambda_s, lambda_c),
           mu = c(mu_s, mu_c),
           Treatment = rep(c('Spurters', 'Controls'),
                           each = sim)
           )

head(simDF)
# 3 vars : lambdda, mu, treatment
  # 2000 rows, 
    #the 1st 1000 is spurters
    # the 2nd 1000 is controls


simDF$stdev <- simDF$lambda^{-1/2}
# -1: change to variance
# sqrt: change to std deviation


#---------------------
# plot the simulations
#---------------------

ggplot(simDF, aes(x = mu, y = stdev,
                  colour = Treatment,
                  shape = Treatment)) +
  geom_point(alpha = 0.2) + 
  labs(x = expression(paste(mu, "(mean change in IQ)")),
       y = expression(paste(lambda^{-1/2}, "(std. dev. of change in IQ)"))) + 
  ggtitle("Posterior samples") + 
  theme_bw()

# for a given mean of the IQ change, 
  # spurters group show larger variations
  # while control groups show much less variations

# for a given std.dev level of IQ change
  # spurters group's mean range is more clustered in  
    # postive change from (0, 50) change in IQ
  # control group's mean range is more spread out
    # ranging from (-50, 50), means there are students
      # whose IQ decreases after one year of study
  

#---------------------------------------
# compare posterior mean between two groups
#---------------------------------------
# calculate P(mu_s > mu_c | data) = ?

sum(mu_s > mu_c) / 1e6
# [1] 0.703124


#---------------------
# prior assumptions
#---------------------

#There are a few ways you can check that the prior conforms with our prior beliefs. 
# Let's go back and checks these. 
# Draw some samples from the prior and look at them---
# this is probably the best general strategy. 
# See Figure \ref{figure:pygmalion-prior}. 
# It's also a good idea to look at sample hypothetical datasets 
# $X_{1:n}$ drawn using these sampled parameter values. 


# sample from gamma distribution for a = 1/2, b = 10^2*a
# then sample from mu |gam ~ N(m = 0, (c*lambda)^{-1})
# with c = 1

lambda_prior = rgamma(1000, shape = 1/2, rate = 10^2*a)
mu_prior = sapply(lambda_prior^{-1/2}, rnorm, n = 1, mean = 0)

sim_Prior_DF = data.frame(Lambda_prior = lambda_prior, 
           mu_prior = mu_prior)

ggplot(sim_Prior_DF, aes(x = mu_prior, y = lambda_prior)) +
  geom_point(alpha = 0.9) + 
  labs(x = expression(paste(mu, "(prior mean of IQ change)")),
       y = expression(paste(lambda^{-1/2}, "(std.dev. of IQ change)")))+
  ggtitle('samples from prior distribution') + 
  theme_bw()



