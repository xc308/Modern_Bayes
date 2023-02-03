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












