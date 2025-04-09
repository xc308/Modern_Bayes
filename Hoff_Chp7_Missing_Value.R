#===========================================
# Missing value (Pima health data)
#===========================================

# Pima.tr data set is in MASS package
install.packages("MASS")
library(MASS)

data("Pima.tr")

head(Pima.tr)
Y <- Pima.tr[2:5]
Y0 <- Pima.tr[, 2:5] #same

n <- dim(Y)[1]; p <- dim(Y)[2]


# self-construct missing value
set.seed(1)
O <- matrix(rbinom(n*p, 1, 0.9), n, p)
#rbinom(n, size, prob)
#interpreted as the number of ‘successes’ in size trials.

Y[O == 0] <- NA


#---------------------------
# pairwise relationship plot
#---------------------------
pdf("pima_pairwise_plot.pdf", family = "Times")
par(mfrow = c(p, p), mar = c(1, 1, 0.5, 0) * 1.75, #botom, left, top, right margin
    mgp = c(1.75, 0.75, 0)) # title mar, x-axis, y-ax margin

for (j1 in 1:p) {
  for (j2 in 1:p) {
    if (j1 == j2) {
      hist(Y[, j1], main = "")
      mtext(colnames(Y)[j1], side = 3, line = -0.1, cex = .7)
    } else {
      plot(Y[, j1], Y[, j2], xlab = "", ylab = "", pch = 16, cex = .7)
    }
  }
}
dev.off()


#---------------------
# parameters in prior
#---------------------

mu0 <- round(apply(Y, 2, function(x) mean(x, na.rm = T)), 0)
mu0<-c(120,64,26,26)

sd0 <- round(apply(Y, 2, function(x) sd(x, na.rm = T)), 0)
sd0<-(mu0/2)

L0 <- matrix(0.1, p, p)
diag(L0) <- 1
L0 <- L0 * outer(sd0, sd0) # var = sd*sd
#        glu    bp skin  bmi
#glu  1024.0  35.2   32 19.2
#bp     35.2 121.0   11  6.6
#skin   32.0  11.0  100  6.0
#bmi    19.2   6.6    6 36.0

eta0 <- p + 2
S0 <- L0


#-----------------
# starting values
#-----------------
Sigma_start <- S0

Y.fll <- Y
for (j in 1:p) {
  Y.fll[is.na(Y.fll[, j]), j] <- mean(Y.fll[, j], na.rm = T)
}


which(O == 0) # index where NA in Y lies

#--------------
# Gibbs sampler
#--------------
library(mvtnorm)
library(MCMCpack)


Sig <- Sigma_start
n.sim <- 1000

THETA <- NULL
SIGMA <- NULL
Y.MISS <- NULL


rsum_NA <- rowSums(Y)
array(rsum_NA)
miss_row_idx <- which(is.na(array(rsum_NA)))

for (s in 1:n.sim) {
  
  # sample theta | Sig, y1..yn ~ MVN(mu_n, L_n)
  iL0 <- solve(L0)
  iSig <- solve(Sig)
  ybar <- apply(Y.fll, 2, mean)
  
  L_n = solve(iL0 + n * iSig)
  mu_n = L_n %*% (iL0 %*% mu0 + n * iSig %*% ybar)
  
  theta <- rmvnorm(1, mean = mu_n, sigma = L_n)
  theta <- t(theta)
  
  # sample Sig | y1...yn, theta ~ InvWishart(eta_n, S_n)
  eta_n = eta0 + n
  
  S_theta = (t(Y.fll) - c(theta)) %*% t(t(Y.fll) - c(theta))
  S_n = S0 + S_theta
  
  Sig <- riwish(v = eta_n, S = S_n)
  #Sig <- solve(rwish(v = eta_n, S = solve(S_n)))

  
  # sample Yi,miss | yi,obs, theta, Sig ~ MVN(theta_ms|ob, Sig_ms|ob)
  for (i in miss_row_idx) {
    
    #str(O) int [1:200, 1:4]
    ms <- (O[i, ] == 0) # logical, only look at T
    ob <- (O[i, ] == 1) 
    
    iSig_obs <- solve(Sig[ob, ob])
    
    #str(Y.fll[1, ob]) # data.frame':	1 obs. of  3 variables
    #str(t(Y.fll[1, ob])) # num [1:3, 1] 68 28 30.2
    
    theta_msob <- theta[ms] + Sig[ms, ob] %*% iSig_obs %*% (t(Y.fll[i, ob]) - theta[ob])
    Sig_msob <- Sig[ms, ms] - Sig[ms, ob] %*% iSig_obs %*% Sig[ob, ms]
    
    Y.fll[i, ms] <- rmvnorm(1, mean = theta_msob, sigma = Sig_msob)
  }
  
  
  # collect
  THETA <- rbind(THETA, t(theta))
  SIGMA <- rbind(SIGMA, c(Sig))
  Y.MISS <- rbind(Y.MISS, Y.fll[O == 0])
  
  # print
  cat(s, theta, "\n")
}


#------
# Res
#-------
par(mfrow = c(2, 2))
str(THETA) # num [1:1000, 1:4]

for (j in 1:p) {
  hist(THETA[, j], breaks = 100, 
       main = bquote("Histogram of "~ theta[.(j)]))
}


THETA_CI <- apply(THETA, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
t(THETA_CI)
#           2.5%       50%     97.5%
#[1,] 119.54286 123.55182 127.95586
#[2,]  69.45147  71.07638  72.71508
#[3,]  27.76889  29.37676  30.98287
#[4,]  31.31825  32.17941  33.06055


#----------------------
# Y.MISS vs True value
#----------------------

Y.true <- Y0
Y.miss.true <- Y.true[O == 0]
str(Y.miss.true) #  num [1:85] 

str(Y.MISS)
# num [1:1000, 1:85]
Y.MISS.pred <- apply(Y.MISS, 2, mean)

V <- matrix(1:p, nrow = n, ncol = p, byrow = T)
v.miss <- V[O == 0] # missing for each variable 1:p
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
#[21] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3
#[41] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
#[61] 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
#[81] 4 4 4 4 4

# 1 : all the missing for variable 1
# 4: all the missing for variable 4


length(which(is.na(Y[, 1]))) # [1] 15

par(mfrow = c(2, 2))
for (j in 1:p) {
  plot(Y.miss.true[v.miss == j], Y.MISS.pred[v.miss == j],
       xlab = paste("True", colnames(Y.true)[j]),
       ylab = paste("Predict", colnames(Y.true)[j]),
       pch = 16)
  abline(a = 0, b = 1) # intercept 0, slope = 1
  cat(j, mean((Y.miss.true - Y.MISS.pred)^2), "\n")
}


Y


#-------------
# Cov to Corr 
#-------------
# often for a set of variables, 
  # correlation is more of interest than covariance

# corr = cov(x1, x2) / sqrt(var(x1), var(x2))

# transform SIGMA into correlation 
CORR <- array(dim = c(p, p, 1000))
for(i in 1:nrow(SIGMA)) {
  SIG <- matrix(SIGMA[i, ], p, p)
  DD <- outer(diag(SIG), diag(SIG))
  CORR[, , i] <- SIG / sqrt(DD)
}

# generates from 1000 posterior SIGMA to 
  # 1000 slice of correlation matrix
  # each slice is 4*4 correlation matrix 

rownames(CORR) <- colnames(CORR) <- colnames(Y)



# the posterior mean of 1000 slices of correlations
  # E[Corr|y1,...,yn]
Corr_pst <- apply(CORR, c(1, 2), mean)
#           [,1]      [,2]      [,3]      [,4]
#[1,] 1.0000000 0.2282779 0.2522819 0.1928715
#[2,] 0.2282779 1.0000000 0.2545687 0.2434800
#[3,] 0.2522819 0.2545687 1.0000000 0.6566952
#[4,] 0.1928715 0.2434800 0.6566952 1.0000000


Corr_CI <- apply(CORR, c(1, 2), function(x) quantile(x, c(0.025, 0.5, 0.975)))

Corr_CI_t <- array(dim = c(4, 3, 4))
for (s in 1:4) {
  Corr_CI_t[, , s] <- t(Corr_CI[, , s])
  colnames(Corr_CI_t[, , s]) <- c("2.5%", "50%", "97.5%")
  
}

t(Corr_CI[, , 1]) # CI for each sigma element in the 1st row of SIGMA


#------------
# Function for posterior quantile intervals for matrices
#------------
install.packages("sbgcop")
library(sbgcop)

par(mfrow = c(4, 2), mgp = c(1.75, 0.75, 0), mar = c(1, 2.75, 1, 1),
    oma = c(1.5, 0, 0, 0))
plotci.sA(CORR)
# skin and bmi are highly correlated

REG<-sR.sC(CORR)
plotci.sA(REG)




plotci.sA(Corr_CI)
plotci.sA

head(Y)





SIGMA_CI <- apply(SIGMA, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
t(SIGMA_CI)






#----------------
# Mistake points
#----------------

# 
t(Y.fll[1, ob])
str(Y.fll[1, ob]) # data.frame':	1 obs. of  3 variables
str(t(Y.fll[1, ob])) # num [1:3, 1] 68 28 30.2

# matrix subset matrix

vals <- outer(1:5, 1:5, FUN = "paste", sep = ",")
idx <- matrix(ncol = 2, byrow = T, 
              c(1, 1,
                3, 1,
                2, 4))

#matrix subset matrix
  # each row in idx decide the location of one value
  # each col in idx corresponds to one dim of array being subset
  # a 2 col matrix subset a 2 dim array(matrix)
  # a 3 col matrix subset a 3 dim array, etc. 
vals[idx]
# [1] "1,1" "3,1" "2,4"


Y.fll[O == 0] #is NOT the same as Y.fll[(O == 0)]
str(O == 0) # logi [1:200, 1:4]
str(Y.fll) # 'data.frame':	200 obs. of  4 variables:
str((O == 0))

str(Y.fll[O == 0]) # num [1:81]
str(Y.fll[(O == 0)]) # num [1:81] 

which(O == 0) # is NOT the same as Y.fll[O == 0] 





rsum_NA <- rowSums(Y)
array(rsum_NA)
miss_row_idx <- which(is.na(array(rsum_NA)))



#cond <- sapply((O[i, ] == 0), function(x) x == F)

#while (((O[i, ] == 0) == F)[1] & ((O[i, ] == 0) == F)[2] &
# ((O[i, ] == 0) == F)[3] & ((O[i, ] == 0) == F)[4] & i < n) {
#i = i + 1
#return(i)
}
