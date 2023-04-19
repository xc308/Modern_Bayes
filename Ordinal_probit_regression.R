#====================================
# Hoff Chp12 Ordinal non-numeric data (Social mobility data)
#====================================

load("socmob.RData")
head(socmob)

dim(socmob)
#  dim(socmob)
# [1] 1002    7


unique(socmob$INCOME)
# [1]  NA  11   8  25 100  40  75  20  50
# [10]   9  15  12  30   5   3   0   7   1
# [19]   4  10   2   6
 
length(unique(socmob$INCOME)) # [1] 22
sort(unique(socmob$INCOME))
# [1]   0   1   2   3   4   5   6   7   8
# [10]   9  10  11  12  15  20  25  30  40
# [19]  50  75 100

yincc <- match(socmob$INC, sort(unique(socmob$INC)))

unique(socmob$DEGREE)
# [1]  1  0  3  4  2 NA
# no deg, high, associate, bacholor, graduate

yDEG <- socmob$DEGREE + 1
sort(unique(yDEG))
# [1] 1 2 3 4 5


yAGE <- socmob$AGE
sort(unique(yAGE))

sort(unique(socmob$CHILDREN))
# [1] 0 1 2 3 4 5 6 7 8
yCHILD <- socmob$CHILDREN

yPDEG <- 1*(socmob$PDEGREE > 2) # from logical to numeric
# > 2 : bacholor degree
unique(yPDEG)
1*F # [1] 0
1*T # [1] 1
head(yPDEG)

tmp <- lm(yDEG ~ yCHILD + yPDEG + yCHILD:yPDEG)
# but this model is problematic:
  # yDEG is non-numeric ordinal variable
  # such simple naivee model assumes Normal residual is violated by yDEG

plot(table(yDEG + 1) / sum(table(yDEG + 1)),
     lwd = 2, type = "h", 
     xlab = "DEG", ylab = "probability")

#"h" for ‘histogram’ like (or ‘high-density’) vertical lines

plot(table(yCHILD) / sum(table(yCHILD)),
     lwd = 2, type = "h", 
     xlab = "CHILD", ylab = "probability")


#--------------------------
# setting up for regression
#--------------------------
# design matrix X

X <-  cbind(yCHILD, yPDEG, yCHILD*yPDEG)
head(X)
#      yCHILD yPDEG  
#[1,]      3     0 0
#[2,]      3     0 0
#[3,]      1     0 0

Y <- yDEG
head(Y) # [1] 2 1 2 4 4 5


keep <- (1:length(Y))[!is.na(apply(cbind(Y, X), 1, sum))]
length(keep) [1] 959 # out of 1002
# non of X and Y has NA in each line

X <- X[keep, ] 
Y <- Y[keep]

sort(unique(Y))
# [1] 1 2 3 4 5
ranks <- match(Y, sort(unique(Y)))
# to ensure every single element in Y is indeed in sort(unique(Y))

head(Y)
head(ranks)

tail(Y)
tail(ranks)

Y == ranks # all TRUE
uranks <- sort(unique(ranks))

n <- dim(X)[1]; p <- dim(X)[2]



#---------------------------
# Ordinal probit regression
#---------------------------
## set up
set.seed(1)

# rank ()
  # some values equal (called ‘ties’), 
    # the argument ties.method determines 
      # the result at the corresponding indices.

    # The "first" method results in :
      # a permutation with increasing values at each index set of ties,

    # "last" with decreasing values. 

    # The "random" method puts these in random order 

    # the default, "average", replaces them by their mean

    # "max" and "min" replaces them by their maximum and minimum

head(rank(Y, ties.method = "random"), 12) 
# order the values of Y so as for creating corresponding z
  # such that Y is stepwise 

# Y == rank(Y, ties.method = "random") all FALSE
#head(Y, 12)

## y = g(z) = stdNorm(z)
  # z = stdNorm^{-1}(y) = qnorm(y)

# observations
y <- Y

# latent process starting values
z <- qnorm(rank(y, ties.method = "random") / (n + 1))

# regression coeff starting values
beta <- rep(0, p)

# g function
uranks # [1] 1 2 3 4 5 sample space of Y has K values, K = 5

# only K-1 can determine g function
g <- rep(NA, length(uranks) - 1) 
K <- length(uranks)


#-------------------------------
# Full conditional distributions
#-------------------------------

# gk | y1,....,yn, z1,...zn ~ N(mu_k, sigma_k)*delta(a, b)

# zi | beta, g1,...,gk-1, y1,...,yn ~ N(t(beta)Xi, 1)

# beta | z1,...,zn ~ N(mu_n, Sig_n)


#---------------------
# parameters in priors
#---------------------

#mu_k = 0 ; sigma_k = 100

mu = rep(0, K-1)
sigma = rep(100, K-1)

beta = rep(0, p)


#----------------------------------------
# starting values of Unknowns of interest
#----------------------------------------

# starting values for latent process z1,..., zn
  # y = g(z) = N(z; 0, 1)
  # Y is essentially a step-wise function against 
    # continuous latent process z e.g., effort/time spend on education
      # a certain amount of effort/time will corresponds to a certain DEG from
        # from 1 (no deg) to 5 (graduate) ordered

  # so to get z = inv_N(y; 0, 1)
  # first need to order y
    # then get the corresponding ordered z


z <- qnorm(rank(y, ties.method = "random") / (n + 1))
#z.try <- qnorm(rank(y, ties.method = "random") / n)
# qnorm(p) : p is probablity (0, 1)

y

#--------------
# Gibbs sampler
#--------------

# collector for interested unknows
g <- rep(NA, K-1)
BETA <- matrix(NA, 1000, p)
Z <- matrix(NA, 1000, n)

iXX <- solve (t(X) %*% X)  # 3 * 3
Sig <- iXX * (n / (n + 1))
ChoSig <- chol(Sig)

ac <- 0
S <- 25000

for (s in 1:S) {
  
  # update g1,...gk ~ N(mu_k, sigma_k)*delta_{ak, bk}(gk)
  
    # mu_k = 0 ; sigma_k = 100
    # mu = rep(0, K-1)
    # sigma = rep(100, K-1)
    
  for (k in 1:(K-1)) {
  
    ak = max(z[y == k])
    bk = min(z[y == (k+1)])
    
    u = runif(1, pnorm((ak - mu[k]) / sigma[k]), 
              pnorm((bk - mu[k]) / sigma[k]))
    g[k] = mu[k] + sigma[k] * qnorm(u)
  }
  
  
  # update beta ~ MVN(mu_n, Sig_n)
  mu_n = Sig %*% (t(X) %*% z)
  beta = ChoSig %*% rnorm(p) + mu_n   # 3*1
  
  
  # update zi | beta, g1...gk-1, y1,...,yn ~ N_{a, b}(t(beta) %*% Xi, 1)
    mu_z = X %*% beta # n * 1
    
    # c(-Inf,g)
    # [1]       -Inf -1.2422886
    # [3]  0.3372440  0.5227847
    # [5]  1.2736327
    
    # [match( y-1, 0:K)] : n elements of y that matich 0:K
    a <- c(-Inf,g)[match(y-1, 0:K)] # n*1
    b <- c(g, Inf)[y]
    
    u <- runif(n, pnorm(a - mu_z), pnorm(b - mu_z))
    z <- mu_z + qnorm(u) # n*1
    
    
 # To help mixing, using Metropolis to update z and g
    # propose z_star, and g_star
    c <- rnorm(1, 0, n^(-1/3))
    z_star = z + c
    g_star = g + c
    
    # logr
    logr = sum(dnorm(z_star, mu_z, 1, log = T) - dnorm(z, mu_z, 1, log = T)) + 
      sum(dnorm(g_star, mu, sigma, log = T) - dnorm(g, mu, sigma, log = T))
  
    # AC 
    if (log(runif(1)) < logr) {
      z <- z_star
      g <- g_star
      
      ac <- ac + 1
    }
  
  
    # collect every 25th iteratioin
    if (s %% (S/1000) == 0) {
      cat(s/S, ac/s, "\n" )
      
      BETA[s/(S/1000), ] <- beta
      Z[s/(S/1000), ] <- z
    }
}
































