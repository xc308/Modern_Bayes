#=========================================
# Multivariate Multiple Linear Regression
#     (Uni of Minnisota)
#=========================================

data(mtcars)
head(mtcars)
str(mtcars)

mtcars$cyl <- factor(mtcars$cyl)
str(mtcars)
# 'data.frame':	32 obs. of  11 variables:
#$ mpg : num  21 21 22.8 21.4 18.7 18.1 14.3 24.4 22.8 19.2 ...
#$ cyl : Factor w/ 3 levels "4","6","8": 2 2 1 2 3 2 3 1 1 2 


Y <- as.matrix(mtcars[, c("mpg", "disp", "hp", "wt")])
head(Y) # multivariate responses

mvmod <- lm(Y ~ cyl + am + carb, data = mtcars)
coef(mvmod)
#                   mpg      disp         hp         wt
#(Intercept) 25.320303 134.32487 46.5201421  2.7612069
#cyl6        -3.549419  61.84324  0.9116288  0.1957229
#cyl8        -6.904637 218.99063 87.5910956  0.7723077
#am           4.226774 -43.80256  4.4472569 -1.0254749
#carb        -1.119855   1.72629 21.2764930  0.1749132



#----------------------------------
# SSCR sum of square cross products
#----------------------------------
n <- nrow(Y)
m <- ncol(Y)
ybar <- colMeans(Y)
#    mpg      disp        hp        wt 
# 20.09062 230.72188 146.68750   3.21725 


Ybar <- matrix(ybar, nrow = n, ncol = m, byrow = T)
# repeat ybar for n rows

# crossprod(x, y) = t(x) %*% y

SSCP.T <- crossprod(Y - Ybar)

Yhat <- mvmod$fitted.values #[1:31, 1:4]
SSCP.R <- crossprod(Yhat - Ybar)

SSCP.E <- crossprod(Y - Yhat)


round(SSCP.T, 2) == round(SSCP.R + SSCP.E, 2)
#       mpg disp   hp   wt
#mpg  TRUE TRUE TRUE TRUE
#disp TRUE TRUE TRUE TRUE
#hp   TRUE TRUE TRUE TRUE
#wt   TRUE TRUE TRUE TRUE


#----------------------------
# Estimated Error Covariance
#----------------------------

n <- nrow(Y)
m <- ncol(Y)
p <- nrow(coef(mvmod)) - 1 # 4 regressors except intercept

SSCP.E <- crossprod(Y - Yhat)

Sigma_hat <- SSCP.E / (n - p - 1) # unbiased estimate of Sigma
Sigma_tilde <- SSCP.E / n # mle estimate

diag(Sigma_hat) > diag(Sigma_tilde)
# mpg disp   hp   wt 
# TRUE TRUE TRUE TRUE


mvsum <- summary(mvmod)
mvsum[[1]]
mvsum[[3]]









