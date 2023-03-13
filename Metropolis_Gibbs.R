#=======================================
# Hoff 10.5 Combing Metropolis and Gibbs
#=======================================

load("icecore.RData")

head(icecore)
#     year   co2  tmp
#1 -419095 277.6 0.94
#2 -416872 286.9 1.36

# data include 200 values of temperature measured at roughly equal time intervals
# time between consecutive measurements being approximately 2,000 years
# For each value of temperature there is a CO2 concentration value 
    # corresponding to a date that
    # is roughly 1,000 years previous to the temperature value, on average.
# Temperature is recorded in terms of its difference from 
    # current present temperature in degrees Celsius,

# CO2 concentration is recorded in parts per million by volume.


max((icecore[, 3] - mean(icecore[, 3])) / sd(icecore[, 3])) # [1] 2.867939
min((icecore[, 3] - mean(icecore[, 3])) / sd(icecore[, 3])) # 
# [1] -1.664384

# less variation than on original scale
min(icecore[, 3])
# [1] -9.39

Mx <- c()
Mn <- c()
for (i in 2:3) {
  Mx[i-1]<- max((icecore[, i] - mean(icecore[, i])) / sd(icecore[, i]))
  Mn[i-1] <- min((icecore[, i] - mean(icecore[, i])) / sd(icecore[, i]))
}

Mx # [1] 2.089509 2.867939
Mn # [1] -1.604304 -1.664384



# mgp
# The margin line (in mex units) for the axis title, 
# axis labels and axis line. Note that mgp[1] affects title 
# whereas mgp[2:3] affect axis. The default is c(3, 1, 0).

pdf("fig10_7.pdf",family="Times") 
# font change from "Helctive to Times"

par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
layout(matrix( c(1,1,2),nrow=1,ncol=3))
layout(matrix( c(1,1,2),nrow=1,ncol=3))
# a matrix object specifying the location of the next NN figures on the output device
# Each value in the matrix must be 0 or a positive integer.
# If N is the largest positive integer in the matrix,
# then the integers {1,…,N−1} must also appear at least once in the matrix.

par(mfrow = c(1, 2))

plot(icecore[,1], (icecore[, 3] - mean(icecore[, 3])) / sd(icecore[, 3]),
     type = "l", col = "black", ylim = c(-2.5, 3), 
     xlab = "year", ylab = "standardized measurement")

lines(icecore[, 1], (icecore[, 2] - mean(icecore[, 2])) / sd(icecore[, 2]),
      type = "l", col = "grey")

legend("topright", legend = c("temp", expression(CO[2])), 
       bty = "n", lwd = c(2, 2), col = c("black", "grey"))

plot(icecore[, 2], icecore[, 3], 
     xlab = expression(paste(CO[2], " (ppmv)")), 
     ylab = "temperature difference (deg C)")

## Remark:
  # plot indicates that the temporal history of 
    # temperature and CO2 follow very similar patterns.
  # 2nd plot indicates that CO2 concentration at a given 
    # time point is predictive of temperature 
      # following that time point.
  # One way to quantify this is by fitting a linear regression model 
    # for temperature (Y ) as a function of CO2 (x).

# Ordinary least squares regression
  # gives an estimated model of ˆE[Y |x] = −23.02+0.08x with a nominal standard
  # error of 0.0038 for the slope term.

# The validity of this standard error relies
    # on the error terms in the regression model being 
    # independent and identically distributed, and
    # standard confidence intervals further rely on the errors being
    # normally distributed.




## Residual correlation visualisation
lmfit <- lm(icecore$tmp ~ icecore$co2)
hist(lmfit$residuals)
acf(lmfit$residuals, xlab = "lag", ci.col = "gray")
acf(lmfit$residuals, plot = F)$acf[2]
# [1] 0.5165516


summary(lmfit)
# Coefficients:
#           Estimate Std. Error t value
#(Intercept) -23.024140   0.879543  -26.18
#icecore$co2   0.079853   0.003834   20.83

# the std error is nominal as it dependes on the assumption
  # the the error being iid, 
# and the CI dependends on the error being normal

# the 1st plot indicate no obvious deviation from normality
# the 2nd plot shows non-trivial auto-correlation 0.52 betw
  # two residuals at time consecutives

# such postive correlation generally means
  # there's less information in the data and 
  # less evidence for relationship btw two variables 
    # than is assumed by OLS regression analysis







