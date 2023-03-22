#==================================
# Hoff Chp11 GL Mix effects Model
#==================================
# Analysis of tumor location data

load("tumorLocation.RData")
head(tumorLocation)

m <- nrow(tumorLocation)
Loc <- seq(1, 20, by = 1)
xs <- Loc / 20

Y <- tumorLocation

plot(c(0, 1), range(Y), xlab = "location", ylab = "number of tumors",
     type = "n")
for (j in 1:m) {
  lines(xs, Y[j, ], col = "gray", type = "l")
  # joint the correponding points
}


# average of 21 mice at each location
AVG_Loci <- apply(Y, 2, mean)


#------------------------
# To fit the log average
#------------------------

# 1st, construct design matrix (DM)


X <- cbind(rep(1, ncol(Y)), poly(xs, 4, raw = T))
# raw = T, each of polynomials are not orthogonal to constant 1 col

fit2 <- lm(AVG_Loci ~ -1 + X[, 1:3]) # quadratic
fit3 <- lm(AVG_Loci ~ -1 + X[, 1:4]) # cubic
fit4 <- lm(AVG_Loci ~ -1 + X[, 1:5]) # quatic


y2 <- X[, 1:3] %*% fit2$coefficients
y3 <- X[, 1:4] %*% fit3$coefficients
y4 <- X[, 1:5] %*% fit4$coefficients

# X[, 1:3]: 20*3, 
# fit2$coefficients : 3*1

# X[, 1:4]: 20*4
# fit3$coefficients: 4*1

plot(xs, AVG_Loci,
     ylim = range(c(AVG_Loci, y2, y3, y4)),
     type = "n", xlab = "locations", 
     ylab = "log avg number of tumors",
     lwd = 3)

lines(xs, AVG_Loci, lwd = 3)
lines(xs, y2)
points(xs, y2, pch = "2", col = "black")

lines(xs, y3, col = "gray")
points(xs, y3, pch = "3", col = "red")

lines(xs, y4, col = "blue")
points(xs, y4, pch = "4", col = "blue")

## 
# poly with degrees of 4 is the best fit to the 
  # log avg of number of tumors

#






