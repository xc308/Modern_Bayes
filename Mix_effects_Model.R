#==========================================
# Hoff Chp11 Linear & GL Mix effects Model
#==========================================

load("nelsSES.RData")
# 1993*4
head(nels)
#  sch_id sch_freelunch    stu_ses    stu_mathscore
#1   1011             6 -0.0599483      52.11
#2   1011             6  1.0516522      57.65
str(nels)
# 'data.frame':	1993 obs. of  4 variables:
#$ sch_id       : num  1011 1011 1011 1011 1011 ...
#$ sch_freelunch: num  6 6 6 6 6 6 6 6 6 6 ...
#$ stu_ses      : num  -0.0599 1.0517 -0.8635 -0.7966 -1.6135 ...
#$ stu_mathscore: num  52.1 57.6 66.4 44.7 40.6 ...

ids <- sort(unique(nels[, 1]))
m <- length(ids)  # how many schools? 100L

# To get the data for each school
Y <- list(); X <- list()
N <- NULL
for (j in 1:m) {
  Y[[j]] <- nels[nels$sch_id == ids[j], 4]  # sch-wise M score
  N[j] <- sum(nels$sch_id == ids[j])      # how many stus in each school
  xj <- nels[nels$sch_id == ids[j], 3]   # get the shool-wise ses
  xj <- xj - mean(xj)    # center for easy interpretation
  X[[j]] <- cbind(rep(1, N[j]), xj)
}

head(Y)
head(X)


# OLS fit for each school
BETA.LS <- NULL
SIG2.LS <- NULL
for (j in 1:m) {
  lmfit <- lm(Y[[j]] ~ -1 + X[[j]])
  
  BETA.LS <- rbind(BETA.LS, lmfit$coeff) # update BETA.LS
  SIG2.LS <- rbind(SIG2.LS, summary(lmfit)$sigma^2)
}



S2.LS<-BETA.LS<-NULL
for(j in 1:m) {
  fit<-lm(Y[[j]]~-1+X[[j]] )
  BETA.LS<-rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
} 


pdf("SES_VS_Mathscore.pdf", family = "Times")
par(mfrow = c(1, 3))
plot(range(nels[, 3]), range(nels[, 4]),
     type = "n", xlab = "SES", ylab = "Maths Score")

for (j in 1:m) {
  abline(a = BETA.LS[j, 1], b = BETA.LS[j, 2], 
        col = "gray")
}
abline(a = mean(BETA.LS[, 1]), b = mean(BETA.LS[, 2]), lwd = 2)


plot(N, BETA.LS[, 1], xlab = "Sample Size", ylab = "Intercept")
abline(h = mean(BETA.LS[, 1]), lwd = 2)

plot(N, BETA.LS[, 2], xlab = "Sample Size", ylab = "Slope")
abline(h = mean(BETA.LS[, 2]), lwd = 2)
dev.off()

## Remark:
  # 1. from the 1st plot Maths score against SES
    # a large majority of schools shows an increasing 
      # relationship betw the M score with the rise of SES

    # but few shows a negative relationship betw
      # the M score with the rise of SES

    
    # 2. from the 2nd and 3rd plots: 
      # regression coeffs vs sample size
      # the schools with higher sample size, 
        # the corresponding coeff are closer to the mean
     # extreme coeffs are with schools with small sample size

      # the smaller the sample size, the more probable
      # the under-representative data are sampled 
      # and extreme OLS coeff estimates are produced

    # remedy to the small sample size problem
      # is by sharing information across different schools (groups)
      # to stablize the estimates of small sample size 
      # via hierachical model




    



