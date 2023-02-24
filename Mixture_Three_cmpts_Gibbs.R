#==========================
# Three-components Mixture
#==========================

Mix_Data <- read.csv("Lab8Mixture.csv")
head(Mix_Data)

mean(Mix_Data)
Mix_Data[500, ] <- 2.45394289200374

names(Mix_Data) <- "X"
head(Mix_Data)
str(Mix_Data)
# 'data.frame':	500 obs. of  1 variable:

Dat_Mix <- Mix_Data[, 1]
str(Dat_Mix)
# num [1:500] 3.41 1.97 2.46 4.88 -1.06 ...
mean(Dat_Mix)
# [1] 2.781925

var(Dat_Mix)
# [1] 4.637941

cdf = 0.5
x <- runif(1)
which(x <= cdf)
[1]












