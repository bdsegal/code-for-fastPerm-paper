library(fastPerm)

# hours of pain relief (Lehman, 1975, p. 37) ----------------------------------
x <- c(6.8, 3.1, 5.8, 4.5, 3.3, 4.7, 4.2, 4.9)
y <- c(4.4, 2.5, 2.8, 2.1, 6.6, 0, 4.8, 2.3)
nx <- length(x)
ny <- length(y)
N <- nx + ny

nRep <- 100
pPred <- rep(NA, nRep)
for (i in 1:nRep){
  pPred[i] <- fastPerm(x, y, testStat = diffMean)$pPred
}

B <- 1e5
z <- c(x, y)
t0 <- abs(mean(x) - mean(y))
tStar <- rep(NA, B)
for (b in 1:B) {
  star <- sample(1:N)
  xStar <- z[star[1:nx]]
  yStar <- z[star[(nx+1):N]]
  tStar[b] <- abs(mean(xStar) - mean(yStar))
}

summary(pPred)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.06607 0.07936 0.08289 0.08305 0.08761 0.09522 

fastPermAsym(x, y, testStat = diffMean)$pNorm
# [1] 0.09227375

mean(c(tStar, t0) >= t0)
# [1] 0.100629


# Effect of analgesia, (Lehman, 1975, p. 92)
x <- c(17.9, 13.3, 10.6, 7.6, 5.7, 5.6, 5.4, 3.3, 3.1, 0.9)
y <- c(7.7, 5, 1.7, 0, -3, -3.1, -10.5)
nx <- length(x)
ny <- length(y)
N <- nx + ny

nRep <- 100
pPred <- rep(NA, nRep)
for (i in 1:nRep){
  pPred[i] <- fastPerm(x, y, testStat = diffMean)$pPred
}

B <- 1e5
z <- c(x, y)
t0 <- abs(mean(x) - mean(y))
tStar <- rep(NA, B)
for (b in 1:B) {
  star <- sample(1:N)
  xStar <- z[star[1:nx]]
  yStar <- z[star[(nx+1):N]]
  tStar[b] <- abs(mean(xStar) - mean(yStar))
}

summary(pPred)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.008789 0.010610 0.011500 0.011570 0.012400 0.015630 

fastPermAsym(x, y, testStat = diffMean)$pNorm
# [1] 0.01282504

mean(c(tStar, t0) >= t0)
# [1] 0.01148989
