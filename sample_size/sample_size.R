library(ggplot2)
library(reshape2)
library(rgl)

# paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

PiLgamma <- function(nx, ny){
  
  nMin <- min(nx, ny)
  out <- rep(NA, nMin + 1)
  names(out) <- 0:min(nx, ny)
  
  for (m in 0:nMin){
    out[m + 1] <- exp(
    	lgamma(nx + 1) + lgamma(ny + 1) + lgamma(nx + ny - nMin + 1) +
    	lgamma(nMin + 1) - lgamma(nx - m + 1) - lgamma(ny - m + 1) -
    	lgamma(nx + ny +1) - 2*lgamma(m+1)
    	)
  }
  
  return(out)
}

# expected value of W in partition m
muw <- function(m, muy, mux) {

  m * (muy - mux)
}

# g function of the expected value of w
g <- function(m, nx, ny, mux, muy){
    
  numerator <- nx * mux + muw(m = m, mux = mux, muy = muy)
  denominator <- ny * muy - muw(m = m, mux = mux, muy = muy)
  return(ny / nx * numerator / denominator)
}


gp <- function(m, nx, ny, mux, muy){
  
  numerator <- ny * muy + nx * mux
  denominator <- (ny * muy - muw(m = m, mux = mux, muy = muy))^2
  return(ny / nx * numerator / denominator)
}

# g function of the expected value of w for linear T = |x - y|
gLin <- function(m, nx, ny, mux, muy){
    
  return(mux - muy + (1 / nx + 1 / ny) * muw(m = m, mux = mux, muy = muy))
}

gpLin <- function(nx, ny){
  
  return((1 / nx + 1 / ny))
}

# gp(1:nx,nx,ny,mux=2,muy=1)

VarW <- function(m, nx, ny, sigmax2, sigmay2){
  
  return( m * (ny - m) / ny * sigmay2 + m * (nx - m) / nx * sigmax2)
}

# VarW(0:nx,nx,ny,sigmax2,sigmay2,mux,muy)

VarR <- function(m,nx,ny,mux,muy,sigmax2,sigmay2){

  gp(m,nx,ny,mux,muy)^2 * VarW(m,nx,ny,sigmax2,sigmay2)
}

VarRLin <- function(m,nx,ny,mux,muy,sigmax2,sigmay2){

  return(gpLin(nx,ny)^2 * VarW(m,nx,ny,sigmax2,sigmay2))
}

# plot(sqrt(VarR(0:nx,nx,ny,mux=1,muy=2,sigmax2,sigmay2)))

pApproxR <- function(m, nx, ny, sigmax2, sigmay2, mux, muy){
  
  E_mean <- g(m, nx, ny, mux, muy)
  E_var <- VarR(m,nx,ny,mux,muy,sigmax2,sigmay2)
  
  t <- max(mux/muy, muy/mux)

  xi <- (t - E_mean)/sqrt(E_var)
   
  p <- pnorm(xi,lower.tail=FALSE)
  
  return(list(p=p, xi=xi))
}

pApproxRlin <- function(m, nx, ny, sigmax2, sigmay2, mux, muy){
  
  E_mean <- gLin(m, nx, ny, mux, muy)
  E_var <- VarRLin(m,nx,ny,mux,muy,sigmax2,sigmay2)
  
  t <- abs(mux - muy)

  xi <- (t - E_mean)/sqrt(E_var)
  p <- pnorm(xi,lower.tail=FALSE)
  
  return(list(p=p, xi=xi))
}

mStopRatioMean <- function(nx, ny, sigmax2, sigmay2, mux, muy, B=1000, 
                           pFun = pApproxR, plot = FALSE){

  # swap x and y if xbar < ybar
  if (mux < muy) {
    muytemp <- muy
    muy <- mux
    mux <- muytemp
    
    sigmay2temp <- sigmay2
    sigmay2 <- sigmax2
    sigmax2 <- sigmay2temp
    
    nytemp <- ny
    ny <- nx
    nx <- nytemp
  }
  
  nMax <- min(nx, ny)
  pmf <- PiLgamma(nx, ny)
  mMax <- min(which(pmf == max(pmf)))

  # get p-value
  if (length(mMax)==1 & mMax[1]==nMax){
    m <- 1:mMax
  } else if (length(mMax)==1 & mMax[1]<nMax){
    m <- c(1:mMax,((mMax-1):0)[1:(nMax-mMax)])
  } else if (length(mMax)==2 & mMax[2]==nMax){
    m <- c(1:mMax[1], mMax[1])
  } else {
    m <- c(1:mMax[1], (mMax[1]:0)[1:(nMax - mMax[1])])
  }

  a1 <- pFun(m, nx, ny, sigmax2, sigmay2, mux, muy)
  a2 <- pFun(m, ny, nx, sigmay2, sigmax2, muy, mux)
  xi1 <- a1$xi
  xi2 <- a2$xi
  if (nx != ny) { 
    pvec <- c(1, pnorm(xi1, lower.tail=FALSE) + pnorm(xi2, lower.tail = FALSE))
  } else { 
    pvec <- c(1, pnorm(xi1[-length(xi1)], lower.tail=FALSE) + 
                  pnorm(xi2[-length(xi2)], lower.tail=FALSE), 1)
  }  
  pVal <- pmf %*% pvec
  
  # get mStop
  m <- 1:min(mMax,(nMax-1))
  u <- qnorm(1-1/B)
  xi <- pFun(m, nx, ny, sigmax2, sigmay2, mux, muy)$xi
  nu <- which(xi > u)
  if (length(nu) > 0) {
    mStop <- m[min(nu)]
  } else {
    mStop <- mMax
  }

  if (plot) {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(x = m, y = xi, pch=19, ylab = expression(paste(xi, "(m)")), 
         main= bquote(paste(m["stop"]^"asym", " = ", .(mStop), sep = "")),
         cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4)
    abline(a = u, b = 0, col="red")
  }
  
  out <- c(mStop, mMax, pVal)
  names(out) <- c("mStop", "mMax", "p")

  return(out)
}

##########
# plots for Ratio statistic
nx=100
ny=100
mux <- sigmax2 <- 3
muy <- sigmay2 <- 2

results <- data.frame(n = 3:100)
results$mStop <- NA
results$mMax <- NA
results$p <- NA
for (i in 1:nrow(results)) {
  n <- results$n[i]
  results[i, -1] <- mStopRatioMean(nx = n, ny = n, sigmax2, sigmay2, mux, muy, B=1000, 
                                   pFun = pApproxR,  plot = FALSE)
}

nHat <- results$n[min(which(results$mStop >= 4))]
with(results, plot(x = n, y = mStop))
abline(v = nHat, lty = 2)
abline(h = 4, lty = 2)

with(results, plot3d(x = n, y = log(p, 10), z = mStop))
qplot(x = n, y = )

cutoff <- c(4, 5, 10)
ratio <- list()
nIterMax <- 1000

for (k in 1:length(cutoff)) {
  print(k)

  ratTab <- data.frame(muSigmaY = 2, muSigmaX = c(3, 4, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.6, 6.7, 6.8, 6.9, 7))
  ratTab$nHat <- NA
  ratTab$mStop <- NA
  ratTab$mMax <- NA
  ratTab$p <- NA

  for (i in 1:nrow(ratTab)) {
    # print(i)
    flag <- FALSE
    n <- 3
    while (n < nIterMax & !flag) {
      out <- mStopRatioMean(nx = n, ny = n, 
                            sigmax2 = ratTab$muSigmaX[i], sigmay2 = ratTab$muSigmaY[i], 
                            mux = ratTab$muSigmaX[i], muy = ratTab$muSigmaY[i], 
                            pFun = pApproxR, B=1000, plot = FALSE)
                            
      if(out["mStop"] >= cutoff[k]) {
        ratTab[i, -(1:2)] <- c(n, out)
        flag <- TRUE
      }
      n <- n + 1
    }
  }
  ratio[[k]] <- ratTab
  # names(ratio[[k]]) <- cutoff[k]
}
names(ratio) <- paste("c = ", cutoff, sep = "")

ratioPrint <- cbind(cbind(ratio[[1]][, c(1, 2, 3, 6)], ratio[[2]][, c(3, 6)]), ratio[[3]][, c(3, 6)])
signif(ratioPrint, 3)
#    muSigmaY muSigmaX nHat           p nHat           p nHat           p
# 1         2     3.00    5 2.40235e-01    8 1.32676e-01   18 2.44091e-02
# 2         2     4.00    6 2.36705e-02   10 2.68493e-03   58 3.72757e-12
# 3         2     5.00   13 2.41629e-05   25 6.55798e-09   NA          NA
# 4         2     5.25   16 1.28342e-06   33 5.59347e-12   NA          NA
# 5         2     5.50   19 6.03878e-08   46 5.52717e-17   NA          NA
# 6         2     5.75   24 4.23241e-10   67 3.23284e-25   NA          NA
# 7         2     6.00   31 4.06009e-13  111 1.63894e-42   NA          NA
# 8         2     6.25   40 4.33295e-17  240 1.42459e-93   NA          NA
# 9         2     6.50   55 1.14308e-23   NA          NA   NA          NA
# 10        2     6.60   63 3.25951e-27   NA          NA   NA          NA
# 11        2     6.70   74 4.54328e-32   NA          NA   NA          NA
# 12        2     6.80   87 7.73300e-38   NA          NA   NA          NA
# 13        2     6.90  105 7.80347e-46   NA          NA   NA          NA
# 14        2     7.00  130 6.02572e-57   NA          NA   NA          NA


##########
# plots for Difference statistic
# plots for Ratio statistic
mux <- 2.5
muy <- 0

sigmax2 <- sigmay2 <- 1

results <- data.frame(n = 3:300)
results$mStop <- NA
results$mMax <- NA
results$p <- NA
for (i in 1:nrow(results)) {
  n <- results$n[i]
  results[i, -1] <- mStopRatioMean(nx = n, ny = n, sigmax2, sigmay2, mux, muy, B=1000, 
                                   pFun = pApproxRlin,  plot = FALSE)
}

nHat <- results$n[min(which(results$mStop >= 4))]
with(results, plot(x = n, y = mStop))
abline(v = nHat, lty = 2)
abline(h = 4, lty = 2)

with(results, plot3d(x = n, y = log(p, 10), z = mStop))
qplot(x = n, y = )

nIterMax <- 1000

diff <- list()
for (k in 1:length(cutoff)) {
  print(k)
  
  diffTab <- data.frame(muy = 0, mux = c(1.5, 2, 2.2, 2.25, 2.3, 2.4, 2.45, 2.475, 2.48, 2.49, 2.5))
  diffTab$nHat <- NA
  diffTab$mStop <- NA
  diffTab$mMax <- NA
  diffTab$p <- NA
  
  for (i in 1:nrow(diffTab)) {
    # print(i)
    flag <- FALSE
    n <- 3
    while (n < nIterMax & !flag) {
      out <- mStopRatioMean(nx = n, ny = n, 
                            sigmax2 = 1, sigmay2 = 1, 
                            mux = diffTab$mux[i], muy = diffTab$muy[i], 
                            pFun = pApproxRlin, B=1000, plot = FALSE)
                            
      if(out["mStop"] >= cutoff[k]) {
        diffTab[i, -(1:2)] <- c(n, out)
        flag <- TRUE
      }
      n <- n + 1
    }
  }
  diff[[k]] <- diffTab
}
names(diff) <- paste("c = ", cutoff, sep = "")

diffPrint <- cbind(cbind(diff[[1]][, c(1, 2, 3, 6)], diff[[2]][, c(3, 6)]), diff[[3]][, c(3, 6)])

diffPrint
   muy   mux nHat            p nHat            p nHat  p
# 1    0 1.500    5 5.373160e-02    8 9.505679e-03   NA NA
# 2    0 2.000    9 7.746931e-04   25 3.228155e-08   NA NA
# 3    0 2.200   13 2.136954e-05   NA           NA   NA NA
# 4    0 2.250   15 3.700476e-06   NA           NA   NA NA
# 5    0 2.300   18 3.119071e-07   NA           NA   NA NA
# 6    0 2.400   32 4.036866e-12   NA           NA   NA NA
# 7    0 2.450   53 2.288307e-19   NA           NA   NA NA
# 8    0 2.475   80 1.283136e-28   NA           NA   NA NA
# 9    0 2.480   89 1.065835e-31   NA           NA   NA NA
# 10   0 2.490  115 1.508698e-40   NA           NA   NA NA
# 11   0 2.500  165 1.373760e-57   NA           NA   NA NA

