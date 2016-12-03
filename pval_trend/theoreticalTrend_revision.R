# Asymptotic trend in p-values across the partitions

library(ggplot2)
paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# Functions derived in the Appendix -----------------------

# expected value of W in partition m
muw <- function(m, ybar, xbar) {

  m * (ybar - xbar)
}

# g function of the expected value of w
g <- function(m, nx, ny, xbar, ybar){
    
  numerator <- nx * xbar + muw(m = m, xbar = xbar, ybar = ybar)
  denominator <- ny * ybar - muw(m = m, xbar = xbar, ybar = ybar)
  return(ny / nx * numerator / denominator)
}

# plot(g(1:nx,nx,ny,xbar=2,ybar=1))

gp <- function(m, nx, ny, xbar, ybar){
  
  numerator <- ny * ybar + nx * xbar
  denominator <- (ny * ybar - muw(m = m, xbar = xbar, ybar = ybar))^2
  return(ny / nx * numerator / denominator)
}

# g function of the expected value of w for linear T = |x - y|
gLin <- function(m, nx, ny, xbar, ybar){
    
  return(xbar - ybar + (1 / nx + 1 / ny) * muw(m = m, xbar = xbar, ybar = ybar))
}

gpLin <- function(nx, ny){
  
  return((1 / nx + 1 / ny))
}

# gp(1:nx,nx,ny,xbar=2,ybar=1)

VarW <- function(m, nx, ny, sigmax2, sigmay2){
  
  return( m * (ny - m) / ny * sigmay2 + m * (nx - m) / nx * sigmax2)
}

# VarW(0:nx,nx,ny,sigmax2,sigmay2,xbar,ybar)

VarR <- function(m,nx,ny,xbar,ybar,sigmax2,sigmay2){

  gp(m,nx,ny,xbar,ybar)^2 * VarW(m,nx,ny,sigmax2,sigmay2)
}

plot(VarR(0:nx,nx,ny,xbar,ybar,sigmax2,sigmay2))

VarRLin <- function(m,nx,ny,xbar,ybar,sigmax2,sigmay2){

  return(gpLin(nx,ny)^2 * VarW(m,nx,ny,sigmax2,sigmay2))
}

# plot(sqrt(VarR(0:nx,nx,ny,xbar=1,ybar=2,sigmax2,sigmay2)))

pApproxR <- function(m, nx, ny, sigmax2, sigmay2, xbar, ybar){
  
  E_mean <- g(m, nx, ny, xbar, ybar)
  E_var <- VarR(m,nx,ny,xbar,ybar,sigmax2,sigmay2)
  
  t <- max(xbar/ybar, ybar/xbar)

  eta_plus <- (t - E_mean)/sqrt(E_var)
   
  p <- pnorm(eta_plus,lower.tail=FALSE)
  
  return(list(p=p, eta_plus=eta_plus))
}

pApproxRlin <- function(m, nx, ny, sigmax2, sigmay2, xbar, ybar){
  
  E_mean <- gLin(m, nx, ny, xbar, ybar)
  E_var <- VarRLin(m,nx,ny,xbar,ybar,sigmax2,sigmay2)
  
  t <- abs(xbar - ybar)

  eta_plus <- (t - E_mean)/sqrt(E_var)
  p <- pnorm(eta_plus,lower.tail=FALSE)
  
  return(list(p=p, eta_plus=eta_plus))
}

# Plots ---------------------------------------------------

nx=100
ny=100

xbar <- sigmax2 <- 4
ybar <- sigmay2 <- 2
       
pm1 <- pApproxR(pmin(1:(nx-1), (nx-1):1), nx = nx, ny = ny,
               sigmax2 = sigmax2, sigmay2 = sigmay2,
               xbar = xbar, ybar = ybar)

pm2 <- pApproxR(pmin(1:(nx-1), (nx-1):1), nx = ny, ny = nx,
               sigmax2 = sigmay2, sigmay2 = sigmax2,
               xbar = ybar, ybar = xbar)

p <- c(1, pm1$p + pm2$p, 1)

data <- data.frame(m=0:nx, p=p)

dev.new(width=5, height=5)
ggplot(aes(x = m, y = log(p,10)), data = data)+
  geom_point(size = 1.5)+
  theme_bw(25)+
  labs(y = expression(paste("lo", g[10], "(p)")), x = "Partition m")

ggsave(file.path(paperPath,"pApprox_ratio_100_100_revision.png"))

data$p[which(data$p <= 10^-3)] <- 10^-3

dev.new(width=5, height=5)
ggplot(aes(x = m, y = log(p,10)), data = data)+
  geom_point(size = 1.5)+
  theme_bw(25)+
  labs(y = expression(paste("lo", g[10], "(p)")), x = "Partition m")
ggsave(file.path(paperPath,"pApprox_100_100_10ToThe3Cutoff_revision.png"))


# difference
nx=100
ny=100

xbar <- 4
ybar <- 2
sigmax2 <- sigmay2 <- 1
# xbar <- sigmax2 <- 4
# ybar <- sigmay2 <- 2
       
pm1 <- pApproxRlin(pmin(1:(nx-1), (nx-1):1), nx = nx, ny = ny,
               sigmax2 = sigmax2, sigmay2 = sigmay2,
               xbar = xbar, ybar = ybar)

pm2 <- pApproxRlin(pmin(1:(nx-1), (nx-1):1), nx = ny, ny = nx,
               sigmax2 = sigmay2, sigmay2 = sigmax2,
               xbar = ybar, ybar = xbar)

p <- c(1, pm1$p + pm2$p, 1)
# p <- c(1, pm1$p, 1)

# all <- as.data.frame(rbind(cbind(m=1:(nx-1), p = pm1$p, type = 1), 
#              cbind(m=1:(nx-1), p = pm2$p, type = 2)))
# ggplot(aes(x = m, y = log(p), color = type), data = all)+geom_point()

data <- data.frame(m=0:nx, p=p)

dev.new(width=5, height=5)
ggplot(aes(x = m, y = log(p,10)), data = data)+
  geom_point(size = 1.5)+
  theme_bw(22)+
  labs(y = expression(paste("lo", g[10], "(p)")), x = "Partition m")
ggsave(file.path(paperPath,"pApprox_Diff_100_100.png"))
