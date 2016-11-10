# Asymptotic trend in p-values across the partitions

library(ggplot2)
paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# Functions derived in the Appendix -----------------------

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
  
  t <- mux/muy

  eta_plus <- (t - E_mean)/sqrt(E_var)
   
  p <- pnorm(eta_plus,lower.tail=FALSE)
  
  return(list(p=p, eta_plus=eta_plus))
}

pApproxRlin <- function(m, nx, ny, sigmax2, sigmay2, mux, muy){
  
  E_mean <- gLin(m, nx, ny, mux, muy)
  E_var <- VarRLin(m,nx,ny,mux,muy,sigmax2,sigmay2)
  
  t <- mux - muy

  eta_plus <- (t - E_mean)/sqrt(E_var)
  p <- pnorm(eta_plus,lower.tail=FALSE)
  
  return(list(p=p, eta_plus=eta_plus))
}

# Plots ---------------------------------------------------

nx=100
ny=100

mux <- sigmax2 <- 4
muy <- sigmay2 <- 2
       
pm <- pApproxR(pmin(1:(nx-1), (nx-1):1), nx = nx, ny = ny,
               sigmax2 = sigmax2, sigmay2 = sigmay2,
               mux = mux, muy = muy)

p <- c(1, pm$p, 1)

data <- data.frame(m=0:nx, p=p)

dev.new(width=5, height=5)
ggplot(aes(x = m, y = log(p,10)), data = data)+
  geom_point(size = 1.5)+
  theme_bw(25)+
  labs(y = expression(paste("lo", g[10], "(p)")), x = "m")
ggsave(file.path(paperPath,"pApprox_ratio_100_100_revision.png"))

data$p[which(data$p <= 10^-3)] <- 10^-3

dev.new(width=5, height=5)
ggplot(aes(x = m, y = log(p,10)), data = data)+
  geom_point(size = 1.5)+
  theme_bw(25)+
  labs(y = expression(paste("lo", g[10], "(p)")), x = "m")
ggsave(file.path(paperPath,"pApprox_100_100_10ToThe3Cutoff_revision.png"))
