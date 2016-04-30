# Asymptotic trend in p-values across the partitions

library(ggplot2)
paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# Functions derived in the Appendix -----------------------

muw <- function(m,lambda,N,muy,mux) {
  # expected value of W in partition m
  m / sqrt(lambda * (1 - lambda) * N) * (muy-mux)
}

g <- function(m, nx, ny, mux, muy){
  # g function of the expected value of W
  
  N <- nx + ny
  lambda <- nx / N
  
  numerator <- sqrt((lambda * N)/(1 - lambda) ) * mux + 
                     muw(m, lambda, N, muy, mux)
  
  denominator <- sqrt((1 - lambda) * N / lambda) * muy - 
                      muw(m, lambda, N, muy, mux)
  
  return((1 - lambda) / lambda * numerator / denominator)
}

gp <- function(m, nx, ny, mux, muy){
  # derivative of the g function of the expected value of W
  
  N <- nx + ny
  lambda <- nx / N
  
  numerator <- sqrt((lambda * N) / (1 - lambda)) * mux + 
               sqrt(((1 - lambda) * N) / lambda) * muy

  denominator <- (sqrt(((1 - lambda) * N) / lambda) * muy - 
                  muw(m, lambda, N, muy, mux))^2
  
  return((1 - lambda) / lambda * numerator / denominator)
}

ER <- function(m,nx,ny,mux,muy){
  # Expected value of the ratio stat in partition m
  
  num <- mux + m / nx * (muy - mux)
  den <- muy + m / ny * (mux - muy)
  
  return(num / den)
}
  
VarW <- function(m, nx, ny, sigmax2, sigmay2, mux, muy){
  # Variance of W
  
  N <- nx + ny
  lambda <- nx / N
  
  # variance terms
  EAy <- 1 / (lambda * (1 - lambda) * N) * m / ny * (1 - m / ny) *
          ny * (sigmay2 + muy^2)
  EAx <- 1 / (lambda * (1 - lambda) * N) * m / nx * (1 - m / nx) *
          nx * (sigmax2 + mux^2)

  # covariance terms
  ECy <- 1 / (lambda * (1 - lambda) * N) * (m * (ny - m)) / 
          (ny^2 * (ny - 1)) * ny * (ny - 1) * muy^2

  ECx <- 1 / (lambda * (1 - lambda) * N) * (m * (nx - m)) /
          (nx^2 * (nx - 1)) * nx * (nx - 1) * mux^2

  return(EAx + EAy - ECx - ECy)
}

VarR <- function(m,nx,ny,mux,muy,sigmax2,sigmay2){
  # Variance of the ratio stat R
  
  gp(m, nx, ny, mux, muy)^2 * 
  VarW(m, nx, ny, sigmax2, sigmay2, mux, muy)
}

pApproxR <- function(m, nx, ny, sigmax2, sigmay2, mux, muy){
  # 
  
  N <- nx + ny
  lambda <- nx / N
  
  Emean <- g(m, nx, ny, mux, muy)
  Evar <- VarR(m, nx, ny, mux, muy, sigmax2, sigmay2)
  
  t <- mux / muy

  eta <- (t - Emean) / sqrt(Evar)
  
  p= pnorm(eta,lower.tail=FALSE)
  
  return(list(p=p, eta=eta))
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
ggsave(file.path(paperPath,"pApprox_ratio_100_100.png"))

data$p[which(data$p <= 10^-3)] <- 10^-3

dev.new(width=5, height=5)
ggplot(aes(x = m, y = log(p,10)), data = data)+
  geom_point(size = 1.5)+
  theme_bw(25)+
  labs(y = expression(paste("lo", g[10], "(p)")), x = "m")
ggsave(file.path(paperPath,"pApprox_100_100_10ToThe3Cutoff.png"))
