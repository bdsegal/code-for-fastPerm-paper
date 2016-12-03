library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"


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

  eta_plus <- (t - E_mean)/sqrt(E_var)
   
  p <- pnorm(eta_plus,lower.tail=FALSE)
  
  return(list(p=p, eta_plus=eta_plus))
}

pApproxRlin <- function(m, nx, ny, sigmax2, sigmay2, mux, muy){
  
  E_mean <- gLin(m, nx, ny, mux, muy)
  E_var <- VarRLin(m,nx,ny,mux,muy,sigmax2,sigmay2)
  
  t <- abs(mux - muy)

  eta_plus <- (t - E_mean)/sqrt(E_var)
  p <- pnorm(eta_plus,lower.tail=FALSE)
  
  return(list(p=p, eta_plus=eta_plus))
}


##########
# plots for Ratio statistic
nx=100
ny=100
mux <- sigmax2 <- 4
muy <- sigmay2 <- 2

pm1 <- pApproxR(pmin(1:(nx-1), (nx-1):1), nx = nx, ny = ny,
               sigmax2 = sigmax2, sigmay2 = sigmay2,
               mux = mux, muy = muy)

pm2 <- pApproxR(pmin(1:(nx-1), (nx-1):1), nx = ny, ny = nx,
               sigmax2 = sigmay2, sigmay2 = sigmax2,
               mux = muy, muy = mux)

p <- c(1, pm1$p + pm2$p, 1)

pDF <- data.frame(m=0:nx, p=p, logp=log(p,10))

M <- min(which(p<10^-3))-1
# M <- nx / 2
pDF$count <- with(pDF, p/min(p[1:M]))


pDF$estimate <- c(rep(TRUE, M), rep(FALSE,nrow(pDF)-M))

dev.new(width=5, height=5)
ggplot(aes(x=m, y=log(p,10), color=estimate), data=pDF)+
  geom_point()+
  theme_bw(25)+
  scale_color_manual(values=c("grey", "black"), guide=FALSE)+
  geom_hline(yintercept=-3)+
  labs(x="m", y=expression(paste("lo",g[10],"(p)", sep="")))+
  scale_y_continuous(breaks=c(0,-3,-10, -20, -30))

ggsave(file.path(paperPath,"algoExample.png"))