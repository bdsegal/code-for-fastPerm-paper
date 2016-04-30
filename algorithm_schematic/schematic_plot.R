library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

pmf <- function(nx, ny, m){
  choose(nx, m) * choose(ny, m) / choose(nx + ny, min(nx, ny))
}

# expected value of W in partition m
muw <- function(m,lambda,N,muy,mux) {
  m/sqrt(lambda*(1-lambda)*N) * (muy-mux)
}

# g function of the expected value of w
g <- function(m, nx, ny, mux, muy){
  
  N <- nx+ny
  lambda <- nx/N
  
  numerator <- sqrt( (lambda*N)/(1-lambda) )*mux + muw(m,lambda,N,muy,mux)
  
  denominator <- sqrt( (1-lambda)*N/lambda ) * muy - muw(m,lambda,N,muy,mux)
  
  return( (1-lambda)/lambda * numerator/denominator )
}

  # check that it works -- at m=0, should be ratio of mux/muy -- ok
  # plot(g(1:nx,nx,ny,mux=2,muy=1))
  
  # plot(g(1:nx,nx,ny,mux=2,muy=1))

ER <- function(m,nx,ny,mux,muy){
  num <- mux + m/nx*(muy-mux)
  den <- muy + m/ny*(mux-muy)
  
  return(num/den)
}
  
gp <- function(m, nx, ny, mux, muy){
  
  N <- nx+ny
  lambda <- nx/N
  
  numerator <- sqrt( (lambda*N)/(1-lambda) ) * mux + sqrt( ((1-lambda)*N)/lambda ) * muy

  denominator <- (sqrt( ((1-lambda)*N)/lambda) * muy - muw(m,lambda,N,muy,mux))^2
  
  return((1-lambda)/lambda * numerator/denominator)
}

# gp(1:nx,nx,ny,mux=2,muy=1)

VarW <- function(m, nx, ny, sigmax2, sigmay2, mux, muy){
  
  N <- nx+ny
  lambda <- nx/N
  
  # variance terms
  EA_y <- 1/(lambda*(1-lambda)*N)*m/ny*(1-m/ny)*ny*(sigmay2 + muy^2)
  EA_x <- 1/(lambda*(1-lambda)*N)*m/nx*(1-m/nx)*nx*(sigmax2 + mux^2)

  # covariance terms
  EB_y <- 1/(lambda*(1-lambda)*N)*(m*(ny-m))/(ny^2*(ny-1))*ny*(ny-1)*muy^2

  EB_x <- 1/(lambda*(1-lambda)*N)*(m*(nx-m))/(nx^2*(nx-1))*nx*(nx-1)*mux^2

  return(EA_x + EA_y - EB_x - EB_y)
}

# VarW(0:nx,nx,ny,sigmax2,sigmay2,mux,muy)

VarR <- function(m,nx,ny,mux,muy,sigmax2,sigmay2){
  gp(m,nx,ny,mux,muy)^2 * VarW(m,nx,ny,sigmax2,sigmay2,mux,muy)
}
# plot(sqrt(VarR(0:nx,nx,ny,mux=1,muy=2,sigmax2,sigmay2)))

pApproxR <- function(m, nx, ny, sigmax2, sigmay2, mux, muy){
  
  N <- nx+ny
  lambda <- nx/N
  
  E_mean <- g(m, nx, ny, mux, muy)
  E_var <- VarR(m,nx,ny,mux,muy,sigmax2,sigmay2)
  
  t <- mux/muy

  eta_plus <- (t - E_mean)/sqrt(E_var)
  
  eta_minus_U <- (1/t - E_mean)/sqrt(E_var)
  
  eta_minus_L <- -E_mean/sqrt(E_var)

  # p=1-pnorm(eta_plus) + pnorm(eta_minus_U) - pnorm(eta_minus_L)
  # p=1-pt(eta_plus,df=df) + pt(eta_minus_U,df=df) - pt(eta_minus_L,df=df)
  
  p= pnorm(eta_plus,lower.tail=FALSE) #+ pnorm(eta_minus_U) #- pnorm(eta_minus_L)
  # p= pt(eta_plus, df=N-2, lower.tail=FALSE) #+ pnorm(eta_minus_U) #- pnorm(eta_minus_L)
  
  return(list(p=p, eta_plus=eta_plus))#, eta_minus_U=eta_minus_U,
        # eta_minus_L=eta_minus_L))

}

##########
# plots for Ratio statistic
nx=100
ny=100
mux <- sigmax2 <- 4
muy <- sigmay2 <- 2

# E_mean <- g(m=1:nx, nx, ny, mux, muy)
# E_var <- VarR(m=1:nx,nx,ny,mux,muy,sigmax2,sigmay2)

  
pm <- pApproxR(pmin(1:(nx-1), (nx-1):1), nx=nx, ny=ny, sigmax2=sigmax2, sigmay2=sigmay2, mux=mux, muy=muy)
p <- c(1,pm$p,1)

pDF <- data.frame(m=0:nx, p=p, logp=log(p,10))

M <- min(which(p<10^-3))-1
pDF$count <- with(pDF, p/min(p[1:M]))

fit <- glm(p~m, data=pDF[1:M,], family=poisson)

pDF$pred <- predict(fit, newdata=data.frame(m=pmin(pDF$m, rev(pDF$m))))

pDF$estimate <- c(rep(TRUE, M), rep(FALSE,nrow(pDF)-M))

dev.new(width=5, height=5)
ggplot(aes(x=m, y=log(p,10), color=estimate), data=pDF)+
  geom_point()+
  theme_bw(25)+
  # geom_point(aes(y=pred), color="red")+
  scale_color_manual(values=c("grey", "black"), guide=FALSE)+
  geom_hline(yintercept=-3)+
  labs(x="m", y=expression(paste("lo",g[10],"(p)", sep="")))+
  scale_y_continuous(breaks=c(0,-3,-10, -20, -30))

ggsave(file.path(paperPath,"algoExample.png"))