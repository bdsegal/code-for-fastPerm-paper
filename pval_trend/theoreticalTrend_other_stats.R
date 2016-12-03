# Asymptotic trend in p-values across the partitions

library(ggplot2)
paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# Functions derived in the Appendix -----------------------

# expected value of W in partition m
muw <- function(m, ybar, xbar) {

  m * (ybar - xbar)
}

# # g function of the expected value of w
# g <- function(m, nx, ny, xbar, ybar){

#   mu <- muw(m = m, xbar = xbar, ybar = ybar)
    
#   numerator <- exp(xbar + (1 / nx) * mu)
#   denominator <- ybar - (1 / ny) * mu
#   return(numerator / denominator)
# }

# plot(g(1:nx,nx,ny,xbar=2,ybar=1))

# gp <- function(m, nx, ny, xbar, ybar){

#   mu <- muw(m = m, xbar = xbar, ybar = ybar)
  
#   num1 <- (1 / nx) * (xbar + (1 / nx) * mu) * exp(xbar + (1 / nx) * mu) * (ybar - (1 / ny) * mu)
#   num2 <- (1 / ny) * exp(xbar + (1 / nx) * mu)
#   den <- (ybar - (1 / ny) * mu)^2
#   return((num1 + num2) / den)
# }

# # plot(gp(1:nx,nx,ny,xbar=2,ybar=1))

# g function of the expected value of w
# g <- function(m, nx, ny, xbar, ybar){

#   mu <- muw(m = m, xbar = xbar, ybar = ybar)
    
#   numerator <- xbar + (1 / nx) * mu
#   denominator <- ybar - (1 / ny) * mu
#   return(log(numerator / denominator))
# }

# # plot(g(1:nx,nx,ny,xbar=2,ybar=1))

# gp <- function(m, nx, ny, xbar, ybar){

#   mu <- muw(m = m, xbar = xbar, ybar = ybar)
  
#   val <- (1 / nx) / (xbar + (1 / nx) * mu) + (1 / ny) / (ybar - (1 / ny) * mu)
#   return(val)
# }

# plot(gp(1:nx,nx,ny,xbar=2,ybar=1))

g <- function(m, nx, ny, xbar, ybar){

  mu <- muw(m = m, xbar = xbar, ybar = ybar)
    
  num <- xbar + (1 / nx) * mu
  den <- ybar - (1 / ny) * mu
  lin <- xbar - ybar + (1 / nx + 1 / ny) * mu
  return((num / den) + lin)
}

# plot(g(1:nx,nx,ny,xbar=2,ybar=1))

gp <- function(m, nx, ny, xbar, ybar){

  mu <- muw(m = m, xbar = xbar, ybar = ybar)
  
  den <- (1 / nx) * (ybar - (1 / ny) * mu) + (1 / ny) * (xbar + (1 / nx) * mu)
  num <- (ybar - (1 / ny) * mu)^2
  lin <- 1 / nx + 1 / ny
  return(den / num + lin)
}

VarW <- function(m, nx, ny, sigmax2, sigmay2){
  
  return( m * (ny - m) / ny * sigmay2 + m * (nx - m) / nx * sigmax2)
}

# plot(VarW(0:nx,nx,ny,sigmax2,sigmay2))

VarR <- function(m,nx,ny,xbar,ybar,sigmax2,sigmay2){

  gp(m,nx,ny,xbar,ybar)^2 * VarW(m,nx,ny,sigmax2,sigmay2)
}

# plot(VarR(0:nx,nx,ny,xbar,ybar,sigmax2,sigmay2))
# plot(sqrt(VarR(0:nx,nx,ny,xbar,ybar,sigmax2,sigmay2)))

pApproxR <- function(m, nx, ny, sigmax2, sigmay2, xbar, ybar){
  
  E_mean <- g(m, nx, ny, xbar, ybar)
  E_var <- VarR(m,nx,ny,xbar,ybar,sigmax2,sigmay2)
  
  # t <- max(xbar/ybar, ybar/xbar)
  # t <- max(exp(xbar) / ybar, exp(ybar) / ybar)
  t <- max(xbar / ybar, ybar / xbar) + abs(xbar - ybar)

  eta_plus <- (t - E_mean)/sqrt(E_var)
   
  p <- pnorm(eta_plus,lower.tail=FALSE)
  
  return(list(p=p, eta_plus=eta_plus))
}

# Plots ---------------------------------------------------

nx=100
ny=100

xbar <- sigmax2 <- 3
ybar <- sigmay2 <- 2
# sigmax2 <- sigmay2 <- 1

m <- pmin(1:(nx-1), (nx-1):1)
pm1 <- pApproxR(pmin(1:(nx-1), (nx-1):1), nx = nx, ny = ny,
               sigmax2 = sigmax2, sigmay2 = sigmay2,
               xbar = xbar, ybar = ybar)

pm2 <- pApproxR(pmin(1:(nx-1), (nx-1):1), nx = ny, ny = nx,
               sigmax2 = sigmay2, sigmay2 = sigmax2,
               xbar = ybar, ybar = xbar)

# plot(log(pm1$p))
# plot(log(pm2$p))

p <- c(1, pm1$p + pm2$p, 1)

data <- data.frame(m=0:nx, p=p)

dev.new(width=5, height=5)
ggplot(aes(x = m, y = log(p,10)), data = data)+
  geom_point(size = 1.5)+
  theme_bw(25)+
  labs(y = expression(paste("lo", g[10], "(p)")), x = "m")
# ggsave(file.path(paperPath,"pApprox_ratio_100_100_other1.png"))



# studentized stat... run simulations...
nx <- 100
ny <- 150
xbar <- 0.75
ybar <- 0
N <- nx + ny
x <- rnorm(n = nx, mean = xbar, sd = 3)
y <- rnorm(n = ny, mean = ybar, sd = 1)

z <- c(x, y)
B <- 1e4
Tstud <- rep(NA, B)
T <- rep(NA, B)

t0 <- mean(x) - mean(y)
tStud0 <- t0 / sqrt(var(x) / nx + var(y) / ny)
for (b in 1:B) {
  pi <- sample(1:N, replace = FALSE)
  xStar <- z[pi[1:nx]]
  yStar <- z[pi[(nx+1):N]]
  T[b] <- mean(xStar) - mean(yStar)
  Tstud[b] <- (mean(xStar) - mean(yStar)) / sqrt(var(xStar) / nx + var(yStar) / ny)
}

library(fastPerm)
fpStud <- fastPerm(x, y, testStat = diffMeanStudent)
fp <- fastPerm(x, y, testStat = diffMean)

data.frame( method = c("t.test", "perm student", "perm", "fp student", "fp"),
  p = c(t.test(x, y)$p.value,
  mean(abs(c(Tstud, tStud0)) >= abs(tStud0)),
  mean(abs(c(T, t0)) >= abs(t0)),
  fpStud$pPred,
  fp$pPred))

fpStud
fp

hist(Tstud)
abline(v = tStud0, col = "red")
dev.new()
hist(T)
abline(v = t0, col = "red")