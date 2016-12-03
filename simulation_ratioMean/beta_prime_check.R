# Visual check to make sure beta prime distribution is correct

library(ggplot2)
library(PearsonDS)

Beta <- function(a, b) {
  exp(lgamma(a) + lgamma(b) - lgamma(a + b))
}
# beta prime density (also called the inverted beta)
dbetaPrime <- function(x,a,b){
  (x^(a-1)*(1+x)^(-a-b))/Beta(a,b)
}

# group sample sizes and common rate parameter
nx <- 10
ny <- 100
lambda <- 1

# number of simulated ratio statistics; r holds the ratio of sums
Niter <- 1e5
r1 <- rep(NA, Niter)
r2 <- rep(NA, Niter)

for (i in 1:Niter){

  x <- rexp(n=nx, rate=lambda)
  y <- rexp(n=ny, rate=lambda)
  
  r1[i] <- sum(x)/sum(y)
  r2[i] <- sum(y)/sum(x)
}

# empirical cdf (blakck) with theoretical overlaid (red)
t0 <- 1.5
t1 <- (nx / ny) * t0
t2 <- (ny / nx) * t0
z <- seq(0,2,.01)

hist(r1)
dev.new()
plot(ecdf(r1), main = "lambda = 1")
lines(z,pbeta(z/(1+z), nx, ny, lower.tail=TRUE), col="red")
lines(z, ppearsonVI(q = z, a = nx, b = ny, location = 0, scale = 1, lower.tail = TRUE),
       col = "blue")
abline(v = c(t1, t2))

dev.new()
z <- seq(0,100,.01)
plot(ecdf(r2), main = "lambda = 1")
lines(z, pbeta(z/(1+z), ny, nx, lower.tail=TRUE), lwd = 2, col="red")
lines(z, ppearsonVI(q = z, a = ny, b = nx, location = 0, scale = 1, lower.tail = TRUE),
       col = "blue")
abline(v = t2)

pbeta(t1/(1+t1), nx, ny, lower.tail = FALSE)+
pbeta(t2/(1+t2), ny, nx, lower.tail = FALSE)

plot(ecdf(r))
points(z,pbeta(z/(z+1), nx, ny, lower.tail=TRUE), type="l", col="blue")

# empirical density (black) with theoretical overlaid (red)
qplot(x=r, geom="density")+
stat_function(fun = dbetaPrime,
  args = list(a=nx, b=ny), 
  color = "red")+
  theme_bw()


#
z <- seq(0,2,.01)
plot(z,beta(z/(1+z), nx, ny, lower.tail=TRUE), type="l", col="red")
plot(z,dbeta(z/(1+z), ny, nx), type="l", col="red")
abline(v = c(1/t, t))





nx <- 200
ny <- 300
lambda <- 1
x <- rexp(n=nx, rate=lambda)
y <- rexp(n=ny, rate=lambda)
t0 <- mean(x) / mean(y)
z <- c(x, y)
N <- length(z)
nIter <- 1e4
tStar <- rep(NA, length = N)
for (i in 1:nIter){
  piStar <- sample(1:N)
  xStar <- z[piStar[1:nx]]
  yStar <- z[piStar[(nx+1):N]]
  tStar[i] <- mean(xStar) / mean(yStar)
}

tStarMax <- pmax(tStar, 1/tStar)
t0Max <- max(t0, 1/t0)
mean(c(tStarMax, t0Max) >= t0Max)
mean(c(tStar, t0) >= t0) + mean(c(tStar, t0) <= 1/t0) 
mean(c(tStar, t0) <= t0) + mean(c(tStar, t0) >= 1/t0) 

hist(z)
dev.new()
hist(tStar)
plot(density(tStar))
abline(v = c(t0, 1/t0), col = "red", lwd = 2)

sd(z)
mean(z)

sd(tStar, na.rm = TRUE)
mean(tStar, na.rm = TRUE)

qplot(x=nx /ny * tStar, geom="density")+
stat_function(fun = dbetaPrime,
  args = list(a = nx, b = ny), 
  color = "red")+
  theme_bw()

qplot(x= tStar, geom="density")+
stat_function(fun = dpearsonVI,
  args = list(a = nx, b = ny, scale = ny/nx, location = 0), 
  color = "red")+
  theme_bw()+
  geom_vline(xintercept = t0)

qplot(x= tStar, geom="density")+
stat_function(fun = dpearsonVI,
  args = list(a = ny, b = nx, scale = nx/ny, location = 0), 
  color = "red")+
stat_function(fun = dpearsonVI,
  args = list(a = nx, b = ny, scale = ny/nx, location = 0), 
  color = "blue")+
  theme_bw()+
  geom_vline(xintercept = t0)


  ppearsonVI

mu <- nx / (ny - 1)
sigma <- sqrt(nx*(nx+ny-1)/((ny-2)*(ny-1)^2))
tStarNorm <- (tStar - mu) / sigma

qplot(x= tStar, geom="density")+
stat_function(fun = dnorm,
  args = list(mean = mu, sd = sigma), 
  color ="red")+
  theme_bw()

qplot(x= tStarNorm, geom="density")+
stat_function(fun = dt,
  args = list(df = 5), 
  color ="red")+
  theme_bw()

# beta prime distribution holds for large n, but not for small n
z <- seq(0,2,.01)
plot(ecdf(tStarMax))
lines(z, pbeta(z/(1+z), nx, ny, lower.tail=TRUE), col="red")