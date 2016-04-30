# Visual check to make sure beta prime distribution is correct

library(ggplot2)

dbetaPrime <- function(x,a,b){
  (x^(a-1)*(1+x)^(-a-b))/beta(a,b)
}

N <- 1e5
nx <- 100
ny <- 500
lambda <- 1

r <- rep(NA, N)
xbar <- rep(NA, N)
ybar <- rep(NA, N)

for (i in 1:N){
  # x <- rexp(n=nx, rate=lambda/nx)
  # y <- rexp(n=ny, rate=lambda/ny)
  
  x <- rexp(n=nx, rate=lambda)
  y <- rexp(n=ny, rate=lambda)
  
  # x <- rnbinom(n = 1000, size = 3, mu = 4)
  # y <- rnbinom(n = 1000, size = 3, mu = 4)

  # x <- rgamma(n=1, shape=nx, rate=lambda)
  # y <- rgamma(n=1, shape=ny, rate=lambda)
  xbar[i] <- mean(x)
  ybar[i] <- mean(y)
  r[i] <- (mean(x)*nx)/(mean(y)*ny)
  # r[i] <- x/y
}

qplot(x=r, geom="density")+
stat_function(fun = dbetaPrime,
  args = list(a=nx, b=ny), 
  color = "red")

z <- seq(0,2,.01)
plot(ecdf(r))
points(z,pbeta(z/(1+z), nx, ny, lower.tail=TRUE), type="l", col="red")