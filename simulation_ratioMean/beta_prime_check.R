# Visual check to make sure beta prime distribution is correct

library(ggplot2)

# beta prime density (also called the inverted beta)
dbetaPrime <- function(x,a,b){
  (x^(a-1)*(1+x)^(-a-b))/beta(a,b)
}

# group sample sizes and common rate parameter
nx <- 250
ny <- 500
lambda <- 5

# number of simulated ratio statistics; r holds the ratio of sums
Niter <- 1e5
r <- rep(NA, Niter)

for (i in 1:Niter){

  x <- rexp(n=nx, rate=lambda)
  y <- rexp(n=ny, rate=lambda)
  
  r[i] <- sum(x)/sum(y)
}

# empirical cdf (blakck) with theoretical overlaid (red)
z <- seq(0,2,.01)
plot(ecdf(r))
points(z,pbeta(z/(1+z), nx, ny, lower.tail=TRUE), type="l", col="red")

# empirical density (black) with theoretical overlaid (red)
qplot(x=r, geom="density")+
stat_function(fun = dbetaPrime,
  args = list(a=nx, b=ny), 
  color = "red")+
  theme_bw()