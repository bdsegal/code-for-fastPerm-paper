# note: mcc isn't available for the version of R on the cluster
# so I'm just running this sequentially on a desktop
# Some of this code is not efficient, but I'm just running it once
# so I'm not too worried about it.

library(mcc)

MCmulti <- function(X, Y, M, Bmc = 1e5) {
  # simultaneously conduct studentized permutation test with simple MC

  Z <- cbind(X, Y)
  nx <- ncol(X)
  ny <- ncol(Y)
  N <- nx + ny

  varx <- apply(X, 1, var)
  vary <- apply(Y, 1, var)
  xbar <- apply(X, 1, mean)
  ybar <- apply(Y, 1, mean)

  T0 <- abs(xbar - ybar) / sqrt(varx / nx + vary / ny)
  T <- matrix(NA, nrow = M, ncol = Bmc)

  for(b in 1:Bmc) {
    star <- sample(1:N)
    Zstar <- Z[, star]
    Xstar <- Zstar[, 1:nx]
    Ystar <- Zstar[, (nx+1):N]

    varxstar <- apply(Xstar, 1, var)
    varystar <- apply(Ystar, 1, var)
    xbarstar <- apply(Xstar, 1, mean)
    ybarstar <- apply(Ystar, 1, mean)

    T[, b] <-  abs(xbarstar - ybarstar) / sqrt(varxstar / nx + varystar / ny)
  }

  TAll <- cbind(T0, T)
  pMC <- apply(TAll, 1, function(x) {mean(x[-1] >= x[1])})
  return(pMC)
}

# Small sample size -----------------------------------------------------------

M <- 100 # number of replicates per parameter combo

# Symmetric sample sizes ------------------------
n <- c(20, 40, 60)
lambdax <- 1

lambday <- c(2.5, 3)
params05 <- expand.grid(n = n, lambday = lambday, lambdax = lambdax, alpha = 0.5)
lambday <- c(1.5, 1.75)
params3 <- expand.grid(n = n, lambday = lambday, lambdax = lambdax, alpha = 3)
lambday <- c(1.25, 1.5)
params5 <- expand.grid(n = n, lambday = lambday, lambdax = lambdax, alpha = 5)

params <- rbind(params05, params3, params5)

twosidedp <- NULL
pdouble <- NULL
pt <- NULL
pMC <- NULL

for (i in 1:nrow(params)) {
  print(i)

  ny <- params$n[i]
  nx <- params$n[i]
  N <- nx + ny

  X <- matrix(rgamma(n = nx*M, shape = params$alpha[i], rate = params$lambdax[i]), 
              ncol = nx, nrow = M)
  Y <- matrix(rgamma(n = ny*M, shape = params$alpha[i], rate = params$lambday[i]),
              ncol = ny, nrow = M)

  Z <- cbind(X, Y)
  y <- c(rep(1/nx, nx), -rep(1/ny, ny))

  output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)
  pdouble <- c(pdouble, output$pdouble)
  twosidedp <- c(twosidedp, output$twosidedp)

  pMC <- c(pMC, MCmulti(X, Y, M))

  pt <- c(pt, apply(Z, 1, function(x){
    x1 <- x[1:nx]
    x2 <- x[(nx+1):N]
    t.test(x1, x2, var.equal = FALSE)$p.value
    }))

}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt, pMC = pMC)
mccOut$n <- rep(params$n, each = M)
mccOut$alpha <- rep(params$alpha, each = M)

save(mccOut, file = "mccOut_symGamma_smallN.Rdata")

# non-symmetric sample sizes --------------------
nx <- c(20, 40, 60)
ny <- 100
lambdax <- 1

lambday <- c(2.5, 3)
params05 <- expand.grid(nx = nx, ny = ny, lambday = lambday, lambdax = lambdax, alpha = 0.5)
lambday <- c(1.5, 1.75)
params3 <- expand.grid(nx = nx, ny = ny, lambday = lambday, lambdax = lambdax, alpha = 3)
lambday <- c(1.25, 1.5)
params5 <- expand.grid(nx = nx, ny = ny, lambday = lambday, lambdax = lambdax, alpha = 5)

params <- rbind(params05, params3, params5)

twosidedp <- NULL
pdouble <- NULL
pt <- NULL
pMC <- NULL

for (i in 1:nrow(params)) {
  print(i)

  ny <- params$nx[i]
  nx <- params$ny[i]
  N <- nx + ny

  X <- matrix(rgamma(n = nx*M, shape = params$alpha[i], rate = params$lambdax[i]), 
              ncol = nx, nrow = M)
  Y <- matrix(rgamma(n = ny*M, shape = params$alpha[i], rate = params$lambday[i]),
              ncol = ny, nrow = M)

  Z <- cbind(X, Y)
  y <- c(rep(1/nx, nx), -rep(1/ny, ny))

  output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)

  pdouble <- c(pdouble, output$pdouble)
  twosidedp <- c(twosidedp, output$twosidedp)

  pMC <- c(pMC, MCmulti(X, Y, M))

  pt <- c(pt, apply(Z, 1, function(x){
    x1 <- x[1:nx]
    x2 <- x[(nx+1):N]
    t.test(x1, x2, var.equal = FALSE)$p.value
    }))

}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt, pMC = pMC)
mccOut$nx <- rep(params$nx, each = M)
mccOut$alpha <- rep(params$alpha, each = M)

save(mccOut, file = "mccOut_nonSymGamma_smallN.Rdata")


# null ------------------------------------------------------------------------

M <- 1000 # number of replicates per parameter combo

# Symmetric sample sizes ------------------------
n <- c(20, 40, 60)
alpha <- c(0.5, 3, 5)
lambda <- c(1, 5)
params <- expand.grid(n = n, lambda = lambda, alpha = alpha)

twosidedp <- NULL
pdouble <- NULL
pt <- NULL
pMC <- NULL

for (i in 1:nrow(params)) {
  print(i)

  ny <- params$n[i]
  nx <- params$n[i]
  N <- nx + ny

  X <- matrix(rgamma(n = nx*M, shape = params$alpha[i], rate = params$lambda[i]), 
              ncol = nx, nrow = M)
  Y <- matrix(rgamma(n = ny*M, shape = params$alpha[i], rate = params$lambda[i]),
              ncol = ny, nrow = M)

  Z <- cbind(X, Y)
  y <- c(rep(1/nx, nx), -rep(1/ny, ny))

  output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)
  pdouble <- c(pdouble, output$pdouble)
  twosidedp <- c(twosidedp, output$twosidedp)

  pMC <- c(pMC, MCmulti(X, Y, M))

  pt <- c(pt, apply(Z, 1, function(x){
    x1 <- x[1:nx]
    x2 <- x[(nx+1):N]
    t.test(x1, x2, var.equal = FALSE)$p.value
    }))

}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt, pMC = pMC)
mccOut$n <- rep(params$n, each = M)
mccOut$alpha <- rep(params$alpha, each = M)

save(mccOut, file = "mccOut_symGamma_null.Rdata")

# non-symmetric sample sizes --------------------
nx <- c(20, 40, 60)
ny <- 100
alpha <- c(0.5, 3, 5)
lambda <- c(1, 5)
params <- expand.grid(nx = nx, ny = ny, lambda = lambda, alpha = alpha)

twosidedp <- NULL
pdouble <- NULL
pt <- NULL
pMC <- NULL

for (i in 1:nrow(params)) {
  print(i)
  
  ny <- params$nx[i]
  nx <- params$ny[i]
  N <- nx + ny

  X <- matrix(rgamma(n = nx*M, shape = params$alpha[i], rate = params$lambda[i]), 
              ncol = nx, nrow = M)
  Y <- matrix(rgamma(n = ny*M, shape = params$alpha[i], rate = params$lambda[i]),
              ncol = ny, nrow = M)

  Z <- cbind(X, Y)
  y <- c(rep(1/nx, nx), -rep(1/ny, ny))

  output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)

  pdouble <- c(pdouble, output$pdouble)
  twosidedp <- c(twosidedp, output$twosidedp)

  pMC <- c(pMC, MCmulti(X, Y, M))

  pt <- c(pt, apply(Z, 1, function(x){
    x1 <- x[1:nx]
    x2 <- x[(nx+1):N]
    t.test(x1, x2, var.equal = FALSE)$p.value
    }))

}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt, pMC = pMC)
mccOut$nx <- rep(params$nx, each = M)
mccOut$alpha <- rep(params$alpha, each = M)

save(mccOut, file = "mccOut_nonSymGamma_null.Rdata")
