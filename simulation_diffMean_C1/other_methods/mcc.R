# Section 5, Figure 2
# Simulated comparison with t-test
library(ggplot2)
library(reshape2)
library(mcc)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# Symmetric sample sizes ------------------------------------------------------
mux <- c(0.75, 1)
muLen <- length(mux)
n <- c(100, 500, 1000)
nLen <- length(n)
M <- 100 # M is number of repetitions for each scenario

twosidedp <- NULL
pdouble <- NULL
pt <- NULL

for (i in 1:nLen) {
  ny <- n[i]
  nx <- n[i]

  for(j in 1:muLen) {
    X <- matrix(rnorm(n = nx*M, mean = mux[j], sd = 1), nrow = M, ncol = nx)
    Y <- matrix(rnorm(n = ny*M, mean = 0, sd = 1), nrow = M, ncol = ny)
    Z <- cbind(X, Y)
    y <- c(rep(1/nx, nx), -rep(1/ny, ny))

    output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)

    pdouble <- c(pdouble, output$pdouble)
    twosidedp <- c(twosidedp, output$twosidedp)

    pt <- c(pt, apply(Z, 1, function(x){
      x1 <- x[1:nx]
      x2 <- x[(nx+1):(nx+ny)]
      t.test(x1, x2, var.equal = TRUE)$p.value
      }))

  }
}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt)
mccOut$n <- rep(n, each = muLen * M)
mccOut$mux <- rep(mux, each = M)

mccOutM <- melt(mccOut, id.vars = c("n", "mux", "pt"))
levels(mccOutM$variable) <- c("Two sided", "Double")

dev.new(width=8, height=5)
ggplot(aes(x=log(pt,10), y=log(value, 10), color=as.factor(n),
    shape=as.factor(n)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[t],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"mcc_sym.png"))
ggsave(file.path(paperPath,"mcc_sym.tiff"))
ggsave(file.path(paperPath,"mcc_sym.pdf"))



# non-symmetric sample sizes ------------------------------------------------------
mux <- c(0.75, 1)
muLen <- length(mux)
nxVec <- c(50, 200, 350)
nLen <- length(nxVec)
M <- 100 # M is number of repetitions for each scenario

twosidedp <- NULL
pdouble <- NULL
pt <- NULL

for (i in 1:nLen) {
  ny <- 500
  nx <- nxVec[i]

  for(j in 1:muLen) {
    X <- matrix(rnorm(n = nx*M, mean = mux[j], sd = 1), nrow = M, ncol = nx)
    Y <- matrix(rnorm(n = ny*M, mean = 0, sd = 1), nrow = M, ncol = ny)
    Z <- cbind(X, Y)
    y <- c(rep(1/nx, nx), -rep(1/ny, ny))

    output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)

    pdouble <- c(pdouble, output$pdouble)
    twosidedp <- c(twosidedp, output$twosidedp)

    pt <- c(pt, apply(Z, 1, function(x){
      x1 <- x[1:nx]
      x2 <- x[(nx+1):(nx+ny)]
      t.test(x1, x2, var.equal = TRUE)$p.value
      }))

  }
}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt)
mccOut$nx <- rep(nxVec, each = muLen * M)
mccOut$mux <- rep(mux, each = M)

mccOutM <- melt(mccOut, id.vars = c("nx", "mux", "pt"))
levels(mccOutM$variable) <- c("Two sided", "Double")

dev.new(width=8, height=5)
ggplot(aes(x=log(pt,10), y=log(value, 10), color=as.factor(nx),
    shape=as.factor(nx)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[t],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"mcc_nonSym.png"))
ggsave(file.path(paperPath,"mcc_nonSym.tiff"))
ggsave(file.path(paperPath,"mcc_nonSym.pdf"))


# Small sample size
# Symmetric sample sizes ------------------------------------------------------
mux <- c(2, 3)
muLen <- length(mux)
n <- c(20, 40, 60)
nLen <- length(n)
M <- 100 # M is number of repetitions for each scenario

twosidedp <- NULL
pdouble <- NULL
pt <- NULL

for (i in 1:nLen) {
  ny <- n[i]
  nx <- n[i]

  for(j in 1:muLen) {
    X <- matrix(rnorm(n = nx*M, mean = mux[j], sd = 1), nrow = M, ncol = nx)
    Y <- matrix(rnorm(n = ny*M, mean = 0, sd = 1), nrow = M, ncol = ny)
    Z <- cbind(X, Y)
    y <- c(rep(1/nx, nx), -rep(1/ny, ny))

    output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)

    pdouble <- c(pdouble, output$pdouble)
    twosidedp <- c(twosidedp, output$twosidedp)

    pt <- c(pt, apply(Z, 1, function(x){
      x1 <- x[1:nx]
      x2 <- x[(nx+1):(nx+ny)]
      t.test(x1, x2, var.equal = TRUE)$p.value
      }))

  }
}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt)
mccOut$n <- rep(n, each = muLen * M)
mccOut$mux <- rep(mux, each = M)

mccOutM <- melt(mccOut, id.vars = c("n", "mux", "pt"))
levels(mccOutM$variable) <- c("Two sided", "Double")

dev.new(width=8, height=5)
ggplot(aes(x=log(pt,10), y=log(value, 10), color=as.factor(n),
    shape=as.factor(n)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[t],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")
ggsave(file.path(paperPath,"mcc_sym_smallN.png"))
ggsave(file.path(paperPath,"mcc_sym_smallN.tiff"))
ggsave(file.path(paperPath,"mcc_sym_smallN.pdf"))


# non-symmetric sample sizes ------------------------------------------------------
mux <- c(2, 3)
muLen <- length(mux)
nxVec <- c(20, 40, 60)
nLen <- length(nxVec)
M <- 100 # M is number of repetitions for each scenario

twosidedp <- NULL
pdouble <- NULL
pt <- NULL

for (i in 1:nLen) {
  ny <- 100
  nx <- nxVec[i]

  for(j in 1:muLen) {
    X <- matrix(rnorm(n = nx*M, mean = mux[j], sd = 1), nrow = M, ncol = nx)
    Y <- matrix(rnorm(n = ny*M, mean = 0, sd = 1), nrow = M, ncol = ny)
    Z <- cbind(X, Y)
    y <- c(rep(1/nx, nx), -rep(1/ny, ny))

    output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)

    pdouble <- c(pdouble, output$pdouble)
    twosidedp <- c(twosidedp, output$twosidedp)

    pt <- c(pt, apply(Z, 1, function(x){
      x1 <- x[1:nx]
      x2 <- x[(nx+1):(nx+ny)]
      t.test(x1, x2, var.equal = TRUE)$p.value
      }))

  }
}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt)
mccOut$nx <- rep(nxVec, each = muLen * M)
mccOut$mux <- rep(mux, each = M)

mccOutM <- melt(mccOut, id.vars = c("nx", "mux", "pt"))
levels(mccOutM$variable) <- c("Two sided", "Double")

dev.new(width=8, height=5)
ggplot(aes(x=log(pt,10), y=log(value, 10), color=as.factor(nx),
    shape=as.factor(nx)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[t],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))
ggsave(file.path(paperPath,"mcc_nonSym_smallN.png"))
ggsave(file.path(paperPath,"mcc_nonSym_smallN.tiff"))
ggsave(file.path(paperPath,"mcc_nonSym_smallN.pdf"))



# null

# Small sample size
# Symmetric sample sizes ------------------------------------------------------
n <- c(20, 40, 60)
nLen <- length(n)
M <- 100 # M is number of repetitions for each scenario

twosidedp <- NULL
pdouble <- NULL
pt <- NULL

for (i in 1:nLen) {
  ny <- n[i]
  nx <- n[i]

    X <- matrix(rnorm(n = nx*M, mean = 0, sd = 1), nrow = M, ncol = nx)
    Y <- matrix(rnorm(n = ny*M, mean = 0, sd = 1), nrow = M, ncol = ny)
    Z <- cbind(X, Y)
    y <- c(rep(1/nx, nx), -rep(1/ny, ny))

    output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)

    pdouble <- c(pdouble, output$pdouble)
    twosidedp <- c(twosidedp, output$twosidedp)

    pt <- c(pt, apply(Z, 1, function(x){
      x1 <- x[1:nx]
      x2 <- x[(nx+1):(nx+ny)]
      t.test(x1, x2, var.equal = TRUE)$p.value
      }))

}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt)
mccOut$n <- rep(n, each = M)

mccOutM <- melt(mccOut, id.vars = c("n", "pt"))
levels(mccOutM$variable) <- c("Two sided", "Double")

dev.new(width=8, height=5)
ggplot(aes(x=pt, y=value, color=as.factor(n),
    shape=as.factor(n)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(p[t]),
    y="p")+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
 theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"mcc_sym_null.png"))
ggsave(file.path(paperPath,"mcc_sym_null.tiff"))
ggsave(file.path(paperPath,"mcc_sym_null.pdf"))


# non-symmetric sample sizes ------------------------------------------------------
nxVec <- c(20, 40, 60)
nLen <- length(nxVec)
M <- 100 # M is number of repetitions for each scenario

twosidedp <- NULL
pdouble <- NULL
pt <- NULL

for (i in 1:nLen) {
  ny <- 100
  nx <- nxVec[i]

    X <- matrix(rnorm(n = nx*M, mean = 0, sd = 1), nrow = M, ncol = nx)
    Y <- matrix(rnorm(n = ny*M, mean = 0, sd = 1), nrow = M, ncol = ny)
    Z <- cbind(X, Y)
    y <- c(rep(1/nx, nx), -rep(1/ny, ny))

    output=getbetap.A(getAmoment(x=Z,y=y),A=NULL,fix.obs=FALSE)

    pdouble <- c(pdouble, output$pdouble)
    twosidedp <- c(twosidedp, output$twosidedp)

    pt <- c(pt, apply(Z, 1, function(x){
      x1 <- x[1:nx]
      x2 <- x[(nx+1):(nx+ny)]
      t.test(x1, x2, var.equal = TRUE)$p.value
      }))

}

mccOut <- data.frame(twosidedp = twosidedp, pdouble = pdouble, pt = pt)
mccOut$nx <- rep(nxVec, each = M)

mccOutM <- melt(mccOut, id.vars = c("nx", "pt"))
levels(mccOutM$variable) <- c("Two sided", "Double")

dev.new(width=8, height=5)
ggplot(aes(x=pt, y=value, color=as.factor(nx),
    shape=as.factor(nx)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(p[t]),
    y="p")+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"mcc_nonSym_null.png"))
ggsave(file.path(paperPath,"mcc_nonSym_null.tiff"))
ggsave(file.path(paperPath,"mcc_nonSym_null.pdf"))



#------------------------------------------------------------------------------
# fix getbetap.A function so that it doesn't give errors -- 
# replaced factorial(effective.n) with exp(lfactorial(effective.n))
# make note to authors
getbetap.A <- function (getAmoment.list, A = NULL, fix.obs = F) 
{
    x = getAmoment.list$x
    y = getAmoment.list$y
    n = getAmoment.list$n
    z = getAmoment.list$z
    if (length(A) == 0) {
        A = getAmoment.list$A
    }
    else {
        A = A - getAmoment.list$mu
    }
    EA = getAmoment.list$EA
    EA2 = getAmoment.list$EA2
    EA3 = getAmoment.list$EA3
    EA4 = getAmoment.list$EA4
    V = EA2 - EA^2
    s = (EA3 - 3 * EA * V - 3 * EA^3)/V^(3/2)
    k = EA4/V^2 - 3
    effective.n = -1 - 6/(-abs(k))
    lowest.alpha = 2/exp(lfactorial(effective.n))
    lowest.alpha[EA2 == 0] = 1
    r = A/sqrt(V * (n - 1))
    alpha1 = (3 * k + 36 * s * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * 
        k * s^2 + 96 * k + 36 * s^4 - 180 * s^2))^(1/2) - 18 * 
        s^3 * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * k * s^2 + 96 * 
        k + 36 * s^4 - 180 * s^2))^(1/2) - 3 * s^2 + 3 * k^2 * 
        s * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * k * s^2 + 96 * 
        k + 36 * s^4 - 180 * s^2))^(1/2) - 3 * k * s^3 * (-1/(-k^2 * 
        s^2 + 32 * k^2 - 84 * k * s^2 + 96 * k + 36 * s^4 - 180 * 
        s^2))^(1/2) + 24 * k * s * (-1/(-k^2 * s^2 + 32 * k^2 - 
        84 * k * s^2 + 96 * k + 36 * s^4 - 180 * s^2))^(1/2) + 
        6)/(2 * k - 3 * s^2) - (-6 * s^2 + 6 * k + 12)/(2 * k - 
        3 * s^2)
    alpha2 = (3 * k - 36 * s * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * 
        k * s^2 + 96 * k + 36 * s^4 - 180 * s^2))^(1/2) + 18 * 
        s^3 * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * k * s^2 + 96 * 
        k + 36 * s^4 - 180 * s^2))^(1/2) - 3 * s^2 - 3 * k^2 * 
        s * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * k * s^2 + 96 * 
        k + 36 * s^4 - 180 * s^2))^(1/2) + 3 * k * s^3 * (-1/(-k^2 * 
        s^2 + 32 * k^2 - 84 * k * s^2 + 96 * k + 36 * s^4 - 180 * 
        s^2))^(1/2) - 24 * k * s * (-1/(-k^2 * s^2 + 32 * k^2 - 
        84 * k * s^2 + 96 * k + 36 * s^4 - 180 * s^2))^(1/2) + 
        6)/(2 * k - 3 * s^2) - (-6 * s^2 + 6 * k + 12)/(2 * k - 
        3 * s^2)
    beta1 = -(3 * k + 36 * s * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * 
        k * s^2 + 96 * k + 36 * s^4 - 180 * s^2))^(1/2) - 18 * 
        s^3 * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * k * s^2 + 96 * 
        k + 36 * s^4 - 180 * s^2))^(1/2) - 3 * s^2 + 3 * k^2 * 
        s * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * k * s^2 + 96 * 
        k + 36 * s^4 - 180 * s^2))^(1/2) - 3 * k * s^3 * (-1/(-k^2 * 
        s^2 + 32 * k^2 - 84 * k * s^2 + 96 * k + 36 * s^4 - 180 * 
        s^2))^(1/2) + 24 * k * s * (-1/(-k^2 * s^2 + 32 * k^2 - 
        84 * k * s^2 + 96 * k + 36 * s^4 - 180 * s^2))^(1/2) + 
        6)/(2 * k - 3 * s^2)
    beta2 = -(3 * k - 36 * s * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * 
        k * s^2 + 96 * k + 36 * s^4 - 180 * s^2))^(1/2) + 18 * 
        s^3 * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * k * s^2 + 96 * 
        k + 36 * s^4 - 180 * s^2))^(1/2) - 3 * s^2 - 3 * k^2 * 
        s * (-1/(-k^2 * s^2 + 32 * k^2 - 84 * k * s^2 + 96 * 
        k + 36 * s^4 - 180 * s^2))^(1/2) + 3 * k * s^3 * (-1/(-k^2 * 
        s^2 + 32 * k^2 - 84 * k * s^2 + 96 * k + 36 * s^4 - 180 * 
        s^2))^(1/2) - 24 * k * s * (-1/(-k^2 * s^2 + 32 * k^2 - 
        84 * k * s^2 + 96 * k + 36 * s^4 - 180 * s^2))^(1/2) + 
        6)/(2 * k - 3 * s^2)
    alpha = alpha1
    beta = beta1
    which.negative = grep(T, (alpha <= 0) | (beta <= 0))
    alpha[which.negative] = alpha2[which.negative]
    beta[which.negative] = beta2[which.negative]
    beta.mean = alpha/(alpha + beta)
    beta.var = (alpha * beta)/((alpha + beta)^2 * (alpha + beta + 
        1))
    c0 = beta.mean
    c1 = sqrt(beta.var * (n - 1))
    rprime = c0 + c1 * r
    rprime.high = c0 + c1 * abs(r)
    rprime.low = c0 - c1 * abs(r)
    twosidedp = pbeta(rprime.high, alpha, beta, lower.tail = F) + 
        pbeta(rprime.low, alpha, beta)
    rightp = pbeta(rprime, alpha, beta, lower.tail = F)
    leftp = pbeta(rprime, alpha, beta)
    doublep = 2 * apply(cbind(rightp, leftp), 1, min)
    t = A/sqrt(EA2)
    r2 = t^2/(n - 2 + t^2)
    pt = pbeta(r2, 1/2, 0.5 * (n - 2), lower.tail = F)
    if (fix.obs == T) {
        which = unique(c(grep(T, rightp == 0), grep(T, leftp == 
            0), grep(T, twosidedp == 0)))
        if (length(which) > 1) {
            x0 = x[which, ]
            result0 = getbetap.A.2(x0, y, z = z)
            leftp[which] = result0$leftp
            rightp[which] = result0$rightp
            twosidedp[which] = result0$twosidedp
            doublep[which] = 2 * apply(cbind(rightp[which], leftp[which]), 
                1, min)
        }
    }
    chebyshev.p = EA2/A^2
    fail.chebyshev = grep(T, (chebyshev.p < twosidedp))
    twosidedp[fail.chebyshev] = chebyshev.p[fail.chebyshev]
    doublep[fail.chebyshev] = chebyshev.p[fail.chebyshev]
    leftp[fail.chebyshev] = chebyshev.p[fail.chebyshev]/2
    rightp[fail.chebyshev] = chebyshev.p[fail.chebyshev]/2
    if (fix.obs == T) {
        which = unique(c(grep(T, rightp == 0), grep(T, leftp == 
            0), grep(T, twosidedp == 0)))
        if (length(which) == 1) {
            current = cor.test(x[which, ], y, method = "spearman")
            doublep[which] = twosidedp[which] = current$p.value
        }
        if (length(which) > 1) {
            for (j in (1:length(which))) {
                current = cor.test(x[which[j], ], y, method = "spearman")
                doublep[which[j]] = twosidedp[which[j]] = current$p.value
            }
        }
    }
    return(list(twosidedp = twosidedp, rightp = rightp, leftp = leftp, 
        pdouble = doublep, chebyshev.p = chebyshev.p, pt = pt, 
        lowest.alpha = lowest.alpha))
}