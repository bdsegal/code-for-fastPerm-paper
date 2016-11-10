# Section 5, Figure 2
# Simulated comparison with t-test

library(doSNOW)
library(itertools)

# Symmetric sample sizes ----------------------------------
n <- c(100, 500, 1000)
nLen <- length(n)
nClust <- 6

mux <- c(0.75, 1)
muLen <- length(mux)

M <- 100 # M is number of repetitions for each scenario

params <- data.frame(n=rep(n, each=muLen*M),
    mux=rep(mux, each=M)
)

simFun <- function(sub){

  library(fastPerm)
  library(EXPERT)
  # library(mcc)

  sub$pt <- NA
  sub$pPred <- NA
  sub$pAsymNorm <- NA
  sub$pAsymT <- NA
  sub$pExpert6 <- NA

  sub$timePred <- NA
  sub$timeExpert6 <- NA

  sub$mStop <- NA
  sub$EmStop <- NA

  sub$maxErrorExpert6 <- NA

  for (i in 1:nrow(sub)){

    n <- sub$n[i]
    x <- rnorm(n = n, mean = sub$mux[i], sd=1)
    y <- rnorm(n = n, mean = 0, sd=1)

    sub$pt[i] <- t.test(x, y, var.equal = TRUE)$p.value

    sub$timePred[i] <- system.time(fp <- 
      fastPerm(x, y, testStat = diffMean))[3]
    sub$pPred[i] <- fp$pPred
    sub$mStop[i] <- fp$mStop

    sub$EmStop[i] <- mStopDiffMean(x,y)

    fpAsym <- fastPermAsym(x,y, testStat=diffMean)
    sub$pAsymNorm[i] <- fpAsym$pNorm
    sub$pAsymT[i] <- fpAsym$pT
    
    # # using EXPERT package
    data.input<-list(x=x, y=y)
    t.obs<-t.test.statistic(data.input)

    try({
      sub$timeExpert6[i] <- 
        system.time(res6 <- SAMC.adapt(data.input, t.obs, t.start=0,
          n.iter.1=5e4, n.iter.2=1e6, 
          n.region.1=101, n.region.2=301,
          prop.change=0.05, gain.factor.t0=1000,
          fun.test.statistic=t.test.statistic,
          fun.proposal=proposal.permute.vector)
        )[3]
      # The estimated p-value
      sub$pExpert6[i] <- res6$p.value
      sub$maxErrorExpert6[i] <- maxError(res6)
      })

    # using MCC

  }

  return(sub)  
}

cl <- makeCluster(nClust, type="SOCK")
registerDoSNOW(cl)

blocks <- isplitIndices(nrow(params), chunks=nClust)

system.time(symResults <- foreach(j=blocks, .combine=rbind) %dopar% {
                   return(simFun(params[j,]))
                 })

stopCluster(cl)

save(symResults, file ="symResultsDiff_parallel.RData")







# using mcc package
library(mcc)
m=300
n=100

X <- matrix(rnorm(n = n*m, mean = 1, sd = 1), 
            nrow = m, ncol = n)
Y <- matrix(rnorm(n = n*m, mean = 0, sd = 1), 
            nrow = m, ncol = n)
Z <- cbind(X, Y)

contrast <- c(rep(1/n, n), rep(-1/n, n))

output <- getbetap.A(getAmoment(Z, contrast), A=NULL, fix.obs = FALSE)
output2 <- getbetap.A.2(Z, contrast)

pt <- rep(NA, nrow(X))
for (i in 1:nrow(X)){
  pt[i] <- t.test(X[i, 1:n], Y[i, 1:n])$p.value
}

plot(y = log(pt, 10), x = log(output$twosidedp, 10))
abline(a = 0, b = 1, col = "red", lty = 2)


hist(output$twosidedp)

getbetap.A

getAmoment.list <- getAmoment(Z, contrast)

# alter function to use lfactorial to remove errors
getbetap.A <- function (getAmoment.list, A = NULL, fix.obs = F) 
{
    x = getAmoment.list$x
    y = getAmoment.list$y
    n = getAmoment.list$n
    z = getAmoment.list$z
    if (length(A) == 0) {
        A = getAmoment.list$A
    } else {
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
    # lowest.alpha = 2/factorial(effective.n)
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