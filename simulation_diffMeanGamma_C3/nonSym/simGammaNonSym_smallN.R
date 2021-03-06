# Simulations comparing proposed algorithm to SAMC
# exponential data, unequal sample size

# arg <- commandArgs(TRUE)

library(doSNOW)
library(itertools)

nClust <- 6

# Nonsymmetric sample sizes ----------------------------------
nx <- c(20, 40, 60)
ny <- 100
lambdax <- 1

lambday <- c(2.5, 3)
params05 <- expand.grid(nx = nx, ny = ny, lambday = lambday, lambdax = lambdax, alpha = 0.5)

lambday <- c(1.5, 1.75)
params3 <- expand.grid(nx = nx, ny = ny, lambday = lambday, lambdax = lambdax, alpha = 3)

lambday <- c(1.25, 1.5)
params5 <- expand.grid(nx = nx, ny = ny, lambday = lambday, lambdax = lambdax, alpha = 5)

M <- 100 # M is number of repetitions for each scenario

params <- rbind(params05[rep(1:nrow(params05), each = M), ],
                params3[rep(1:nrow(params3), each = M), ],
                params5[rep(1:nrow(params5), each = M), ])

simFun <- function(sub){

  library(fastPerm)
  library(EXPERT)
  library(PearsonDS)
  library(gammaDist)

  # test stat with unequal variance for SAMC
  t.test.stat.uneqVar <- function (data.input) 
  {
      x <- data.input$x
      y <- data.input$y
      abs(t.test(x = x, y = y, var.equal = FALSE)$statistic)
  }

  sub$ptUneq <- NA
  sub$ptEq <- NA

  sub$gammaDiff <- NA
  sub$gammaDiffSaddle <- NA

  sub$pPred <- NA
  sub$pPredStud <- NA
  sub$mStop <- NA
  sub$mStopStud <- NA
  sub$EmStop <- NA

  sub$pAsymNorm <- NA
  sub$pAsymT <- NA
  # sub$pAsymNormStud <- NA
  # sub$pAsymTStud <- NA
  
  sub$pPerm <- NA
  sub$pPermAdj <- NA
  sub$pPermStud <- NA
  sub$pPermStudAdj <- NA
  
  sub$pExpert6 <- NA
  sub$timePred <- NA
  sub$timePredStud <- NA
  sub$timeExpert6 <- NA
  sub$maxErrorExpert6 <- NA

  Bmc <- 1e5

  for (i in 1:nrow(sub)){

    nx <- sub$nx[i]
    ny <- sub$ny[i]
    x <- rgamma(n = nx, shape = sub$alpha[i], rate = sub$lambdax[i])
    y <- rgamma(n = ny, shape = sub$alpha[i], rate = sub$lambday[i])

    # t-test
    sub$ptUneq[i] <- t.test(x, y, var.equal = FALSE)$p.value
    sub$ptEq[i] <- t.test(x, y, var.equal = TRUE)$p.value

    #
    mle <- gammaFit(c(x,y))$MLE

    # difference of gammas
    t0 <- abs(mean(x) - mean(y))
    try({
      sub$gammaDiff[i] <- 1-pgammaDif(t0, nx = nx, ny = ny, alpha = sub$alpha[i], lambda = mle[2]) +
                            pgammaDif(-t0, nx = nx, ny = ny, alpha = sub$alpha[i], lambda = mle[2])
    })
    try({
      sub$gammaDiffSaddle[i] <- 
          pgammaDifSaddle(t0, nx = nx, ny = ny, alpha = sub$alpha[i], lambda = mle[2], lower.tail = FALSE) +
          pgammaDifSaddle(-t0, nx = nx, ny = ny, alpha = sub$alpha[i], lambda = mle[2], lower.tail = TRUE)
    })

    # fastPerm
    try({
      sub$timePred[i] <- system.time(fp <- 
        fastPerm(x, y, testStat = diffMean))[3]
      sub$pPred[i] <- fp$pPred
      sub$mStop[i] <- fp$mStop
    sub$EmStop[i] <- mStopDiffMean(x,y)
    })

    try({
      sub$timePredStud[i] <- system.time(fpStud <- 
        fastPerm(x, y, testStat = diffMeanStudent))[3]
      sub$pPredStud[i] <- fpStud$pPred
      sub$mStopStud[i] <- fpStud$mStop
    })

    # asymptotic method
    fpAsym <- fastPermAsym(x,y, testStat=diffMean)
    sub$pAsymNorm[i] <- fpAsym$pNorm
    sub$pAsymT[i] <- fpAsym$pT

    # Monte Carlo
    z <- c(x, y)
    N <- length(z)
    Tstud <- rep(NA, Bmc)
    T <- rep(NA, Bmc)

    t0 <- mean(x) - mean(y)
    tStud0 <- t0 / sqrt(var(x) / nx + var(y) / ny)
    for (b in 1:Bmc) {
      star <- sample(1:N, replace = FALSE)
      xStar <- z[star[1:nx]]
      yStar <- z[star[(nx+1):N]]
      T[b] <- mean(xStar) - mean(yStar)
      Tstud[b] <- (mean(xStar) - mean(yStar)) / sqrt(var(xStar) / nx + var(yStar) / ny)
    }
    sub$pPerm[i] <- mean(abs(T) >= abs(t0))
    sub$pPermAdj[i] <- mean(abs(c(T, t0)) >= abs(t0))
    sub$pPermStud[i] <- mean(abs(Tstud) >= abs(tStud0))
    sub$pPermStudAdj[i] <- mean(abs(c(Tstud, tStud0)) >= abs(tStud0))

    # # using EXPERT package
    # data.input <- list(x=x, y=y)
    # t.obs <- t.test.stat.uneqVar(data.input)

    # try({
    #   sub$timeExpert6[i] <- 
    #     system.time(res6 <- SAMC.adapt(data.input, t.obs, t.start=0,
    #       n.iter.1=5e4, n.iter.2=1e6, 
    #       n.region.1=101, n.region.2=301,
    #       prop.change=0.05, gain.factor.t0=1000,
    #       fun.test.statistic=t.test.stat.uneqVar,
    #       fun.proposal=proposal.permute.vector)
    #     )[3]
    #   # The estimated p-value
    #   sub$pExpert6[i] <- res6$p.value
    #   sub$maxErrorExpert6[i] <- maxError(res6)
    #   })

  }

  return(sub)  
}

cl <- makeCluster(nClust, type="SOCK")
registerDoSNOW(cl)

blocks <- isplitIndices(nrow(params), chunks=nClust)

nonSymResults <- foreach(j=blocks, .combine=rbind) %dopar% {
                   return(simFun(params[j,]))
                 }

stopCluster(cl)

save(nonSymResults, file ="nonSymResultsGammaDiff_smallN.RData")
