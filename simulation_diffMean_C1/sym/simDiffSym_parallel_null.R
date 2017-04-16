# Section 5, Figure 2
# Simulated comparison with t-test

library(doSNOW)
library(itertools)

# Symmetric sample sizes ----------------------------------
n <- c(20, 40, 60)
nLen <- length(n)
nClust <- 6

mux <- 0
muLen <- length(mux)

M <- 1000 # M is number of repetitions for each scenario

params <- data.frame(n=rep(n, each=muLen*M),
    mux=rep(mux, each=M)
)

simFun <- function(sub){

  library(fastPerm)
  library(EXPERT)

  sub$pt <- NA
  sub$pPred <- NA
  sub$pAsymNorm <- NA
  sub$pAsymT <- NA
  sub$pExpert6 <- NA
  sub$pPerm <- NA
  sub$pPermAdj <- NA

  sub$timePred <- NA
  sub$timeExpert6 <- NA

  sub$mStop <- NA
  sub$EmStop <- NA

  sub$maxErrorExpert6 <- NA

  nIter <- 1e5

  for (i in 1:nrow(sub)){

    n <- sub$n[i]
    x <- rnorm(n = n, mean = sub$mux[i], sd=1)
    y <- rnorm(n = n, mean = 0, sd=1)

    sub$pt[i] <- t.test(x, y, var.equal = TRUE)$p.value

    # MC simulation
    z <- c(x, y)
    N <- length(z)
    tStar <- rep(NA, length = N)
    t0 <- abs(mean(x) - mean(y))
    for (j in 1:nIter){
      piStar <- sample(1:N)
      xStar <- z[piStar[1:n]]
      yStar <- z[piStar[(n+1):N]]
      tStar[j] <- abs(mean(xStar) - mean(yStar))
    }
    sub$pPerm[i] <- mean(tStar >= t0)
    sub$pPermAdj[i] <- mean(c(tStar, t0) >= t0)


    try({
      sub$timePred[i] <- system.time(fp <- 
        fastPerm(x, y, testStat = diffMean))[3]
      sub$pPred[i] <- fp$pPred
      sub$mStop[i] <- fp$mStop
    })

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

save(symResults, file ="symResultsDiff_parallel_null.RData")


