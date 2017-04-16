# Section 5, Figure 2
# Simulated comparison with t-test

library(doSNOW)
library(itertools)

# Symmetric sample sizes ----------------------------------

# non-symmetric sample sizes ------------------------------
nx <- c(50, 200, 350)
nLen <- length(nx)
ny <- 500

mux <- c(0.75, 1)
muLen <- length(mux)

M <- 100 # M is number of repetitions for each scenario

nClust <- 6

params <- data.frame(nx=rep(nx, each=muLen*M),
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

  sub$timePred <- NA
  sub$timeExpert6 <- NA

  sub$mStop <- NA
  sub$EmStop <- NA

  sub$maxErrorExpert6 <- NA

  for (i in 1:nrow(sub)){

    nx <- sub$nx[i]
    x <- rnorm(n = nx, mean = sub$mux[i], sd=1)
    y <- rnorm(n = ny, mean = 0, sd=1)

    sub$pt[i] <- t.test(x, y, var.equal = TRUE)$p.value

    sub$timePred[i] <- system.time(fp <- 
      fastPerm(x, y, testStat = diffMean))[3]
    sub$pPred[i] <- fp$pPred
    sub$mStop[i] <- fp$mStop

    sub$EmStop[i] <- mStopDiffMean(x,y)

    fpAsym <- fastPermAsym(x,y, testStat=diffMean)
    sub$pAsymNorm[i] <- fpAsym$pNorm
    sub$pAsymT[i] <- fpAsym$pT
    
    # using EXPERT package
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

system.time(nonSymResults <- foreach(j=blocks, .combine=rbind) %dopar% {
                   return(simFun(params[j,]))
                 })

stopCluster(cl)

save(nonSymResults, file ="nonSymResultsDiff_parallel.RData")
