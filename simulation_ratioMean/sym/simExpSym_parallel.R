# Simulations comparing proposed algorithm to SAMC
# exponential data, equal sample size

# arg <- commandArgs(TRUE)

library(doSNOW)
library(itertools)

# Symmetric sample sizes ----------------------------------
n <- c(100, 500, 1000)
nLen <- length(n)

ratey <- c(1.75, 2.25)
rateLen <- length(ratey)

M <- 100 # M is number of repetitions for each scenario

nClust <- 6

params <- data.frame(n=rep(n, each=rateLen*M),
    ratey=rep(ratey, each=M)
)

simFun <- function(sub){

  library(fastPerm)
  library(EXPERT)

  ratio.test.statistic <- function(data.input){
      xbar <- mean(data.input$x)
      ybar <- mean(data.input$y)

      max(xbar/ybar, ybar/xbar)
  }

  sub$pBeta <- NA
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
    x <- rexp(n = n, rate = 1)
    y <- rexp(n = n, rate = sub$ratey[i])

  # R ~ B'(nx, ny) and 1/R ~ B'(ny, nx)
    R <- mean(x)/mean(y)
    sub$pBeta[i] <- pbeta(R/(1+R), n, n, lower.tail=FALSE)+
      pbeta(1/(1+R), n, n, lower.tail=TRUE)

    sub$timePred[i] <- system.time(fp <- 
      fastPerm(x, y, testStat = ratioMean))[3]
    sub$pPred[i] <- fp$pPred
    sub$mStop[i] <- fp$mStop

    sub$EmStop[i] <- mStopRatioMean(x,y)

    fpAsym <- fastPermAsym(x,y, testStat=ratioMean)
    sub$pAsymNorm[i] <- fpAsym$pNorm
    sub$pAsymT[i] <- fpAsym$pT
    
    # using EXPERT package
    data.input<-list(x=x, y=y)
    t.obs<-ratio.test.statistic(data.input)

    try({
      sub$timeExpert6[i] <- 
        system.time(res6 <- SAMC.adapt(data.input, t.obs, t.start=0,
          n.iter.1=5e4, n.iter.2=1e6, 
          n.region.1=101, n.region.2=301,
          prop.change=0.05, gain.factor.t0=1000,
          fun.test.statistic=ratio.test.statistic,
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

symResults <- foreach(j=blocks, .combine=rbind) %dopar% {
                   return(simFun(params[j,]))
                 }

stopCluster(cl)

save(symResults, file ="symResultsExp_parallel.RData")
