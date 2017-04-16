# Section 5, Figure 2
# Simulated comparison with t-test

library(doSNOW)
library(itertools)

# test stat with unequal variance for SAMC
t.test.stat.uneqVar <- function (data.input) 
{
    x <- data.input$x
    y <- data.input$y
    abs(t.test(x = x, y = y, var.equal = FALSE)$statistic)
}

# non-symmetric sample sizes ------------------------------
nx <- c(20, 40, 60)
nLen <- length(nx)
ny <- 100

mux <- 0
muLen <- length(mux)

M <- 1000 # M is number of repetitions for each scenario

nClust <- 6

params <- data.frame(nx=rep(nx, each=muLen*M),
    mux=rep(mux, each=M)
)

simFun <- function(sub){

  library(fastPerm)
  library(EXPERT)

  sub$pt <- NA
  sub$pPermStud <- NA
  sub$pPerm <- NA
  sub$pPredStud <- NA
  sub$pPred <- NA

  sub$pAsymNorm <- NA
  sub$pAsymT <- NA

  sub$timePredStud <- NA
  sub$timePred <- NA

  sub$EmStop <- NA
  sub$mStopStud <- NA
  sub$mStop <- NA

  sub$timeExpert6 <- NA
  sub$pExpert6 <- NA
  sub$maxErrorExpert6 <- NA

  sub$timeExpert6uneqVar <- NA
  sub$pExpert6uneqVar <- NA
  sub$maxErrorExpert6uneqVar <- NA

  for (i in 1:nrow(sub)){

    nx <- sub$nx[i]
    N <- nx + ny
    x <- rnorm(n = nx, mean = sub$mux[i], sd=3)
    y <- rnorm(n = ny, mean = 0, sd=1)

    sub$pt[i] <- t.test(x, y, var.equal = FALSE)$p.value

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

    # Monte Carlo
    z <- c(x, y)
    B <- 1e5
    Tstud <- rep(NA, B)
    T <- rep(NA, B)

    t0 <- mean(x) - mean(y)
    tStud0 <- t0 / sqrt(var(x) / nx + var(y) / ny)
    for (b in 1:B) {
      pi <- sample(1:N, replace = FALSE)
      xStar <- z[pi[1:nx]]
      yStar <- z[pi[(nx+1):N]]
      T[b] <- mean(xStar) - mean(yStar)
      Tstud[b] <- (mean(xStar) - mean(yStar)) / sqrt(var(xStar) / nx + var(yStar) / ny)
    }
    
    sub$pPermStud[i] <- mean(abs(Tstud) >= abs(tStud0))
    sub$pPerm[i] <- mean(abs(T) >= abs(t0))

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

    t.obs.uneqVar<-t.test.stat.uneqVar(data.input)
    try({
      sub$timeExpert6uneqVar[i] <- 
        system.time(res6uneqVar <- SAMC.adapt(data.input, t.obs.uneqVar, t.start=0,
          n.iter.1=5e4, n.iter.2=1e6, 
          n.region.1=101, n.region.2=301,
          prop.change=0.05, gain.factor.t0=1000,
          fun.test.statistic=t.test.stat.uneqVar,
          fun.proposal=proposal.permute.vector)
        )[3]
      # The estimated p-value
      sub$pExpert6uneqVar[i] <- res6uneqVar$p.value
      sub$maxErrorExpert6uneqVar[i] <- maxError(res6uneqVar)
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

save(nonSymResults, file ="nonSymResultsDiff_parallel_uneqVar_null.RData")
