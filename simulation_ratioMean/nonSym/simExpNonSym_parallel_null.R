# Simulations comparing proposed algorithm to SAMC
# exponential data, unequal sample size

# arg <- commandArgs(TRUE)

library(doSNOW)
library(itertools)

# Nonsymmetric sample sizes ----------------------------------
nx <- c(20, 40, 60)
nLen <- length(nx)
ny <- 100

ratey <- 1
rateLen <- length(ratey)

M <- 100 # M is number of repetitions for each scenario

nClust <- 6

params <- data.frame(nx=rep(nx, each=rateLen*M),
    ratey=rep(ratey, each=M)
)

simFun <- function(sub){

  library(fastPerm)
  library(EXPERT)
  library(PearsonDS)

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
  sub$pPerm <- NA

  sub$timePred <- NA
  sub$timeExpert6 <- NA

  sub$mStop <- NA
  sub$EmStop <- NA

  sub$maxErrorExpert6 <- NA

  nIter <- 1e5

  for (i in 1:nrow(sub)){
    print(i)

    nx <- sub$nx[i]
    x <- rexp(n = nx, rate = 1)
    y <- rexp(n = ny, rate = sub$ratey[i])

    # beta prime distribution
    t0 <- max(mean(x) / mean(y), mean(y) / mean(x))
    sub$pBeta[i] <- ppearsonVI(t0, a = nx, b = ny, scale = ny/nx, location = 0,
                               lower.tail = FALSE)+
                    ppearsonVI(t0, a = ny, b = nx, scale = nx/ny, location = 0,
                               lower.tail = FALSE)

    # MC simulation
    z <- c(x, y)
    N <- length(z)
    tStar <- rep(NA, length = N)
    t0 <- mean(x) / mean(y)
    for (j in 1:nIter){
      piStar <- sample(1:N)
      xStar <- z[piStar[1:nx]]
      yStar <- z[piStar[(nx+1):N]]
      tStar[j] <- mean(xStar) / mean(yStar)
    }
    tStarMax <- pmax(tStar, 1/tStar)
    t0Max <- max(t0, 1/t0)
    sub$pPerm[i] <- mean(c(tStarMax, t0Max) >= t0Max)

    try({
      sub$timePred[i] <- system.time(fp <- 
        fastPerm(x, y, testStat = ratioMean, adjusted = TRUE)
        )[3]
      sub$pPred[i] <- fp$pPred
      sub$mStop[i] <- fp$mStop
      sub$EmStop[i] <- mStopRatioMean(x,y)
      })

    try({
      fpAsym <- fastPermAsym(x,y, testStat=ratioMean)
      sub$pAsymNorm[i] <- fpAsym$pNorm
      sub$pAsymT[i] <- fpAsym$pT
    })

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

# sub$t0GT1 <- sub$t0 >=1
# plot(sub$pPerm, sub$pPred, col = "blue")
# plot(sub$pPerm, sub$pBeta, col = "green")
# ggplot(aes(x = pPerm, y = pBeta, color = t0GT1), data= sub)+
#   geom_point()
# points(sub$pPerm, sub$pAsymNorm)
# # points(sub$pPerm, 2*sub$pAsymNorm)
# abline(a = 0, b = 1, col = "red")
# legend("topleft", legend = c("Alg", "pBeta", "pAsym"),
#   col = c("blue", "green", "black"), lty = 1)

cl <- makeCluster(nClust, type="SOCK")
registerDoSNOW(cl)

blocks <- isplitIndices(nrow(params), chunks=nClust)

system.time(nonSymResults <- foreach(j=blocks, .combine=rbind) %dopar% {
                   return(simFun(params[j,]))
                 }
           )

stopCluster(cl)

save(nonSymResults, file ="nonSymResultsExp_parallel_null.RData")

