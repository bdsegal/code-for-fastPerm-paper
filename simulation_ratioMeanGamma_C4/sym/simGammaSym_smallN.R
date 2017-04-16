# Simulations comparing proposed algorithm to SAMC
# exponential data, equal sample size

# arg <- commandArgs(TRUE)

library(doSNOW)
library(itertools)

nClust <- 6

# Symmetric sample sizes ----------------------------------
n <- c(20, 40, 60)
lambdax <- 1

lambday <- c(7, 12.5)
params05 <- expand.grid(n = n, lambday = lambday, lambdax = lambdax, alpha = 0.5)
params05

lambday <- c(5, 10)
params3 <- expand.grid(n = n, lambday = lambday, lambdax = lambdax, alpha = 3)

lambday <- c(5, 7.5)
params5 <- expand.grid(n = n, lambday = lambday, lambdax = lambdax, alpha = 5)

M <- 100 # M is number of repetitions for each scenario

params <- rbind(params05[rep(1:nrow(params05), each = M), ],
                params3[rep(1:nrow(params3), each = M), ],
                params5[rep(1:nrow(params5), each = M), ])

simFun <- function(sub){

  library(fastPerm)
  library(EXPERT)
  library(PearsonDS)

  pDeltaFun <- function(x, y, varEqual = FALSE) {

    nx <- length(x)
    ny <- length(y)
    xbar <- mean(x)
    ybar <- mean(y)

    if (varEqual) {
      sdx <- sdy <- sd(c(x, y))
    } else {
      sdx <- sd(x)
      sdy <- sd(y)
    }

    r <- xbar / ybar
    vx <- sqrt((sdx^2 / nx) / ybar^2 + (sdy^2 / ny) * (xbar^2 / ybar^4))
    vy <- sqrt((sdy^2 / ny) / xbar^2 + (sdx^2 / nx) * (ybar^2 / xbar^4))
  
    if(r >= 1) {
      pDelta <- pnorm(r, mean = 1, sd = vx, lower.tail = FALSE) +
                   pnorm(1 / r, mean = 1, sd = vy, lower.tail = TRUE)
    } else {
      pDelta <- pnorm(r, mean = 1, sd = vx, lower.tail = TRUE) +
                   pnorm(1 / r, mean = 1, sd = vy, lower.tail = FALSE)
    }
    return(pDelta)
  }

  ratio.test.statistic <- function(data.input){
      xbar <- mean(data.input$x)
      ybar <- mean(data.input$y)

      max(xbar/ybar, ybar/xbar)
  }

  sub$pDeltaUneq <- NA
  sub$pDeltaEq <- NA
  sub$pBeta <- NA
  sub$pPred <- NA
  sub$pAsymNorm <- NA
  sub$pAsymT <- NA
  sub$pExpert6 <- NA

  # sub$timePred <- NA
  # sub$timeExpert6 <- NA

  sub$mStop <- NA

  sub$EmStop <- NA

  sub$maxErrorExpert6 <- NA

  for (i in 1:nrow(sub)){

    n <- sub$n[i]
    x <- rgamma(n = n, shape = sub$alpha[i], rate = sub$lambdax[i])
    y <- rgamma(n = n, shape = sub$alpha[i], rate = sub$lambday[i])

    # delta method with plug-in estimates 
    sub$pDeltaUneq[i] <- pDeltaFun(x, y, varEqual = FALSE)
    sub$pDeltaEq[i] <- pDeltaFun(x, y, varEqual = TRUE)

    # beta prime distribution
    t0 <- max(mean(x) / mean(y), mean(y) / mean(x))
    sub$pBeta[i] <- ppearsonVI(t0, a = n * sub$alpha[i] , b = n * sub$alpha[i],
                               scale = 1, location = 0, lower.tail = FALSE)+
                    ppearsonVI(t0, a = n*sub$alpha[i] , b = n*sub$alpha[i],
                               scale = 1, location = 0, lower.tail = FALSE)

    try({
      sub$timePred[i] <- system.time(fp <- 
        fastPerm(x, y, testStat = ratioMean))[3]
      sub$pPred[i] <- fp$pPred
      sub$mStop[i] <- fp$mStop
    })

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

save(symResults, file = "symResultsGamma_smallN.RData")
