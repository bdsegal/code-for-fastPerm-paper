library(ggplot2)
library(fastPerm)

### Analysis after computing p-values
load("cancerAnalysis.RData")

# get top genes
resultSort <- result[order(result$pPred),]

# Prep for latex table, use excel2LaTeX macro to make table
geneNames <- gsub("\\|[0-9]*", "", as.character(resultSort[, "geneNames"]))

# replicate test for top 15 to get variability ------------
# load("cancerAnalysisImage.RData")
indexVar <- resultSort[1:15, "index"]
geneNameVar <- geneNames[1:15]

  N <- 1000
  
  pPred <- matrix(nrow = N, ncol = length(indexVar))
  mStop <- matrix(nrow = N, ncol = length(indexVar))
  deviance <- matrix(nrow = N, ncol = length(indexVar))
  aic <- matrix(nrow = N, ncol = length(indexVar))
  
  colnames(pPred) <- geneNameVar
  colnames(mStop) <- geneNameVar
  colnames(deviance) <- geneNameVar
  colnames(aic) <- geneNameVar

  time  <- rep(NA, N*length(indexVar))

  count <- 1
  for (n in 1:N){
    print(paste(n, " of ", N, sep = ""))

    for (i in 1:length(indexVar)){
    
      x <- LUADmatAbove[indexVar[i],]
      y <- LUSCmatAbove[indexVar[i],]
      
      time[count] <- system.time(fp <- fastPerm(x, y, testStat=ratioMean, B=1000, adjusted=FALSE))[3]
      pPred[n,i] <- fp$pPred
      mStop[n,i] <- fp$m
      deviance[n,i] <- fp$deviance
      aic[n,i] <- fp$aic

      count <- count + 1
    }
  }

# average time to calculate pvals for top 15
summary(time)

pPredQuants <- signif(t(apply(pPred, 2, quantile, probs = c(0.1, 0.5, 0.9))), 3)

signif(log(pPredQuants,10),3)

save(pPred, time, mStop, deviance, aic, file="replicateTop15.RData")
