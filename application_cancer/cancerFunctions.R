# Functions for the cancer analysis in Section 7-----------

importFun <- function(path){
  # Function for importing data
  
  files <- list.files(path)
  
  tables <- list()
  for (i in 1:length(files)){
    print(paste(i, " of ", length(files), sep=""))
    tables[[i]] <- read.table(paste(path, files[i], sep = ""),
            header = TRUE, stringsAsFactors = FALSE)
  }

  uniqueID <- list()
  for (i in 1:length(tables)){
    uniqueID[[i]] <- unique(tables[[i]]$gene_id)
  }

  allGeneID <- unique(do.call(c, uniqueID))

  # rows=genes, columns=subjects
  mat <- matrix(nrow = length(allGeneID),
            ncol = length(tables))
  rownames(mat) <- allGeneID

  # if 0 then ok to combine the tables
  namesOff <- sum(sapply(tables, function(x){sum(x$gene_id != allGeneID)}))

  # combine the tables
  for (i in 1:length(tables)){
    mat[,i] <- tables[[i]]$normalized_count
  }
  
  return(list(namesOff=namesOff, mat=mat))
}

ratioMeanFun <- function(mat1, mat2){
  # Utility function for ratioPerm
  
  mat1_mean<- apply(mat1, 1, mean)
  mat2_mean <- apply(mat2, 1, mean)  

  return(mat1_mean / mat2_mean)
}

ratioMedianFun <- function(mat1, mat2){
  # Utility function for ratioPerm

  mat1_median <- apply(mat1, 1, median)
  mat2_median <- apply(mat2, 1, median)  

  return(mat1_median / mat2_median)
}

ratioPerm <- function(mat1, mat2, B = 1000){
  # Function for running initial permutation test

  # meanRatio0 <- ratioMeanFun(mat1, mat2)
  
  all_mat <- cbind(mat1, mat2)  
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  n <- n1 + n2
  
  meanRatio <- matrix(nrow = B, ncol = nrow(mat1))
  
  for (b in 1:B){
    print(paste(b, " of ", B, sep=""))

    star <- sample(1:n)
    
    mat1_new <- all_mat[, star[1:n1]]
    mat2_new <- all_mat[, star[(n1 + 1):n]]
    
    meanRatio[b,] <- ratioMeanFun(mat1_new, mat2_new)

  }
  
  # return(list(meanRatio = meanRatio, meanRatio0 = meanRatio0))
  return(meanRatio)
}

pFunRatio <- function(x, t0){
  # Function for calculating p-values from initial permutation test
  
  if(!is.nan(t0)){
  
  B <- length(x)

    # Note: the following two-sided calculation with t=x/y is equivalent
    # to using a one sided calculation with t=max(x/y, y/x)
    if (t0 > 1){
      (sum(x >= t0) + sum(x <= 1/t0)) / B
    } else if (t0 < 1) {
      (sum(x <= t0) +  sum(x >= 1/t0)) / B
    } else {
      return(1)
    }
  } else { NA }
}
  
pRatioBeta <- function(x,y){
  
  nx <- length(x)
  ny <- length(y)

  t <- mean(x)/mean(y)
  r1 <- nx/ny*t
  r2 <- ny/nx*t
  return(pbeta(r1/(1+r1), nx, ny, lower.tail=FALSE)+
    pbeta(r2/(1+r2), ny, nx, lower.tail=FALSE)
    )
}

parFun <- function(genesToTestSub){
# Function for computing pPred in parallel
  library(fastPerm)

  sub <- data.frame(geneNames = rownames(LUADmatAbove)[genesToTestSub],
        index = genesToTestSub,
        pPred = NA,
        mStop = NA,
        deviance = NA,
        aic = NA,
        df.residual = NA,
        EmStop = NA,
        pAsymNorm = NA, 
        pAsymT = NA,
        pBeta = NA)

  for (i in 1:length(genesToTestSub)){
    x <- LUADmatAbove[genesToTestSub[i],]
    y <- LUSCmatAbove[genesToTestSub[i],]
    
    try({
      fp <- fastPerm(x, y, testStat = ratioMean)

      sub$pPred[i] <- fp$pPred
      sub$deviance[i] <- fp$deviance
      sub$aic[i] <- fp$aic
      sub$df.residual[i] <- fp$df.residual
      sub$mStop[i] <- fp$mStop
      sub$EmStop[i] <- mStopRatioMean(x,y)

      fpAsym <- fastPermAsym(x,y, testStat=ratioMean)
      sub$pAsymNorm[i] <- fpAsym$pNorm
      sub$pAsymT[i] <- fpAsym$pT

      sub$pBeta[i] <- pRatioBeta(x,y)
    })
  }

  sub$t0LUADoverLUSC <- ratioMean0[genesToTestSub]

  return(sub)
}
