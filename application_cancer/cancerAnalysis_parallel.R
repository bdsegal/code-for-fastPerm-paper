
library(doSNOW)
library(itertools)

source("cancerFunctions.R")

nClust <- 6

# Import data ---------------------------------------------

LUADpath <- "data/unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.Level_3.1.12.0/"
LUSCpath <- "data/unc.edu_LUSC.IlluminaHiSeq_RNASeqV2.Level_3.1.8.0/"

LUADfiles <- list.files(LUADpath)
LUSCfiles <- list.files(LUSCpath)

# Get data
imptLUAD <- importFun(LUADpath)
imptLUAD$namesOff # namesOff == 0 if imported correctly
LUADmat <- imptLUAD$mat

imptLUSC <- importFun(LUSCpath)
imptLUSC$namesOff # namesOff == 0 if imported correctly
LUSCmat <- imptLUSC$mat

# check that genes are the same in the 2 samples
# (ok if setdiff is empty)
setdiff(rownames(LUADmat), rownames(LUSCmat))

# check that genes are in the same order in the 2 samples
# (ok if sum is 0)
sum(rownames(LUADmat) != rownames(LUSCmat))  

# if the two checks above are ok, cbind samples together
allMat <- cbind(LUSCmat, LUADmat)

# Minimum cutoff for testing ------------------------------

expressionCutoff <- quantile(allMat, prob = 0.25)
expressionCutoff

percAboveCut <- apply(allMat, 1, function(x){mean(x >= expressionCutoff)})

# number of genes with at least 50% of values above cutoff
sum(percAboveCut >= 0.5)

# total number of genes in the sample
nrow(allMat)

# limit remaining analysis to genese with at least 50%
# of measurements above cutoff
genesAbove <- which(percAboveCut >= 0.5)
LUADmatAbove <- LUADmat[genesAbove, ]
LUSCmatAbove <- LUSCmat[genesAbove, ]

# Initial permutation tests -------------------------------
# insure integer-valued number of permutations >=1000
Bsub <- ceiling(1000/nClust)

cl <- makeCluster(nClust, type="SOCK")
registerDoSNOW(cl)

# ratioMeanPerms is a matrix with B rows
# Each row contains test statistics for all genes from a single permutations
system.time(ratioMeanPerms <- foreach(j=1:nClust, .combine=rbind) %dopar% {
        return(ratioPerm(LUADmatAbove, LUSCmatAbove, B = Bsub))
       })

stopCluster(cl)

# vector of observed test statistics, one element per gene
ratioMean0 <- ratioMeanFun(LUADmatAbove, LUSCmatAbove)

# get all p-values
pMean <- rep(NA,ncol(ratioMeanPerms))
for (i in 1:ncol(ratioMeanPerms)){
  pMean[i] <- pFunRatio(ratioMeanPerms[,i], ratioMean0[i])
}

# Number of genes with p-value too small to calculate
sum(pMean == 0, na.rm=TRUE)

# index of genes to test with the pPred algorithm
genesToTest <- which(pMean == 0)
genesNotToTest <- which(pMean > 0)

# Results from pPred algorithm, computed in parallel
cl <- makeCluster(nClust, type="SOCK")
registerDoSNOW(cl)

blocks <- isplitIndices(length(genesToTest), chunks=nClust)

system.time(result <- foreach(j=blocks, .combine=rbind) %dopar% {
        return(parFun(genesToTest[j]))
       })

stopCluster(cl)

# save everything just in case
# save.image(file="cancerAnalysisImage.RData")

# save quantities for subsequent analyses
save(result,
     genesAbove,
     genesToTest, genesNotToTest,
     ratioMeanPerms, ratioMean0, pMean,
     LUADmatAbove, LUSCmatAbove,
     file="cancerAnalysis.RData")