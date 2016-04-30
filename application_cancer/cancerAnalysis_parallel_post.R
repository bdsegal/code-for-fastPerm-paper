library(ggplot2)
paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

### Analysis after computing p-values
load("cancerAnalysis.RData")

# get number of p-values less than the Bonferroni cutoff
bon <- 0.05/length(genesAbove)
sigMean <- which(result$pPred <= bon)
length(sigMean)
#[1] 7846

# Distribution of pvals from algorithm
dev.new(width=5,height=5)
qplot(log(result[, "pPred"], 10), geom = "histogram")+
  labs(x = expression(paste("lo", g[10], "(", tilde(p)["pred"], ")", sep = "")),
    y = "count", title = "")+
  theme_bw(24)+
  geom_vline(xintercept = c(-3), color = "red", linetype = "dashed", size = 1.5)
ggsave(file.path(paperPath,"pPred_hist.png"))

# combined p-values
pAll <- c(result$pPred, pMean[genesNotToTest])

# Distribution of all pvals
dev.new(width=5,height=4.5)
qplot(pAll,geom="histogram")+
  theme_bw(24)+
  labs(x="p")
ggsave(file.path(paperPath,"p_all_cancer.png"))

# compare pPred to pAsym
plot(x=log(result$pAsymNorm,10), y=log(result$pPred,10))
abline(a=0,b=1, col="red")  

ggplot(aes(x=log(pAsymNorm,10), y=log(pPred,10)), data=result)+
  geom_point()+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+
  theme_bw(24)+
  labs(x=expression(paste("lo",g[10],"(",hat(p)["asym"],")")),
       y=expression(paste("lo",g[10],"(",tilde(p)["pred"],")")))

# compare pPred to pBeta
plot(x=log(result$pAsymNorm,10), y=log(result$pBeta,10))
abline(a=0,b=1, col="red")  

sum(result$pPred > 10^-3)
#[1] 351

# number of iterations
dev.new(width=5,height=5)
qplot(x=log(pPred,10), y=log(mStop*1000,10), data=result)+
theme_bw(24)+
labs(x=expression(paste("lo",g[10],"(",tilde(p)["pred"],")")),
     y=expression(paste("lo",g[10],"(total iterations)")))
ggsave(file.path(paperPath,"cancer_iter.png"))

# m_stop vs m_stop^asym
dev.new(width=5,height=5)
qplot(x=EmStop, y=mStop, data=result)+
  geom_abline(intercept=0, slope=1, color="red", size=1, linetype="dashed")+
  theme_bw(24)+
  labs(y=expression(m["stop"]),
       x=expression(m["stop"]^"asym"))
ggsave(file.path(paperPath,"mStopPlot.png"))


# get top genes
resultSort <- result[order(result$pPred),]

# Prep for latex table, use excel2LaTeX macro to make table
geneNames <- gsub("\\|[0-9]*", "", as.character(resultSort[, "geneNames"]))
pPredResults <- cbind(geneNames, logpPred = signif(log(resultSort$pPred,10), 3),
          signif(resultSort[, c(12,4,5,6)],3))
write.csv(pPredResults, file = "pPredResults.csv", row.names=FALSE)

pPredResults[1:15,]
#      geneNames    logpPred t0LUADoverLUSC mStop deviance   aic
# 2708      DSG3        -212         0.0100     5    40.10  68.1
# 4763      KRT5        -210         0.0107     4    12.50  38.2
# 2701      DSC3        -197         0.0175     6    41.50  72.1
# 1480    CALML3        -195         0.0138     6    57.80  90.0
# 9189      TP63        -193         0.0308     6    24.20  55.1
# 651     ATP1B3        -193         0.2250     5    28.60  57.7
# 7704     S1PR5        -190         0.0775     6    98.40 131.0
# 4765     KRT6B        -185         0.0173     5    45.40  76.1
# 9255    TRIM29        -183         0.0788     6    39.30  72.0
# 4475      JAG1        -180         0.1700     5    60.70  92.2
# 7155     PVRL1        -180         0.1100     6     8.33  39.2
# 1988     CLCA2        -178         0.0138     7    51.60  86.8
# 848       BNC1        -178         0.0244     7    76.80 112.0
# 3261      FAT2        -177         0.0339     7    53.50  89.0
# 4766     KRT6C        -177         0.0183     6    84.80 119.0