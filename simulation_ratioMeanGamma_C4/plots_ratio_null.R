# summarize results
library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

load("sym/symResultsGamma_null.RData")

nVec <- unique(symResults$n)
alphaVec <- unique(symResults$alpha)
signifLevel <- c(0.01, 0.05, 0.1)

head(symResults)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  for (k in 1:length(alphaVec)) {
    sub <- symResults[which(symResults$n == nVec[i] & symResults$alpha == alphaVec[k]),
                     c("pPerm", "pPred", "pAsymNorm", "pDeltaUneq", "pBeta")]
    for (j in 1:length(signifLevel)) {
      falseAlarmTab <- 
       rbind(falseAlarmTab,
        c(alphaVec[k], signifLevel[j], nVec[i], apply(sub, 2, function(x) {mean(x <= signifLevel[j])})))
    }
  }
}
colnames(falseAlarmTab)[1:3]<- c("alpha", "signifLevel", "n")
falseAlarmTab <- as.data.frame(falseAlarmTab)
falseAlarmTab <- falseAlarmTab[with(falseAlarmTab, order(alpha, signifLevel, n)), ]
signif(falseAlarmTab, 3)
#    alpha signifLevel  n pPerm pPred pAsymNorm pDeltaUneq pBeta
# 1    0.5        0.01 20 0.013 0.018     0.093      0.002 0.015
# 2    0.5        0.05 20 0.050 0.076     0.182      0.026 0.053
# 3    0.5        0.10 20 0.110 0.136     0.243      0.106 0.108
# 4    3.0        0.01 20 0.007 0.012     0.027      0.003 0.006
# 5    3.0        0.05 20 0.043 0.067     0.088      0.046 0.044
# 6    3.0        0.10 20 0.095 0.126     0.143      0.103 0.090
# 7    5.0        0.01 20 0.009 0.015     0.023      0.009 0.009
# 8    5.0        0.05 20 0.046 0.063     0.082      0.054 0.052
# 9    5.0        0.10 20 0.093 0.130     0.139      0.106 0.099
# 10   0.5        0.01 40 0.007 0.014     0.055      0.001 0.007
# 11   0.5        0.05 40 0.050 0.072     0.135      0.037 0.055
# 12   0.5        0.10 40 0.106 0.135     0.196      0.114 0.104
# 13   3.0        0.01 40 0.012 0.016     0.025      0.010 0.010
# 14   3.0        0.05 40 0.053 0.062     0.073      0.052 0.051
# 15   3.0        0.10 40 0.098 0.133     0.147      0.104 0.103
# 16   5.0        0.01 40 0.008 0.013     0.025      0.008 0.011
# 17   5.0        0.05 40 0.048 0.063     0.066      0.050 0.043
# 18   5.0        0.10 40 0.091 0.134     0.138      0.094 0.093
# 19   0.5        0.01 60 0.007 0.010     0.047      0.002 0.011
# 20   0.5        0.05 60 0.048 0.068     0.114      0.043 0.050
# 21   0.5        0.10 60 0.096 0.127     0.178      0.101 0.097
# 22   3.0        0.01 60 0.012 0.015     0.025      0.012 0.008
# 23   3.0        0.05 60 0.059 0.075     0.080      0.061 0.049
# 24   3.0        0.10 60 0.095 0.115     0.116      0.097 0.093
# 25   5.0        0.01 60 0.012 0.012     0.019      0.012 0.013
# 26   5.0        0.05 60 0.055 0.078     0.079      0.057 0.057
# 27   5.0        0.10 60 0.115 0.138     0.136      0.116 0.112

plotData <- melt(symResults, id.vars=c("n", "alpha", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "pDeltaUneq", "pBeta"))

levels(plotData$variable) <- c("Alg 1","Asym", "Delta", "Beta prime")

dev.new(width=9, height=6.5)
ggplot(aes(x=pPerm, y=value, color=as.factor(n),
    shape=as.factor(n)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(alpha ~ variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)), y= "p")+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  # scale_x_continuous(breaks = seq(0, 1, 0.25), 
  #                    labels = c("0", "0.25", "0.5", "0.75", "1"))
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simGamma_sym_null.png"))

# nonsymmetric sample sizes
load("nonSym/nonSymResultsGamma_null.RData")

nVec <- unique(nonSymResults$nx)
alphaVec <- unique(nonSymResults$alpha)
signifLevel <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  for (k in 1:length(alphaVec)) {
    sub <- nonSymResults[which(nonSymResults$nx == nVec[i] & nonSymResults$alpha == alphaVec[k]),
                     c("pPerm", "pPred", "pAsymNorm", "pDeltaUneq", "pBeta")]
    for (j in 1:length(signifLevel)) {
      falseAlarmTab <- 
       rbind(falseAlarmTab,
        c(alphaVec[k], signifLevel[j], nVec[i], apply(sub, 2, function(x) {mean(x <= signifLevel[j])})))
    }
  }
}
colnames(falseAlarmTab)[1:3]<- c("alpha", "signifLevel", "nx")
falseAlarmTab <- as.data.frame(falseAlarmTab)
falseAlarmTab <- falseAlarmTab[with(falseAlarmTab, order(alpha, signifLevel, nx)), ]
signif(falseAlarmTab, 3)
#    alpha signifLevel nx pPerm pPred pAsymNorm pDeltaUneq pBeta
# 1    0.5        0.01 20 0.011 0.015     0.065      0.006 0.011
# 10   0.5        0.01 40 0.015 0.018     0.053      0.003 0.013
# 19   0.5        0.01 60 0.008 0.011     0.042      0.003 0.012
# 2    0.5        0.05 20 0.043 0.069     0.128      0.047 0.053
# 11   0.5        0.05 40 0.057 0.072     0.133      0.048 0.056
# 20   0.5        0.05 60 0.052 0.071     0.112      0.045 0.050
# 3    0.5        0.10 20 0.098 0.121     0.179      0.109 0.091
# 12   0.5        0.10 40 0.113 0.141     0.195      0.119 0.108
# 21   0.5        0.10 60 0.106 0.126     0.172      0.109 0.098
# 4    3.0        0.01 20 0.011 0.016     0.023      0.012 0.011
# 13   3.0        0.01 40 0.005 0.011     0.027      0.005 0.009
# 22   3.0        0.01 60 0.011 0.013     0.017      0.011 0.011
# 5    3.0        0.05 20 0.047 0.070     0.073      0.059 0.039
# 14   3.0        0.05 40 0.058 0.065     0.069      0.057 0.054
# 23   3.0        0.05 60 0.053 0.066     0.070      0.050 0.052
# 6    3.0        0.10 20 0.088 0.128     0.135      0.104 0.087
# 15   3.0        0.10 40 0.094 0.124     0.124      0.101 0.089
# 24   3.0        0.10 60 0.094 0.119     0.117      0.097 0.097
# 7    5.0        0.01 20 0.010 0.014     0.022      0.007 0.009
# 16   5.0        0.01 40 0.011 0.011     0.017      0.011 0.009
# 25   5.0        0.01 60 0.015 0.020     0.025      0.015 0.018
# 8    5.0        0.05 20 0.058 0.074     0.085      0.066 0.054
# 17   5.0        0.05 40 0.046 0.057     0.059      0.048 0.052
# 26   5.0        0.05 60 0.059 0.081     0.085      0.061 0.062
# 9    5.0        0.10 20 0.110 0.145     0.143      0.121 0.114
# 18   5.0        0.10 40 0.081 0.114     0.108      0.085 0.088
# 27   5.0        0.10 60 0.113 0.145     0.138      0.118 0.115

plotData <- melt(nonSymResults, id.vars=c("nx", "alpha", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "pDeltaUneq", "pBeta"))

levels(plotData$variable) <- c("Alg 1","Asym", "Delta", "Beta prime")

dev.new(width=9, height=6.5)
ggplot(aes(x=pPerm, y=value, color=as.factor(nx),
    shape=as.factor(nx)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(alpha ~ variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)),
    y="p")+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simGamma_nonSym_null.png"))
