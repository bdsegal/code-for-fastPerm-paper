# summarize results
library(ggplot2)
library(reshape2)
library(dplyr)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

load("sym/symResultsGammaDiff_null.RData")

head(symResults)

nVec <- unique(symResults$n)
alphaVec <- unique(symResults$alpha)
signifLevel <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  for (k in 1:length(alphaVec)) {
    sub <- symResults[which(symResults$n == nVec[i] & symResults$alpha == alphaVec[k]),
                     c("pPerm", "gammaDiffSaddle", "pPred", "pAsymNorm", "ptUneq")]
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
#    alpha signifLevel nx  pPerm gammaDiffSaddle  pPred pAsymNorm  ptUneq
# 1    0.5        0.01 20 0.0110          0.0100 0.0165    0.0060  0.0045
# 10   0.5        0.01 40 0.0125          0.0110 0.0150    0.0090  0.0085
# 19   0.5        0.01 60 0.0115          0.0085 0.0140    0.0105  0.0105
# 2    0.5        0.05 20 0.0495          0.0560 0.0665    0.0460  0.0410
# 11   0.5        0.05 40 0.0515          0.0490 0.0660    0.0520  0.0485
# 20   0.5        0.05 60 0.0455          0.0450 0.0595    0.0435  0.0425
# 3    0.5        0.10 20 0.1000          0.1020 0.1280    0.1020  0.0945
# 12   0.5        0.10 40 0.0995          0.0950 0.1260    0.1020  0.0975
# 21   0.5        0.10 60 0.0980          0.0950 0.1230    0.0990  0.0965
# 4    3.0        0.01 20 0.0115          0.0070 0.0165    0.0095  0.0095
# 13   3.0        0.01 40 0.0120          0.0115 0.0150    0.0120  0.0120
# 22   3.0        0.01 60 0.0075          0.0075 0.0080    0.0070  0.0070
# 5    3.0        0.05 20 0.0510          0.0465 0.0715    0.0515  0.0495
# 14   3.0        0.05 40 0.0545          0.0575 0.0680    0.0560  0.0525
# 23   3.0        0.05 60 0.0470          0.0475 0.0665    0.0480  0.0475
# 6    3.0        0.10 20 0.0940          0.0990 0.1280    0.0980  0.0940
# 15   3.0        0.10 40 0.0990          0.1000 0.1320    0.0990  0.0980
# 24   3.0        0.10 60 0.0980          0.0985 0.1230    0.0980  0.0980
# 7    5.0        0.01 20 0.0115          0.0095 0.0175    0.0115  0.0115
# 16   5.0        0.01 40 0.0090          0.0065 0.0130    0.0080  0.0080
# 25   5.0        0.01 60 0.0045          0.0055 0.0085    0.0040  0.0040
# 8    5.0        0.05 20 0.0525          0.0525 0.0675    0.0525  0.0505
# 17   5.0        0.05 40 0.0525          0.0545 0.0715    0.0535  0.0520
# 26   5.0        0.05 60 0.0460          0.0445 0.0580    0.0470  0.0470
# 9    5.0        0.10 20 0.0965          0.0960 0.1220    0.0980  0.0955
# 18   5.0        0.10 40 0.1070          0.1060 0.1370    0.1080  0.1080
# 27   5.0        0.10 60 0.0925          0.0905 0.1300    0.0940  0.0915

plotData <- melt(symResults, id.vars=c("n", "alpha", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "ptUneq", "gammaDiffSaddle"))
levels(plotData$variable)
levels(plotData$variable) <- c("Alg 1","Asym", "t-test", "Saddle")

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
ggsave(file.path(paperPath,"simGammaDiff_sym_null.png"))


# nonsymmetric sample sizes ---------------------------------------------------
load("nonSym/nonSymResultsGammaDiff_null.RData")

nVec <- unique(nonSymResults$nx)
alphaVec <- unique(nonSymResults$alpha)
signifLevel <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  for (k in 1:length(alphaVec)) {
    sub <- nonSymResults[which(nonSymResults$nx == nVec[i] & nonSymResults$alpha == alphaVec[k]),
                     c("pPerm", "gammaDiffSaddle", "pPred", "pAsymNorm", "ptUneq")]
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
#    alpha signifLevel nx  pPerm gammaDiffSaddle  pPred pAsymNorm  ptUneq
# 1    0.5        0.01 20 0.0095          0.0095 0.0105    0.0085  0.0245
# 10   0.5        0.01 40 0.0090          0.0060 0.0105    0.0070  0.0140
# 19   0.5        0.01 60 0.0130          0.0160 0.0170    0.0105  0.0135
# 2    0.5        0.05 20 0.0460          0.0465 0.0675    0.0440  0.0740
# 11   0.5        0.05 40 0.0455          0.0470 0.0620    0.0445  0.0540
# 20   0.5        0.05 60 0.0505          0.0500 0.0670    0.0495  0.0530
# 3    0.5        0.10 20 0.0915          0.0930 0.1260    0.0845  0.1220
# 12   0.5        0.10 40 0.0980          0.0945 0.1280    0.0960  0.1040
# 21   0.5        0.10 60 0.1100          0.1080 0.1410    0.1100  0.1080
# 4    3.0        0.01 20 0.0085          0.0095 0.0155    0.0085  0.0135
# 13   3.0        0.01 40 0.0135          0.0120 0.0185    0.0135  0.0140
# 22   3.0        0.01 60 0.0070          0.0055 0.0090    0.0070  0.0070
# 5    3.0        0.05 20 0.0440          0.0440 0.0665    0.0435  0.0480
# 14   3.0        0.05 40 0.0480          0.0555 0.0695    0.0485  0.0530
# 23   3.0        0.05 60 0.0470          0.0495 0.0635    0.0485  0.0460
# 6    3.0        0.10 20 0.0875          0.0885 0.1260    0.0885  0.1000
# 15   3.0        0.10 40 0.1050          0.1040 0.1350    0.1060  0.0975
# 24   3.0        0.10 60 0.1040          0.1080 0.1370    0.1040  0.1040
# 7    5.0        0.01 20 0.0140          0.0110 0.0200    0.0140  0.0145
# 16   5.0        0.01 40 0.0090          0.0100 0.0155    0.0090  0.0100
# 25   5.0        0.01 60 0.0105          0.0090 0.0120    0.0110  0.0075
# 8    5.0        0.05 20 0.0540          0.0535 0.0845    0.0540  0.0620
# 17   5.0        0.05 40 0.0530          0.0525 0.0730    0.0525  0.0555
# 26   5.0        0.05 60 0.0520          0.0510 0.0635    0.0520  0.0500
# 9    5.0        0.10 20 0.1140          0.1160 0.1520    0.1140  0.1130
# 18   5.0        0.10 40 0.0995          0.1000 0.1300    0.0995  0.1040
# 27   5.0        0.10 60 0.1040          0.0985 0.1320    0.1050  0.1060

plotData <- melt(nonSymResults, id.vars=c("nx", "alpha", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "ptUneq", "gammaDiffSaddle"))

levels(plotData$variable)
levels(plotData$variable) <- c("Alg 1","Asym", "t-test", "Saddle")

plotData$value[which(plotData$value > 1)] <- NA

dev.new(width=9, height=6.5)
ggplot(aes(x=pPerm, y=value, color=as.factor(nx),
  shape=as.factor(nx)),  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(alpha ~ variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)),
    y="p")+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simGammaDiff_nonSym_null.png"))