# Plot results

library(reshape2)
library(ggplot2)
library(dplyr)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symmetric sample size
load("sym/symResultsDiff_parallel_uneqVar_null.RData")

nVec <- unique(symResults$n)
alpha <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  sub <- symResults[symResults$n == nVec[i],
                   c("pPredStud", "pPred", "pPermStud", "pPerm", "pt", "pAsymNorm")]
  for (j in 1:length(alpha)) {
    falseAlarmTab <- 
     rbind(falseAlarmTab,
           c(nVec[i], alpha[j], apply(sub, 2, function(x) {mean(x <= alpha[j])})))
  }
}
colnames(falseAlarmTab)[1:2] <- c("n", "alpha")
falseAlarmTab <- as.data.frame(falseAlarmTab)
falseAlarmTab <- falseAlarmTab[with(falseAlarmTab, order(alpha, n)), ]
falseAlarmTab[, c(2,1,3:ncol(falseAlarmTab))]
#   alpha  n pPredStud pPred pPermStud pPerm    pt pAsymNorm
# 1  0.01 20     0.018 0.018     0.013 0.013 0.010     0.012
# 4  0.01 40     0.013 0.013     0.010 0.010 0.009     0.009
# 7  0.01 60     0.014 0.014     0.011 0.011 0.009     0.011
# 2  0.05 20     0.078 0.078     0.056 0.056 0.051     0.057
# 5  0.05 40     0.071 0.070     0.056 0.056 0.054     0.055
# 8  0.05 60     0.070 0.070     0.055 0.055 0.051     0.056
# 3  0.10 20     0.137 0.136     0.107 0.107 0.102     0.113
# 6  0.10 40     0.130 0.131     0.096 0.096 0.095     0.096
# 9  0.10 60     0.127 0.130     0.106 0.106 0.104     0.106

plotData <- melt(symResults, id.vars=c("n", "pPermStud"), 
  measure.vars=c("pPredStud", "pt", "pPred", "pPerm"))

levels(plotData$variable) <- c("Alg 1 student", "t-test", "Alg 1", "MC")

dev.new(width=10, height=5)
ggplot(aes(x=pPermStud, y=value, color=as.factor(n),
           shape=as.factor(n)), data=plotData)+
  geom_point(size=1.5)+
  theme_bw(21)+
  facet_grid( ~ variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x="MC student", y="p")+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simDiff_sym_uneqVar_null.png"))


# plotData <- melt(symResults, id.vars=c("n", "pt"), 
#   measure.vars=c("pPredStud", "pPred", "pPermStud", "pPerm"))

# plotData$student <- "Not studentized"
# plotData$student[grep("Stud", plotData$variable)] <- "Studentized"
# plotData$student <- factor(plotData$student, levels = c("Studentized", "Not studentized"))
# plotData$method <- "Alg 1"
# plotData$method[grep("Perm", plotData$variable)] <- "Monte Carlo"

# # levels(plotData$variable) <- c("Alg 1","Asym","SAMC")

# dev.new(width=9, height=6.5)
# ggplot(aes(x=pt, y=value, color=as.factor(n),
#     shape=as.factor(n)),
#   data=plotData)+
#   geom_point(size=1.5)+
#   theme_bw(24)+
#   facet_grid(student~method)+
#   geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
#   labs(x=expression(p[t]), y="p")+
#   scale_color_discrete("n")+
#   scale_shape_discrete("n")+
#   # scale_x_continuous(breaks=seq(0,-120,-30))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
# ggsave(file.path(paperPath,"simDiff_sym_uneqVar_null.png"))


# non-symmetric sample size
load(file="nonSym/nonSymResultsDiff_parallel_uneqVar_null.RData")

nVec <- unique(nonSymResults$nx)
alpha <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  sub <- nonSymResults[nonSymResults$nx == nVec[i],
                   c("pPredStud", "pPred", "pPermStud", "pPerm", "pt", "pAsymNorm")]
  for (j in 1:length(alpha)) {
    falseAlarmTab <- 
     rbind(falseAlarmTab,
           c(nVec[i], alpha[j], apply(sub, 2, function(x) {mean(x <= alpha[j])})))
  }
}
colnames(falseAlarmTab)[1:2] <- c("nx", "alpha")
falseAlarmTab <- as.data.frame(falseAlarmTab)
falseAlarmTab <- falseAlarmTab[with(falseAlarmTab, order(alpha, nx)), ]
falseAlarmTab[, c(2,1,3:ncol(falseAlarmTab))]
#   alpha nx pPredStud pPred pPermStud pPerm    pt pAsymNorm
# 1  0.01 20     0.024 0.129     0.011 0.175 0.008     0.167
# 4  0.01 40     0.017 0.062     0.016 0.078 0.012     0.075
# 7  0.01 60     0.008 0.027     0.008 0.029 0.008     0.027
# 2  0.05 20     0.090 0.258     0.058 0.284 0.053     0.279
# 5  0.05 40     0.071 0.165     0.060 0.170 0.059     0.164
# 8  0.05 60     0.049 0.109     0.040 0.100 0.035     0.095
# 3  0.10 20     0.152 0.347     0.125 0.370 0.124     0.366
# 6  0.10 40     0.120 0.247     0.102 0.242 0.101     0.235
# 9  0.10 60     0.113 0.193     0.086 0.170 0.083     0.169

plotData <- melt(nonSymResults, id.vars=c("nx", "pPermStud"), 
  measure.vars=c("pPredStud", "pt", "pPred", "pPerm"))

levels(plotData$variable) <- c("Alg 1 student", "t-test", "Alg 1", "MC")

dev.new(width=10, height=5)
ggplot(aes(x=pPermStud, y=value, color=as.factor(nx),
  shape=as.factor(nx)
  ), data=plotData)+
  geom_point(size=1.5)+
  theme_bw(21)+
  facet_grid( ~ variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x="MC student", y="p")+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simDiff_nonSym_uneqVar_null.png"))