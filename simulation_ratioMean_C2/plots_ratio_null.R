# summarize results
library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

load("sym/symResultsExp_parallel_null_revision.RData")

nVec <- unique(symResults$n)
alpha <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  sub <- symResults[symResults$n == nVec[i],
                   c("pPermAdj", "pPerm", "pPred", "pAsymNorm", "pDeltaUneq", "pBeta")]
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
#   alpha  n pPermAdj pPerm pPred pAsymNorm pDeltaUneq pBeta
# 1  0.01 20    0.010 0.010 0.016     0.066      0.003 0.009
# 4  0.01 40    0.010 0.010 0.018     0.050      0.002 0.008
# 7  0.01 60    0.013 0.013 0.013     0.031      0.006 0.015
# 2  0.05 20    0.064 0.064 0.084     0.143      0.045 0.058
# 5  0.05 40    0.061 0.061 0.079     0.112      0.054 0.061
# 8  0.05 60    0.051 0.051 0.063     0.091      0.050 0.047
# 3  0.10 20    0.109 0.109 0.151     0.207      0.116 0.107
# 6  0.10 40    0.106 0.106 0.137     0.172      0.110 0.111
# 9  0.10 60    0.093 0.093 0.111     0.139      0.095 0.092

falseAlarmTabM <- melt(falseAlarmTab, id.vars = c("n", "alpha"))
ggplot(aes(y = value, x = variable), data = falseAlarmTabM)+
  geom_bar(stat = "identity")+
  facet_grid(alpha ~ n)+
  theme_bw(15)+
  geom_hline(aes(yintercept = alpha), linetype = "dashed",
             color = "red", size = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))+
  labs(x = "", y = "Type I error rate")


plotData <- melt(symResults, id.vars=c("n", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "pDeltaUneq", "pBeta"))

levels(plotData$variable) <- c("Alg 1","Asym", "Delta", "Beta prime")

dev.new(width=9, height=5)
ggplot(aes(x=pPerm, y=value, color=as.factor(n),
    shape=as.factor(n)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)), y= "p")+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  # scale_x_continuous(breaks = seq(0, 1, 0.25), 
  #                    labels = c("0", "0.25", "0.5", "0.75", "1"))
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simExp_sym_null.png"))

dev.new(width=6, height=5)
ggplot(aes(x=pBeta, y=log(mStop*1000,10), color=as.factor(n),
  shape=as.factor(n)),
  data=symResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(p[beta]))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simExp_sym_iter_null.png"))

summary(symResults$pExpert6)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 # 0.0003  0.0005  0.0009  0.0012  0.0015  0.0031    2989 

# nonsymmetric sample sizes
load("nonSym/nonSymResultsExp_parallel_null_revision.RData")


nVec <- unique(nonSymResults$nx)
alpha <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  sub <- nonSymResults[nonSymResults$nx == nVec[i],
                   c("pPermAdj", "pPerm", "pPred", "pAsymNorm", "pDeltaUneq", "pBeta")]
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
#   alpha nx pPermAdj pPerm pPred pAsymNorm pDeltaUneq pBeta
# 1  0.01 20    0.011 0.011 0.016     0.054      0.008 0.012
# 4  0.01 40    0.008 0.008 0.012     0.033      0.004 0.006
# 7  0.01 60    0.012 0.012 0.016     0.035      0.007 0.014
# 2  0.05 20    0.061 0.061 0.082     0.127      0.065 0.056
# 5  0.05 40    0.048 0.048 0.062     0.097      0.047 0.050
# 8  0.05 60    0.047 0.047 0.065     0.083      0.044 0.051
# 3  0.10 20    0.118 0.118 0.161     0.190      0.138 0.116
# 6  0.10 40    0.102 0.102 0.141     0.171      0.107 0.104
# 9  0.10 60    0.091 0.091 0.118     0.136      0.093 0.088

falseAlarmTabM <- melt(falseAlarmTab, id.vars = c("nx", "alpha"))
ggplot(aes(y = value, x = variable), data = falseAlarmTabM)+
  geom_bar(stat = "identity")+
  facet_grid(alpha ~ nx)+
  theme_bw(15)+
  geom_hline(aes(yintercept = alpha), linetype = "dashed",
             color = "red", size = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))+
  labs(x = "", y = "Type I error rate")


plotData <- melt(nonSymResults, id.vars=c("nx", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "pDeltaUneq", "pBeta"))

levels(plotData$variable) <- c("Alg 1","Asym", "Delta", "Beta prime")

dev.new(width=9, height=5)
ggplot(aes(x=pPerm, y=value, color=as.factor(nx),
    shape=as.factor(nx)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)),
    y="p")+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simExp_nonSym_null.png"))

dev.new(width=6, height=5)
ggplot(aes(x=pBeta, y=log(mStop*1000,10), color=as.factor(nx),
  shape=as.factor(nx)),
  data=nonSymResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(p[beta]))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simExp_nonSym_iter_null.png"))

summary(nonSymResults$pExpert6)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 # 0.0001  0.0005  0.0009  0.0016  0.0023  0.0058    2986 