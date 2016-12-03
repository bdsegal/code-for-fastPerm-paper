# summarize results
library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

load("sym/symResultsExp_parallel_null.RData")

nVec <- unique(symResults$n)
alpha <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  sub <- symResults[symResults$n == nVec[i],
                   c("pPermAdj", "pPerm", "pBeta", "pPred", "pAsymNorm")]
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
 #      alpha  n pPermAdj pPerm pBeta pPred pAsymNorm
 # [1,]  0.01 20    0.009 0.009 0.008 0.016     0.054
 # [2,]  0.05 20    0.047 0.047 0.054 0.071     0.139
 # [3,]  0.10 20    0.107 0.107 0.109 0.141     0.192
 # [4,]  0.01 40    0.009 0.009 0.008 0.011     0.031
 # [5,]  0.05 40    0.048 0.048 0.043 0.069     0.105
 # [6,]  0.10 40    0.099 0.099 0.094 0.125     0.164
 # [7,]  0.01 60    0.007 0.007 0.007 0.010     0.029
 # [8,]  0.05 60    0.044 0.044 0.042 0.057     0.083
 # [9,]  0.10 60    0.091 0.091 0.093 0.119     0.142

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
  measure.vars=c("pPred", "pAsymNorm","pBeta"))

levels(plotData$variable) <- c("Alg 1","Asym", "Beta prime")

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
 # 0.0003  0.0015  0.0016  0.0019  0.0021  0.0048    2989 

# nonsymmetric sample sizes
load("nonSym/nonSymResultsExp_parallel_null.RData")


nVec <- unique(nonSymResults$nx)
alpha <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  sub <- nonSymResults[nonSymResults$nx == nVec[i],
                   c("pPermAdj", "pPerm", "pBeta", "pPred", "pAsymNorm")]
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
#   alpha nx pPermAdj pPerm pBeta pPred pAsymNorm
# 1  0.01 20    0.015 0.015 0.016 0.022     0.053
# 4  0.01 40    0.015 0.015 0.013 0.014     0.033
# 7  0.01 60    0.010 0.010 0.013 0.014     0.033
# 2  0.05 20    0.056 0.056 0.056 0.076     0.124
# 5  0.05 40    0.045 0.045 0.042 0.064     0.094
# 8  0.05 60    0.059 0.059 0.054 0.069     0.089
# 3  0.10 20    0.105 0.105 0.110 0.150     0.186
# 6  0.10 40    0.097 0.097 0.096 0.120     0.143
# 9  0.10 60    0.095 0.095 0.101 0.121     0.135

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
  measure.vars=c("pPred", "pAsymNorm", "pBeta"))

levels(plotData$variable) <- c("Alg 1","Asym", "Beta prime")

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
 # 0.0002  0.0006  0.0012  0.0019  0.0027  0.0054    2984 