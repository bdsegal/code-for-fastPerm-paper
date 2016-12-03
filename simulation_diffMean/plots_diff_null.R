# Plot results

library(reshape2)
library(ggplot2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symmetric sample size
load("sym/symResultsDiff_parallel_null.RData")

nVec <- unique(symResults$n)
alpha <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  sub <- symResults[symResults$n == nVec[i],
                   c("pPermAdj", "pPerm", "pt", "pPred", "pAsymNorm")]
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
#   alpha  n pPermAdj pPerm    pt pPred pAsymNorm
# 1  0.01 20    0.010 0.010 0.010 0.015     0.010
# 4  0.01 40    0.013 0.013 0.013 0.015     0.013
# 7  0.01 60    0.010 0.010 0.010 0.011     0.010
# 2  0.05 20    0.048 0.048 0.050 0.064     0.050
# 5  0.05 40    0.055 0.055 0.055 0.075     0.056
# 8  0.05 60    0.049 0.049 0.050 0.061     0.050
# 3  0.10 20    0.098 0.098 0.098 0.137     0.105
# 6  0.10 40    0.109 0.109 0.110 0.140     0.112
# 9  0.10 60    0.103 0.103 0.100 0.124     0.102

falseAlarmTabM <- melt(falseAlarmTab, id.vars = c("n", "alpha"))
ggplot(aes(y = value, x = variable), data = falseAlarmTabM)+
  geom_bar(stat = "identity")+
  facet_grid(alpha ~ n)+
  theme_bw(15)+
  geom_hline(aes(yintercept = alpha), linetype = "dashed",
             color = "red", size = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))+
  labs(x = "", y = "Type I error rate")

plotData <- melt(symResults, id.vars=c("mux","n", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "pt"))

levels(plotData$variable) <- c("Alg 1","Asym","t-test")

dev.new(width=9, height=5)
ggplot(aes(x=pPerm, y=value, color=as.factor(n),
    shape=as.factor(n)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)), y="p")+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simDiff_sym_null.png"))

symResults %>% group_by(n, mux) %>%
  summarize(mStopMean = mean(mStop))

dev.new(width=6, height=5)
ggplot(aes(x=pt, y=log(mStop*1000,10), color=as.factor(n),
  shape=as.factor(n)),
  data=symResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(p[t]))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simDiff_sym_iter_null.png"))

summary(symResults$pExpert6)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 # 0.0001  0.0012  0.0020  0.0016  0.0021  0.0028    2989 

# non-symmetrblic sample size
load(file="nonSym/nonSymResultsDiff_parallel_null.RData")

nVec <- unique(nonSymResults$nx)
alpha <- c(0.01, 0.05, 0.1)

falseAlarmTab <- NULL
for (i in 1:length(nVec)) {
  sub <- nonSymResults[nonSymResults$nx == nVec[i],
                   c("pPermAdj", "pPerm", "pt", "pPred", "pAsymNorm")]
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
#   alpha nx pPermAdj pPerm    pt pPred pAsymNorm
# 1  0.01 20    0.013 0.013 0.013 0.018     0.013
# 4  0.01 40    0.016 0.016 0.016 0.018     0.016
# 7  0.01 60    0.010 0.010 0.010 0.013     0.010
# 2  0.05 20    0.049 0.049 0.049 0.075     0.049
# 5  0.05 40    0.047 0.047 0.047 0.066     0.047
# 8  0.05 60    0.044 0.044 0.044 0.057     0.044
# 3  0.10 20    0.090 0.090 0.090 0.135     0.092
# 6  0.10 40    0.103 0.103 0.104 0.135     0.105
# 9  0.10 60    0.090 0.090 0.090 0.133     0.090

falseAlarmTabM <- melt(falseAlarmTab, id.vars = c("nx", "alpha"))
ggplot(aes(y = value, x = variable), data = falseAlarmTabM)+
  geom_bar(stat = "identity")+
  facet_grid(alpha ~ nx)+
  theme_bw()+
  geom_hline(aes(yintercept = alpha), linetype = "dashed",
             color = "red", size = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))+
  labs(x = "", y = "Type I error rate")


plotData <- melt(nonSymResults, id.vars=c("nx", "pPerm"), 
  measure.vars=c("pPred","pAsymNorm", "pt"))

levels(plotData$variable) <- c("Alg 1","Asym","t-test")

dev.new(width=9, height=5)
ggplot(aes(x=pPerm, y=value, color=as.factor(nx),
  shape=as.factor(nx)
  ), data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_wrap(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)), y= "p")+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simDiff_nonSym_null.png"))

nonSymResults %>% group_by(nx) %>%
  summarize(mStopMean = mean(mStop))

dev.new(width=6, height=5)
ggplot(aes(x=pt, y=log(mStop*1000,10), color=as.factor(nx),
  shape=as.factor(nx)), data=nonSymResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(p[t]))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simDiff_nonSym_iter_null.png"))

summary(nonSymResults$pExpert6)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 # 0.0002  0.0006  0.0010  0.0013  0.0019  0.0030    2986 
