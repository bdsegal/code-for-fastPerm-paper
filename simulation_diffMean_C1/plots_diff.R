# Plot results

library(reshape2)
library(ggplot2)
library(dplyr)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symmetric sample size
load("sym/symResultsDiff_parallel.RData")

plotData <- melt(symResults, id.vars=c("mux","n", "pt"), 
  measure.vars=c("pPred", "pAsymNorm", "pExpert6"))

levels(plotData$variable) <- c("Alg 1","Asym","SAMC")

dev.new(width=9, height=5)
ggplot(aes(x=log(pt,10), y=log(value,10), color=as.factor(n),
    shape=as.factor(n)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[t],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  # scale_x_continuous(breaks=seq(0,-120,-30))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simDiff_sym.png"))

symResults %>% group_by(n, mux) %>%
  summarize(mStopMean = mean(mStop))

dev.new(width=6, height=5)
ggplot(aes(x=log(pt,10), y=log(mStop*1000,10), color=as.factor(n),
  shape=as.factor(n)),
  data=symResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(paste("lo",g[10],"(",p[t],")")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simDiff_sym_iter.png"))

summary(log(symResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -32.230 -11.600  -8.948 -10.090  -6.423  -2.650     385 

# non-symmetrblic sample size
load(file="nonSym/nonSymResultsDiff_parallel.RData")

plotData <- melt(nonSymResults, id.vars=c("nx", "pt"), 
  measure.vars=c("pPred","pAsymNorm", "pExpert6"))

levels(plotData$variable) <- c("Alg 1","Asym","SAMC")

dev.new(width=9, height=5)
ggplot(aes(x=log(pt,10), y=log(value,10), color=as.factor(nx),
  shape=as.factor(nx)
  ), data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_wrap(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[t],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simDiff_nonSym.png"))

nonSymResults %>% group_by(nx) %>%
  summarize(mStopMean = mean(mStop))

dev.new(width=6, height=5)
ggplot(aes(x=log(pt,10), y=log(mStop*100,10), color=as.factor(nx),
  shape=as.factor(nx)), data=nonSymResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(paste("lo",g[10],"(",p[t],")")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simDiff_nonSym_iter.png"))

summary(log(nonSymResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -36.780 -22.710 -14.860 -16.370  -8.698  -2.541     179 
