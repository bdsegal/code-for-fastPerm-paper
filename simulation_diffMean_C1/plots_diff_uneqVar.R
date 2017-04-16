# Plot results

library(reshape2)
library(ggplot2)
library(dplyr)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symmetric sample size
load("sym/symResultsDiff_parallel_uneqVar.RData")

head(symResults)

plotData <- melt(symResults, id.vars=c("mux","n", "pt"), 
  measure.vars=c("pPred", "pPredStud", "pAsymNorm", "pExpert6", "pExpert6uneqVar"))

levels(plotData$variable) <- c("Alg 1", "Alg Student", "Asym","SAMC", "SAMC student")

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
ggsave(file.path(paperPath,"simDiff_sym_smallN.png"))

# non-symmetrblic sample size
load(file="nonSym/nonSymResultsDiff_parallel_uneqVar.RData")

head(nonSymResults)

# add SAMC with studentized statistic
plotData <- melt(nonSymResults, id.vars=c("nx", "pt"), 
  measure.vars=c("pPredStud", "pPred", "pAsymNorm", "pExpert6"))

levels(plotData$variable) <- c("Alg 1 studentized", "Alg 1", "Asymptotic", "SAMC")

dev.new(width=10, height=7)
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
ggsave(file.path(paperPath,"simDiff_nonSym_uneqVar.png"))