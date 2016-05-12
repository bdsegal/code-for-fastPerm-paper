# summarize results
library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symResultsAll <- symResults
# load("symResultsExp_parallel.RData")
# symResultsAll <- rbind(symResultsAll, symResults)

load("sym/symResultsExp_parallel.RData")

head(symResults)

plotData <- melt(symResults, id.vars=c("n", "pBeta"), 
  measure.vars=c("pPred", "pAsymNorm","pExpert6"))

levels(plotData$variable) <- c("Alg 1","Asym", "SAMC")

dev.new(width=9, height=5)
ggplot(aes(x=log(pBeta,10), y=log(value,10), color=as.factor(n),
    shape=as.factor(n)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[beta],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")
  # theme(axis.text.x = element_text(angle = 45, hjust = 0.9))

ggsave(file.path(paperPath,"simExp_sym.png"))

dev.new(width=6, height=5)
ggplot(aes(x=log(pBeta,10), y=log(mStop*1000,10), color=as.factor(n),
  shape=as.factor(n)),
  data=symResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(paste("lo",g[10],"(",p[beta],")")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  geom_hline(yintercept=log(1e6+5e4,10))

ggsave(file.path(paperPath,"simExp_sym_iter.png"))

summary(log(symResults$pExpert6,10))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -64.730 -22.000 -10.870 -16.570  -5.857  -2.522     239 

# nonsymmetric sample sizes
load("nonSym/nonSymResultsExp_parallel.RData")

plotData <- melt(nonSymResults, id.vars=c("nx", "pBeta"), 
  measure.vars=c("pPred", "pAsymNorm","pExpert6"))
  # , "pExpert", 
    # "pSAMC6", "pSAMC5", "pSAMC4","pSAMC3"))

levels(plotData$variable) <- c("Alg 1","Asym","SAMC")

dev.new(width=9, height=5)
ggplot(aes(x=log(pBeta,10), y=log(value,10), color=as.factor(nx),
    shape=as.factor(nx)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(~variable)+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[beta],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))
  # theme(axis.text.x = element_text(angle = 45, hjust = 0.9))

ggsave(file.path(paperPath,"simExp_nonSym.png"))

dev.new(width=6, height=5)
ggplot(aes(x=log(pBeta,10), y=log(mStop*1000,10), color=as.factor(nx),
  shape=as.factor(nx)),
  data=nonSymResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(paste("lo",g[10],"(",p[beta],")")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  geom_hline(yintercept=log(1e6+5e4,10))

ggsave(file.path(paperPath,"simExp_nonSym_iter.png"))

summary(log(nonSymResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -51.410 -17.870 -11.700 -14.030  -6.419  -2.443      30  

plot(y=log(nonSymResults$pExpert6,10), x=log(nonSymResults$pBeta,10))
abline(a=0, b=1)