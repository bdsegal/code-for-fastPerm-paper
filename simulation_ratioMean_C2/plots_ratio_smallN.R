# summarize results
library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symResultsAll <- symResults
# load("symResultsExp_parallel.RData")
# symResultsAll <- rbind(symResultsAll, symResults)

load("sym/symResultsExp_parallel_smallN_revision.RData")

head(symResults)

plotData <- melt(symResults, id.vars=c("n", "pBeta"), 
  measure.vars=c("pPred", "pAsymNorm", "pDeltaUneq", "pExpert6"))

levels(plotData$variable) <- c("Alg 1","Asym", "Delta", "SAMC")

dev.new(width=9, height=5)
ggplot(aes(x=log(pmin(pBeta, 1),10), y=log(value,10), color=as.factor(n),
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
ggsave(file.path(paperPath,"simExp_sym_smallN.png"))
ggsave(file.path(paperPath,"simExp_sym_smallN.tiff"))
ggsave(file.path(paperPath,"simExp_sym_smallN.pdf"))

dev.new(width=6, height=5)
ggplot(aes(x=log(pBeta,10), y=log(mStop*1000,10), color=as.factor(n),
  shape=as.factor(n)),
  data=symResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total resamples)")),
      x=expression(paste("lo",g[10],"(",p[beta],")")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simExp_sym_iter_smallN.png"))
ggsave(file.path(paperPath,"simExp_sym_iter_smallN.tiff"))
ggsave(file.path(paperPath,"simExp_sym_iter_smallN.pdf"))

summary(log(symResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -34.740 -13.770  -9.515 -10.800  -6.396  -2.449      15 

# nonsymmetric sample sizes
load("nonSym/nonSymResultsExp_parallel_smallN_revision.RData")

plotData <- melt(nonSymResults, id.vars=c("nx", "pBeta"), 
  measure.vars=c("pPred", "pAsymNorm","pDeltaUneq", "pExpert6"))
  # , "pExpert", 
    # "pSAMC6", "pSAMC5", "pSAMC4","pSAMC3"))

levels(plotData$variable) <- c("Alg 1","Asym", "Delta", "SAMC")

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
ggsave(file.path(paperPath,"simExp_nonSym_smallN.png"))
ggsave(file.path(paperPath,"simExp_nonSym_smallN.tiff"))
ggsave(file.path(paperPath,"simExp_nonSym_smallN.pdf"))

dev.new(width=6, height=5)
ggplot(aes(x=log(pBeta,10), y=log(mStop*1000,10), color=as.factor(nx),
  shape=as.factor(nx)),
  data=nonSymResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total resamples)")),
      x=expression(paste("lo",g[10],"(",p[beta],")")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simExp_nonSym_iter_smallN.png"))
ggsave(file.path(paperPath,"simExp_nonSym_iter_smallN.tiff"))
ggsave(file.path(paperPath,"simExp_nonSym_iter_smallN.pdf"))

summary(log(nonSymResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -34.610 -13.790 -10.340 -11.060  -7.482  -3.008 

# plot(y=log(nonSymResults$pExpert6,10), x=log(nonSymResults$pBeta,10))
# abline(a=0, b=1)