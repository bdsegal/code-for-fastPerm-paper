# summarize results
library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symResultsAll <- symResults
# load("symResultsExp_parallel.RData")
# symResultsAll <- rbind(symResultsAll, symResults)

load("sym/symResultsGamma_smallN.RData")

head(symResults)

plotData <- melt(symResults, id.vars=c("n", "alpha", "pBeta"), 
  measure.vars=c("pPred", "pAsymNorm", "pDeltaUneq", "pExpert6"))

levels(plotData$variable) <- c("Alg 1","Asym", "Delta", "SAMC")
# plotData$alpha <- paste("alpha = ", plotData$alpha, sep = "")

dev.new(width=9, height=6.5)
ggplot(aes(x=log(pmin(pBeta, 1),10), y=log(value,10), color=as.factor(n),
    shape=as.factor(n)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(alpha ~ variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[beta],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")
  # theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simGamma_sym_smallN.png"))
ggsave(file.path(paperPath,"simGamma_sym_smallN.tiff"))
ggsave(file.path(paperPath,"simGamma_sym_smallN.pdf"))

dev.new(width=6, height=5)
ggplot(aes(x=log(pBeta,10), y=log(mStop*1000,10), color=as.factor(n),
  shape=as.factor(n)),
  data=symResults)+
  geom_point(size=1.5)+
  theme_bw(28)+
  labs(y=expression(paste("lo",g[10],"(total resamples)")),
      x=expression(paste("lo",g[10],"(",p[beta],")")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simGamma_sym_iter_smallN.png"))
ggsave(file.path(paperPath,"simGamma_sym_iter_smallN.tiff"))
ggsave(file.path(paperPath,"simGamma_sym_iter_smallN.pdf"))

summary(log(symResults$pExpert6,10))
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 # -39.10  -11.01  -10.03  -10.79   -7.23   -2.45     652 

# nonsymmetric sample sizes
load("nonSym/nonSymResultsGamma_smallN.RData")

head(nonSymResults)

plotData <- melt(nonSymResults, id.vars=c("nx", "alpha", "pBeta"), 
  measure.vars=c("pPred", "pAsymNorm","pDeltaUneq", "pExpert6"))
  # , "pExpert", 
    # "pSAMC6", "pSAMC5", "pSAMC4","pSAMC3"))

levels(plotData$variable) <- c("Alg 1","Asym", "Delta", "SAMC")

dev.new(width=9, height=6.5)
ggplot(aes(x=log(pBeta,10), y=log(value,10), color=as.factor(nx),
    shape=as.factor(nx)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(alpha ~ variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p[beta],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))
  # theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simGamma_nonSym_smallN.png"))
ggsave(file.path(paperPath,"simGamma_nonSym_smallN.tiff"))
ggsave(file.path(paperPath,"simGamma_nonSym_smallN.pdf"))

dev.new(width=6, height=5)
ggplot(aes(x=log(pBeta,10), y=log(mStop*1000,10), color=as.factor(nx),
  shape=as.factor(nx)),
  data=nonSymResults)+
  geom_point(size=1.5)+
  theme_bw(28)+
  labs(y=expression(paste("lo",g[10],"(total resamples)")),
      x=expression(paste("lo",g[10],"(",p[beta],")")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simGamma_nonSym_iter_smallN.png"))
ggsave(file.path(paperPath,"simGamma_nonSym_iter_smallN.tiff"))
ggsave(file.path(paperPath,"simGamma_nonSym_iter_smallN.pdf"))

summary(log(nonSymResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -36.500 -12.040  -8.506  -9.183  -5.062  -2.212     304 
