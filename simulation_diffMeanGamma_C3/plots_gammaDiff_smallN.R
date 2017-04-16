# summarize results
library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symResultsAll <- symResults
# load("symResultsExp_parallel.RData")
# symResultsAll <- rbind(symResultsAll, symResults)

load("sym/symResultsGammaDiff_smallN.RData")

symResults$pPerm[which(symResults$pPerm == 0)] <- NA

B <- 10^5
keep <- which(symResults$pPerm > 1/(B*10^-2))
length(keep)
# [1] 1023

# using MC as baseline
plotData <- melt(symResults[keep, ], id.vars=c("n", "alpha", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "ptUneq", "gammaDiffSaddle"))#, "pExpert6"))

levels(plotData$variable) <- c("Alg 1","Asym", "t-test", "Saddle")#, "SAMC")

dev.new(width=9, height=6.5)
ggplot(aes(x=log(pmin(pPerm, 1),10), y=log(value,10), color=as.factor(n),
    shape=as.factor(n)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(alpha ~ variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",tilde(p),")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simGammaDiff_sym_smallN.png"))

# using saddle point as baseline
plotData <- melt(symResults, id.vars=c("n", "alpha", "gammaDiffSaddle"), 
  measure.vars=c("pPred", "pAsymNorm", "ptUneq", "pPerm", "pExpert6"))

levels(plotData$variable) <- c("Alg 1","Asym", "t-test", "MC", "SAMC")

dev.new(width=9, height=6.5)
ggplot(aes(x=log(pmin(gammaDiffSaddle, 1),10), y=log(value,10), color=as.factor(n),
    shape=as.factor(n)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(alpha ~ variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",p["saddle"],")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simGammaDiff_sym_smallN_saddle.png"))

dev.new(width=6, height=5)
ggplot(aes(x=log(gammaDiffSaddle,10), y=log(mStop*1000,10), color=as.factor(n),
  shape=as.factor(n)),
  data=symResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(paste("lo",g[10],"(",p["saddle"],")")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simGammaDiff_sym_iter_smallN.png"))

summary(log(symResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -11.810  -5.271  -4.037  -4.438  -3.191  -2.322     851 

# nonsymmetric sample sizes
load("nonSym/nonSymResultsGammaDiff_smallN.RData")

nonSymResults$pPerm[which(nonSymResults$pPerm == 0)] <- NA

B <- 10^5
keep <- which(nonSymResults$pPerm > 1/(B*10^-2))
length(keep)
# [1] 573

plotData <- melt(nonSymResults[keep, ], id.vars=c("nx", "alpha", "pPerm"), 
  measure.vars=c("pPred", "pAsymNorm", "ptUneq", "gammaDiffSaddle"))#, "pExpert6"))

levels(plotData$variable) <- c("Alg 1","Asym", "t-test", "Saddle")#, "SAMC")

dev.new(width=9, height=6.5)
ggplot(aes(x=log(pPerm,10), y=log(value,10), color=as.factor(nx),
    shape=as.factor(nx)),
  data=plotData)+
  geom_point(size=1.5)+
  theme_bw(24)+
  facet_grid(alpha ~ variable, labeller = label_bquote(alpha == .(alpha)))+  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",tilde(p),")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))
  # theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simGammaDiff_nonSym_smallN.png"))

dev.new(width=6, height=5)
ggplot(aes(x=log(pPerm,10), y=log(mStop*1000,10), color=as.factor(nx),
  shape=as.factor(nx)),
  data=nonSymResults)+
  geom_point(size=1.5)+
  theme_bw(24)+
  labs(y=expression(paste("lo",g[10],"(total iterations)")),
      x=expression(paste("lo",g[10],"(",p[beta],")")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  geom_hline(yintercept=log(1e6+5e4,10))
ggsave(file.path(paperPath,"simGammaDiff_nonSym_iter_smallN.png"))

summary(log(nonSymResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -11.500  -5.228  -3.830  -4.372  -3.100  -2.275     922 

# plot(y=log(nonSymResults$pExpert6,10), x=log(nonSymResults$pBeta,10))
# abline(a=0, b=1)