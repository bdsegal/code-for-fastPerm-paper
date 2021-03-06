# summarize results
library(ggplot2)
library(reshape2)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"
# paperPath <- "C:/Users/bdsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# symResultsAll <- symResults
# load("symResultsExp_parallel.RData")
# symResultsAll <- rbind(symResultsAll, symResults)

load("sym/symResultsExp_parallel_revision.RData")

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
# ggsave("ratio_means.png")
  # theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"simExp_sym.png"))
ggsave(file.path(paperPath,"simExp_sym.tiff"))
ggsave(file.path(paperPath,"simExp_sym.pdf"))
ggsave(file.path(paperPath,"figure_5a.png"))
ggsave(file.path(paperPath,"figure_5a.tiff"))
ggsave(file.path(paperPath,"figure_5a.pdf"))

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
ggsave(file.path(paperPath,"simExp_sym_iter.png"))
ggsave(file.path(paperPath,"simExp_sym_iter.tiff"))
ggsave(file.path(paperPath,"simExp_sym_iter.pdf"))
ggsave(file.path(paperPath,"figure_5b.png"))
ggsave(file.path(paperPath,"figure_5b.tiff"))
ggsave(file.path(paperPath,"figure_5b.pdf"))

summary(log(symResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -63.530 -24.910 -11.770 -17.610  -5.795  -2.518     246 


# # for presentation
# plotData <- melt(symResults, id.vars=c("n", "pBeta"), 
#   measure.vars=c("pPred", "pExpert6"))

# levels(plotData$variable) <- c("Our method", "SAMC (competitor)")

# dev.new(width=9, height=5)
# ggplot(aes(x=log(pmin(pBeta, 1),10), y=log(value,10), color=as.factor(n),
#     shape=as.factor(n)),
#   data=plotData)+
#   geom_point(size=1.5)+
#   theme_bw(28)+
#   facet_grid(~variable)+
#   geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
#   labs(x=expression(paste("lo",g[10],"(",p[beta],")",sep="")),
#     y=expression(paste("lo",g[10],"(p)",sep="")))+
#   scale_color_discrete("n")+
#   scale_shape_discrete("n")
#   # theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
# ggsave(file.path("/home/bsegal/Dropbox/Research/job_talk_2017/plots","simExp_sym.png"))

# dev.new(width=6, height=5)
# ggplot(aes(x=log(pBeta,10), y=log(mStop*1000,10), color=as.factor(n),
#   shape=as.factor(n)),
#   data=symResults)+
#   geom_point(size=1.5)+
#   theme_bw(28)+
#   labs(y=expression(paste("lo",g[10],"(total resamples)")),
#       x=expression(paste("lo",g[10],"(",p[beta],")")))+
#   scale_color_discrete("n")+
#   scale_shape_discrete("n")+
#   geom_hline(yintercept=log(1e6+5e4,10))
# ggsave(file.path("/home/bsegal/Dropbox/Research/job_talk_2017/plots","simExp_sym_iter.png"))


# nonsymmetric sample sizes
load("nonSym/nonSymResultsExp_parallel_revision.RData")

plotData <- melt(nonSymResults, id.vars=c("nx", "pBeta"), 
  measure.vars=c("pPred", "pAsymNorm", "pDeltaUneq","pExpert6"))
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
ggsave(file.path(paperPath,"simExp_nonSym.png"))
ggsave(file.path(paperPath,"simExp_nonSym.tiff"))
ggsave(file.path(paperPath,"simExp_nonSym.pdf"))
ggsave(file.path(paperPath,"figure_6a.png"))
ggsave(file.path(paperPath,"figure_6a.tiff"))
ggsave(file.path(paperPath,"figure_6a.pdf"))

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
ggsave(file.path(paperPath,"simExp_nonSym_iter.png"))
ggsave(file.path(paperPath,"simExp_nonSym_iter.tiff"))
ggsave(file.path(paperPath,"simExp_nonSym_iter.pdf"))
ggsave(file.path(paperPath,"figure_6b.png"))
ggsave(file.path(paperPath,"figure_6b.tiff"))
ggsave(file.path(paperPath,"figure_6b.pdf"))

summary(log(nonSymResults$pExpert6,10))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -57.960 -18.080 -11.900 -13.880  -6.552  -2.556      33 
