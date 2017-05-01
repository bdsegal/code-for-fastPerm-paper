
# Simulated comparison with t-test
library(ggplot2)
library(reshape2)
library(mcc)

paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# small N
Bmc <- 1e5

# symmetric
load("mccOut_symGamma_smallN.Rdata")
mccOut$pMC[which(mccOut$pMC == 0)] <- NA
length(which(mccOut$pMC > 1/(Bmc*1e-2)))
# [1] 1019

mccOutM <- melt(mccOut[which(mccOut$pMC > 1/(Bmc*1e-2)),], 
                id.vars = c("n", "alpha", "pMC"),
                measure.vars = c("twosidedp", "pdouble"))
levels(mccOutM$variable) <- c("Two sided", "Double")#, "t-test")

dev.new(width=6.5, height=6.5)
ggplot(aes(x=log(pMC,10), y=log(value, 10), color=as.factor(n),
    shape=as.factor(n)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(alpha~variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",tilde(p),")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete("n")+
  scale_shape_discrete("n")
ggsave(file.path(paperPath,"mcc_symGammaDiff_smallN.png"))
ggsave(file.path(paperPath,"mcc_symGammaDiff_smallN.tiff"))
ggsave(file.path(paperPath,"mcc_symGammaDiff_smallN.pdf"))

# unsymmetric
load("mccOut_nonSymGamma_smallN.Rdata")

mccOut$pMC[which(mccOut$pMC == 0)] <- NA
length(which(mccOut$pMC > 1/(Bmc*1e-2)))
# [1] 705

mccOutM <- melt(mccOut[which(mccOut$pMC > 1/(Bmc*1e-2)),], 
                id.vars = c("nx", "alpha", "pMC"),
                measure.vars = c("twosidedp", "pdouble"))
levels(mccOutM$variable) <- c("Two sided", "Double")#, "t-test")

dev.new(width=6.5, height=6.5)
ggplot(aes(x=log(pMC,10), y=log(value, 10), color=as.factor(nx),
    shape=as.factor(nx)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(alpha~variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(paste("lo",g[10],"(",tilde(p),")",sep="")),
    y=expression(paste("lo",g[10],"(p)",sep="")))+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))
ggsave(file.path(paperPath,"mcc_nonSymGammaDiff_smallN.png"))
ggsave(file.path(paperPath,"mcc_nonSymGammaDiff_smallN.tiff"))
ggsave(file.path(paperPath,"mcc_nonSymGammaDiff_smallN.pdf"))

# null ------------------------------------------------------------------------
# symmetric
load("mccOut_symGamma_null.Rdata")
mccOut$pMC[which(mccOut$pMC == 0)] <- NA
length(which(mccOut$pMC > 1/(Bmc*1e-2)))
# [1] 17981
mccOutM <- melt(mccOut[which(mccOut$pMC > 1/(Bmc*1e-2)),], 
                id.vars = c("n", "alpha", "pMC"),
                measure.vars = c("twosidedp", "pdouble"))
levels(mccOutM$variable) <- c("Two sided", "Double")#, "t-test")

dev.new(width=6.5, height=6.5)
ggplot(aes(x=pMC, y=value, color=as.factor(n),
    shape=as.factor(n)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(alpha~variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)),
    y="p")+
  scale_color_discrete("n")+
  scale_shape_discrete("n")+
 theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"mcc_symGammaDiff_null.png"))
ggsave(file.path(paperPath,"mcc_symGammaDiff_null.tiff"))
ggsave(file.path(paperPath,"mcc_symGammaDiff_null.pdf"))

load("mccOut_nonSymGamma_null.Rdata")
mccOut$pMC[which(mccOut$pMC == 0)] <- NA
length(which(mccOut$pMC > 1/(Bmc*1e-2)))
# [1] 17980

mccOutM <- melt(mccOut[which(mccOut$pMC > 1/(Bmc*1e-2)),], 
                id.vars = c("nx", "alpha", "pMC"),
                measure.vars = c("twosidedp", "pdouble"))
levels(mccOutM$variable) <- c("Two sided", "Double")#, "t-test")

dev.new(width=6.5, height=6.5)
ggplot(aes(x=pMC, y=value, color=as.factor(nx),
    shape=as.factor(nx)), data=mccOutM)+
  geom_point()+
  theme_bw(26)+
  facet_grid(alpha~variable, labeller = label_bquote(alpha == .(alpha)))+
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
  labs(x=expression(tilde(p)),
    y="p")+
  scale_color_discrete(expression(n[x]))+
  scale_shape_discrete(expression(n[x]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
ggsave(file.path(paperPath,"mcc_nonSymGammaDiff_null.png"))
ggsave(file.path(paperPath,"mcc_nonSymGammaDiff_null.tiff"))
ggsave(file.path(paperPath,"mcc_nonSymGammaDiff_null.pdf"))
