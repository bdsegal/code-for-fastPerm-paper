library(ggplot2)

if(length(grep("bdsegal",getwd()))>0 ) {
    computer <- "C:/Users/bdsegal/"
} else {
  computer <- "/home/bsegal/"
}
paperPath <- file.path(computer, "Dropbox/Research/PermTest/neighborhoods/paper")

marDefault <- c(5, 4, 4, 2) + 0.1 

source("cancerFunctions.R")
load("cancerAnalysis.RData")

# delta method with plug in estimates -----------------------------------------
pDelta <- rep(NA, nrow(LUADmatAbove))
names(pDelta) <- rownames(LUADmatAbove)
pDelta <- data.frame(geneNames = rownames(LUADmatAbove), r = NA, pDelta = NA)
for (i in 1:nrow(LUADmatAbove)) {
  x <- LUADmatAbove[i, ]
  y <- LUSCmatAbove[i, ]

  nx <- length(x)
  ny <- length(y)
  xbar <- mean(x)
  ybar <- mean(y)
  # Assuming unequal variance
  sdx <- sd(x)
  sdy <- sd(y)

  r <- xbar / ybar
  pDelta$r[i] <- r
  vx <- sqrt((sdx^2 / nx) / ybar^2 + (sdy^2 / ny) * (xbar^2 / ybar^4))
  vy <- sqrt((sdy^2 / ny) / xbar^2 + (sdx^2 / nx) * (ybar^2 / xbar^4))
  
  if(r >= 1) {
    pDelta$pDelta[i] <- pnorm(r, mean = 1, sd = vx, lower.tail = FALSE) +
                 pnorm(1 / r, mean = 1, sd = vy, lower.tail = TRUE)
  } else {
    pDelta$pDelta[i] <- pnorm(r, mean = 1, sd = vx, lower.tail = TRUE) +
                 pnorm(1 / r, mean = 1, sd = vy, lower.tail = FALSE)
  }
}

# fastPerm results ------------------------------------------------------------
fp <- result[, c("geneNames", "pPred")]
all <- merge(fp, pDelta, by = "geneNames")

# all <- merge(rbind(fp, MC), pDelta, by = "geneNames")
# names(pMean) <- rownames(LUADmatAbove)
# MC <- data.frame(geneNames = names(pMean[genesNotToTest]), pPred = pMean[genesNotToTest])
# with(all, plot(y = pDelta, x = pPred)) 
# abline(a = 0, b = 1, col = "red", lty = 2)

png(file.path(paperPath, "pDelta_cancer_app.png"))
par(mar = marDefault + c(0, 1, 0, 0))
with(all, plot(y = log(pDelta, 10), x = log(pPred, 10), 
               ylab = expression(paste("lo", g[10], "(", p[Delta], ")", sep = "")), 
               xlab = expression(paste("lo", g[10], "(", tilde(p)["pred"], ")", sep = "")),
               cex.axis = 1.8, cex.lab = 1.8))
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()

tiff(file.path(paperPath, "pDelta_cancer_app.tiff"))
par(mar = marDefault + c(0, 1, 0, 0))
with(all, plot(y = log(pDelta, 10), x = log(pPred, 10), 
               ylab = expression(paste("lo", g[10], "(", p[Delta], ")", sep = "")), 
               xlab = expression(paste("lo", g[10], "(", tilde(p)["pred"], ")", sep = "")),
               cex.axis = 1.8, cex.lab = 1.8))
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()

pdf(file.path(paperPath, "pDelta_cancer_app.pdf"))
par(mar = marDefault + c(0, 1, 0, 0))
with(all, plot(y = log(pDelta, 10), x = log(pPred, 10), 
               ylab = expression(paste("lo", g[10], "(", p[Delta], ")", sep = "")), 
               xlab = expression(paste("lo", g[10], "(", tilde(p)["pred"], ")", sep = "")),
               cex.axis = 1.8, cex.lab = 1.8))
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()

all <- all[order(all$pDelta), ]
zhanTop <- c("DSG3", "KRT5", "DSC3", "CALML3", "SERPINB13",
             "KRT6B", "KRT6C", "KRT6A", "PVRL1", "LOC642587",
             "PERP", "TP63", "TRIM29", "ATP1B3", "FAT2",
             "MLPH", "SFTA2", "TMC5", "SFTA3", "DDAH1",
             "RORC", "TMEM125", "SMPDL3B", "ALDH3B1", "ACSL5",
             "NKX2", "ATP11A", "CGN", "FMO5", "MUC1")

overlap <- unlist(sapply(zhanTop, function(x) {grep(x, all$geneNames[1:100])}))
overlap
 # PVRL1   PERP ATP1B3 
 #    83     86     18 

all[overlap, ]

# vs MC for large p-values ----------------------------------------------------
names(pMean) <- rownames(LUADmatAbove)
MC <- data.frame(geneNames = names(pMean[genesNotToTest]), pMC = pMean[genesNotToTest])
largeP <- merge(MC, pDelta, by = "geneNames")

nrow(MC)
nrow(fp)
png(file.path(paperPath, "pDelta_pMC_cancer_app.png"))
par(mar = marDefault + c(0, 1, 0, 0))
with(largeP, plot(pMC, pDelta, xlab = expression(tilde(p)), ylab = expression(p[Delta]),
                  cex.axis = 1.8, cex.lab = 1.8))
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()

pdf(file.path(paperPath, "pDelta_pMC_cancer_app.pdf"))
par(mar = marDefault + c(0, 1, 0, 0))
with(largeP, plot(pMC, pDelta, xlab = expression(tilde(p)), ylab = expression(p[Delta]),
                  cex.axis = 1.8, cex.lab = 1.8))
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()

tiff(file.path(paperPath, "pDelta_pMC_cancer_app.tiff"))
par(mar = marDefault + c(0, 1, 0, 0))
with(largeP, plot(pMC, pDelta, xlab = expression(tilde(p)), ylab = expression(p[Delta]),
                  cex.axis = 1.8, cex.lab = 1.8))
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()


all$rGT1 <- all$r >= 1
qplot(x = pMC, y = pDelta, color = rGT1, data = all)+
geom_abline(intercept = 0, slope = 1)


# # example of how wide the normal is, and whether the tails are cut off at 0
# example <- data.frame(x = seq(-1, 3, 0.01))
# example$p <- dnorm(example$x, mean = 1, sd = vy)
# plot(example, type = "l")
# hist(pDelta)
# pDelta <- sort(pDelta)