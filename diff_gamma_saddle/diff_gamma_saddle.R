library(ggplot2)
library(rmutil)
library(gammaDist)

if(length(grep("bdsegal",getwd()))>0 ) {
    computer <- "C:/Users/bdsegal/"
} else {
  computer <- "/home/bsegal/"
}
paperPath <- file.path(computer, "Dropbox/Research/PermTest/neighborhoods/paper")

marDefault <- c(5, 4, 4, 2) + 0.1 

# compare true against saddlepoint approximation ------------------------------
ny <- 100
nx <- 100
alpha <- 1
lambda <- 4

u <- seq(-5, 5, 0.01)
FHat <- pgammaDifSaddle(u, nx, ny, alpha, lambda)
Ftrue <- pgammaDif(u, nx, ny, alpha, lambda)

cbind(u = u, Ftrue = Ftrue, FHat = FHat)

png(file.path(paperPath, "FHat_Ftrue_comparison.png"))
par(mar = marDefault + c(0, 1.5, 0, 0))
plot(x = log(Ftrue, 10), y = log(FHat, 10),
     ylab = expression(paste("lo", g[10], "(", hat(G), ")", sep = "")),
     xlab = expression(paste("lo", g[10], "(G)", sep = "")),
     cex.axis = 1.7, cex.lab = 1.7, cex.main = 1.7)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()

tiff(file.path(paperPath, "FHat_Ftrue_comparison.tiff"))
par(mar = marDefault + c(0, 1.5, 0, 0))
plot(x = log(Ftrue, 10), y = log(FHat, 10),
     ylab = expression(paste("lo", g[10], "(", hat(G), ")", sep = "")),
     xlab = expression(paste("lo", g[10], "(G)", sep = "")),
     cex.axis = 1.7, cex.lab = 1.7, cex.main = 1.7)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()

pdf(file.path(paperPath, "FHat_Ftrue_comparison.pdf"))
par(mar = marDefault + c(0, 1.5, 0, 0))
plot(x = log(Ftrue, 10), y = log(FHat, 10),
     ylab = expression(paste("lo", g[10], "(", hat(G), ")", sep = "")),
     xlab = expression(paste("lo", g[10], "(G)", sep = "")),
     cex.axis = 1.7, cex.lab = 1.7, cex.main = 1.7)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
dev.off()