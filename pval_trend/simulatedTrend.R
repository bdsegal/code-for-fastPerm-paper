# Section 3, Figure 1a
# Simulated trend in p-values across the partitions

library(ggplot2)
paperPath <- "/home/bsegal/Dropbox/Research/PermTest/neighborhoods/paper"

# Functions -----------------------------------------------

pFun <- function(x,t0){
  mean(abs(x) >= abs(t0))
}

ratio_mean <- function(x,y){
  max(mean(x)/mean(y),mean(y)/mean(x))
}

pValPlotRatioFacet <- function(x,y){

  ratioMeanPartition <- list()
    
  B <- 10^3
  for (m in 1:nx){
    print(m)
    
    ratioMeanVec <- rep(NA,B)
    
    for (b in 1:B){
    
      indx <- sample(1:nx, size = m, replace = FALSE)
      indy <- sample(1:ny, size = m, replace = FALSE)

      xNew <- c(x[-indx], y[indy])
      yNew <- c(y[-indy], x[indx])
    
      ratioMeanVec[b] <- ratio_mean(xNew, yNew)
      
    }
    
    ratioMeanPartition[[m]] <- ratioMeanVec

  }

  # calculate p-values

  t0_ratio_mean=ratio_mean(x,y)
  
  p_ratio_mean <- c(1, 
                  sapply(ratioMeanPartition, pFun, t0_ratio_mean))
  
  return(p_ratio_mean)
}
  
# Plots ---------------------------------------------------

nx = 100
ny = 100

out <- list()

x <- rpois(n = nx, lambda = 4)
y <- rpois(n = ny, lambda = 2)
out[[1]] <- data.frame(p = pValPlotRatioFacet(x,y), dist="Poisson")

x <- rnbinom(n = nx, size = 3, mu = 2)
y <- rnbinom(n = ny, size = 3, mu = 4)
out[[2]] <- data.frame(p = pValPlotRatioFacet(x, y),
                    dist = "Negative binomial")
x <- rexp(n = nx, rate = 2)
y <- rexp(n = ny, rate = 1)
out[[3]] <- data.frame(p = pValPlotRatioFacet(x, y),
                    dist = "Exponential")

x <- rlnorm(n = nx, meanlog = 2, sdlog = 1)
y <- rlnorm(n = ny, meanlog = 1, sdlog = 1)
out[[4]] <- data.frame(p = pValPlotRatioFacet(x, y),
                    dist = "Log normal")

plotData <- do.call(rbind, out)
plotData$m <- 0:100

plotData$dist <- factor(plotData$dist, levels=c("Poisson",
                  "Negative binomial", "Exponential", "Log normal"))

dev.new(width = 14, height = 5)
ggplot(aes(x = m, log(p, 10)), data = plotData)+
  geom_point(size = 1.5)+
  theme_bw(24)+
  facet_wrap(~dist, nrow = 1)+
  labs(y = expression(paste("lo", g[10], "(p)", sep = "")),
       x = "Partition m")
ggsave(file.path(paperPath,"logp_ratio_mean_facet.png"))





