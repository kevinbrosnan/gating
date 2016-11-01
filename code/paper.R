#### Code used for Table and Figures in Gating Paper ####

  library(flowClust)

## Rituximab Data (Figure 1-2, Table 1) ##
  
  rituximab <- read.csv("https://raw.githubusercontent.com/significantstats/automatedgating/master/rituximab.csv")
  
  # Figure 1 - FSC v SSC
  Lo.initial <- flowClust(rituximab, varNames = c('FSC.H', 'SSC.H'), K = 1, 
                            B = 100, z.cutoff = 0.5)
  
  MRF.initial <- mrf_gating(rituximab[, c("FSC.H", "SSC.H")])
  
  par(mfrow = c(1,2))
  par(pty = 's')
  plot(Lo.initial, data = rituximab, xlab = 'FSC-Height', ylab = 'SSC-Height', 
       las = 1, col = 'white', main = "(a) t mixture with Box-Cox", show.outliers = FALSE)
  points(rituximab[which(Lo.initial@flagOutliers), c(1,2)], col = gray(3/4), pch = 20)
  
  par(pty = 's')
  plot(MRF.initial$x, type = "n", las = 1)
  points(MRF.initial$x[which(MRF.initial$removals == 1),], pch = "*")
  
  par(mfrow = c(1, 1))

  
  # Figure 2
  
  rit.reduced <- rituximab[rituximab %in% Lo.initial, c('FL3.H', 'FL1.H')]
  Lo.7AAD_antiBrdU <- flowClust(rit.reduced, K = 3, B = 100)
  
  par(mfrow = c(1,2))
  par(pty = 's')
  plot(Lo.7AAD_antiBrdU, data = rit.reduced, ylab = 'Anti-BrdU FITC', 
       xlab = '7 AAD', main = "(a) t mixture with Box-Cox", col = 'white', 
       las = 1, show.outliers = FALSE)
  points(rit.reduced[which(Lo.7AAD_antiBrdU@flagOutliers),], col = gray(3/4), pch = 20)
  
  par(pty = 's')
  
  par(mfrow = c(1, 1))
  

## GvHD Data (Figure 3-4, Table 2) ##
  
  GvHD.pos <- read.csv("https://raw.githubusercontent.com/significantstats/automatedgating/master/GvHD_positive.csv")
  GvHD.con <- read.csv("https://raw.githubusercontent.com/significantstats/automatedgating/master/GvHD_control.csv")  
  