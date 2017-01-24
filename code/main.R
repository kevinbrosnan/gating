#------------------------------------------------------------------------------#
#                  Automated Gating of Flow Cytometry Data via                 #
#                        Markov Random Fields Clustering                       #
#                                                                              #
# Authors: Mr. Kevin Brosnan, Maths & Stats, University of Limerick            #
#          Dr. Norma Bargary, Maths & Stats, University of Limerick            #
#          Dr. Kevin Hayes,   Maths & Stats, University of Limerick            #
#                                                                              #
# Cytometry Part A:                                                            #
#          Submitted: dd/mm/yyyy                                               #
#          Accepted:  dd/mm/yyyy                                               #
#          Published: dd/mm/yyyy                                               #
#                                                                              #
#------------------------------------------------------------------------------#

####--- Set Up -------------------------------------------------------------####

  # Required Libraries
  # install.packages(c('stargazer', 'SDMTools'))
  library(SDMTools)
  library(stargazer)
  library(flowClust)

  # GitHub Access
  gh.proj <- "https://raw.githubusercontent.com/significantstats/gating/master/"
  gh.code <- paste0(gh.proj, "code/functions/")
  gh.data <- paste0(gh.proj, "data/")
  
  # Required Data
  GvHD.pos <- read.csv(paste0(gh.data, "GvHD/GvHD_positive.csv"))
  GvHD.con <- read.csv(paste0(gh.data, "GvHD/GvHD_control.csv"))
  rituximab <- read.csv(paste0(gh.data, "rituximab/rituximab.csv"))
  
  # Required Functions
  source(paste0(gh.code, "energy.R"))
  source(paste0(gh.code, "energy_system.R"))
  source(paste0(gh.code, "explore_flow.R"))
  source(paste0(gh.code, "grid_inc.R"))
  source(paste0(gh.code, "grid_red.R"))
  source(paste0(gh.code, "ising_model.R"))
  source(paste0(gh.code, "make_grid.R"))
  source(paste0(gh.code, "mrf_gating.R"))
  source(paste0(gh.code, "neigh_system.R"))
  source(paste0(gh.code, "neighbours.R"))
  source(paste0(gh.code, "unmake_grid.R"))
  
  # Colour Schemes
  cluster.colours <- c("red", "green", "blue", "purple", "orange", "yellow")
  
####--- Exploratory Analysis (Rituximab Data) ------------------------------####
  
  # Manual Gating Example
  par(pty = "s")
  plot(rituximab[ , c("FSC.H", "SSC.H")], xlab = "FSC-Height", 
       ylab = "SSC-Height", las = 1, col = "gray")
  polygon(x = c(50,50,440,380), y = c(0,200,580,80), 
          border = "red", lwd = 2, lty = 1)
  polygon(x = c(50,50,420,600), y = c(0,200,610,110), 
          border = "blue", lwd = 2, lty = 2)
  legend("topright", legend = c("Expert 1", "Expert 2"), 
         col = c("red", "blue"), lty = c(1,2), lwd = c(2,2))
  par(pty = "m")
  
  # t-mixtures approach
  model <- flowClust(rituximab, varNames = c("FSC.H", "SSC.H"), K = 1:6)
  best <- which.max(criterion(object = model, "BIC"))
  plot(model[[best]], data = rituximab, las = 1,
       xlab = "FSC-Height", ylab = "SSC-Height")
  
  # Descriptive Table & Plots
  explore_flow(rituximab)
  stargazer(rituximab, summary = TRUE, digits = 2, label = "tab:ritdatasum", 
            median = TRUE, iqr = TRUE, mean.sd = FALSE,
            title = "Summary Statistics for Rituximab Data", align = TRUE)
  
####--- Figures & Tables from Paper ----------------------------------------####
  
  # Figure 1 - Rituximab Data FSC-Height v SSC-Height
  Lo.initial <- flowClust(rituximab, varNames = c("FSC.H", "SSC.H"), K = 1, 
                          B = 100, z.cutoff = 0.5)
  
  MRF.initial <- mrf_gating(rituximab[, c("FSC.H", "SSC.H")], temperature = 4)
  
  par(mfrow = c(1,2))
  par(pty = "s")
  plot(Lo.initial, data = rituximab, xlab = "FSC-Height", 
       ylab = "SSC-Height", las = 1, col = 'white', 
       main = "(a) t mixture with Box-Cox", show.outliers = FALSE)
  points(rituximab[which(Lo.initial@flagOutliers), c(1,2)], 
         col = gray(3/4), pch = 20)
  
  par(pty = "s")
  
  plot(MRF.initial$x, type = "n", las = 1)
  points(MRF.initial$x[which(MRF.initial$removals == 1),], pch = "*")
  points(MRF.initial$x[which(MRF.initial$groups != 0),], col = MRF.initial$groups)
  
  par(mfrow = c(1, 1))
  
  # Figure 2 - Rituximab Data 7-AAD v Anti-BrdU
  Lo.gate <- rituximab[rituximab %in% Lo.initial, c("FL3.H", "FL1.H")]
  Lo.7AAD_antiBrdU <- flowClust(Lo.gate, K = 3, B = 100)
  
  MRF.gate <- rituximab[which(MRF.initial$groups == 1), c("FL3.H", "FL1.H")]
  MRF.7AAD_antiBrdU <- mrf_gating(MRF.gate, temperature = 4)
  
  par(mfrow = c(1,2))
  par(pty = "s")
  plot(Lo.7AAD_antiBrdU, data = Lo.gate, ylab = 'Anti-BrdU FITC', 
       xlab = '7 AAD', main = "(a) t mixture with Box-Cox", col = 'white', 
       las = 1, show.outliers = FALSE)
  points(Lo.gate[which(Lo.7AAD_antiBrdU@flagOutliers),], 
         col = gray(3/4), pch = 20)
  
  par(pty = 's')
  plot(MRF.7AAD_antiBrdU$x, type = "n", las = 1)
  points(MRF.7AAD_antiBrdU$x[which(MRF.7AAD_antiBrdU$removals == 1),], pch = "*")
  points(MRF.7AAD_antiBrdU$x[which(MRF.7AAD_antiBrdU$groups == 9),], col = 'red')
  par(mfrow = c(1, 1))
  
  
  # Table 1 - Rituximab Data
  
  
  # Figure 3 - GvHD Data Control Group CD4 v CD8B, CD4 v CD3
  
  
  # Figure 4 - GvHD Data Positive Group CD4 v CD8B, CD4 v CD3
  
  
  # Table 2 - GvHD Data
  
####--- End of Script ------------------------------------------------------####  