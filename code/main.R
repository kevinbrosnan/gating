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
  source(paste0(gh.code, "probability_plot.R"))
  source(paste0(gh.code, "spatial_smooth.R"))
  
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
  
  png('paper/figures/rituximab.png', height = 500, width = 1500)
  par(mfrow = c(1,3))
  par(pty = "s")
  plot(Lo.initial, data = rituximab, xlab = "FSC-Height", 
       ylab = "SSC-Height", las = 1, 
       main = "(a) t mixture with Box-Cox", show.outliers = TRUE,
       pch.outliers = 20)
  plot(MRF.initial, xlab = "FSC-Height", ylab = "SSC-Height",
  		main = "(b) Markov random field probabilities")
  plot(MRF.initial$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(c) Markov random field clusters",
       col = ifelse(MRF.initial$groups == 0, "grey", MRF.initial$groups))
  par(mfrow = c(1, 1))
  dev.off()
  
  # Figure 2 - Rituximab Data 7-AAD v Anti-BrdU
  Lo.gate <- rituximab[rituximab %in% Lo.initial, c("FL3.H", "FL1.H")]
  Lo.7AAD_antiBrdU <- flowClust(Lo.gate, K = 3, B = 100)
  
  MRF.gate <- rituximab[which(MRF.initial$groups == 1), c("FL3.H", "FL1.H")]
  MRF.7AAD_antiBrdU <- mrf_gating(MRF.gate, temperature = 4)
  
  png('paper/figures/rituximab_stagetwo.png', height = 500, width = 1500)
  par(mfrow = c(1,3))
  par(pty = "s")
  plot(Lo.7AAD_antiBrdU, data = Lo.gate, xlab = "FSC-Height", 
       ylab = "SSC-Height", las = 1, 
       main = "(a) t mixture with Box-Cox", show.outliers = TRUE,
       pch.outliers = 20)
  plot(MRF.7AAD_antiBrdU, xlab = "FSC-Height", ylab = "SSC-Height",
       main = "(b) Markov random field probabilities")
  plot(MRF.7AAD_antiBrdU$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(c) Markov random fields clusters",
       col = ifelse(MRF.7AAD_antiBrdU$groups == 0, "grey", MRF.7AAD_antiBrdU$groups))
  par(mfrow = c(1, 1))
  dev.off()
  
  # Table 1 - Rituximab Data
  
  
  # Figure 3 - GvHD Data Control Group CD4 v CD8B, CD4 v CD3
  
  Lo.GvHDcon.CD4.CD8B <- flowClust(GvHD.con, varNames = c("CD4", "CD8b"), B = 100, K = 5)
  Lo.GvHDcon.CD4.CD3 <- flowClust(GvHD.con, varNames = c("CD4", "CD3"), B = 100, K = 8)
  
  MRF.GvHDcon.CD4.CD8B <- mrf_gating(GvHD.con[,c("CD4", "CD8b")], temperature = 4)
  MRF.GvHDcon.CD4.CD3 <- mrf_gating(GvHD.con[,c("CD4", "CD3")], temperature = 4)
  
  png('paper/figures/GvHD_control.png', width = 1500, height = 1000)
  par(mfrow = c(2,3))
  par(pty = "s")
  plot(Lo.GvHDcon.CD4.CD8B, data = GvHD.con, xlab = "CD4", 
       ylab = "CD8b", las = 1, pch.outliers = 20,
       main = "(a) t mixture with Box-Cox", show.outliers = TRUE)
  plot(MRF.GvHDcon.CD4.CD8B, main = "(b) Markov random field probabilities")
  plot(MRF.GvHDcon.CD4.CD8B$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(c) Markov random field clusters",
       col = ifelse(MRF.GvHDcon.CD4.CD8B$groups == 0, "grey", MRF.GvHDcon.CD4.CD8B$groups))
 
  plot(Lo.GvHDcon.CD4.CD3, data = GvHD.con, xlab = "CD4", 
       ylab = "CD3", las = 1, pch.outliers = 20,
       main = "(d) t mixture with Box-Cox", show.outliers = TRUE)
  plot(MRF.GvHDcon.CD4.CD3, main = "(e) Markov random field probabilities")
  plot(MRF.GvHDcon.CD4.CD3$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(f) Markov random field clusters",
       col = ifelse(MRF.GvHDcon.CD4.CD3$groups == 0, "grey", MRF.GvHDcon.CD4.CD3$groups))
  par(mfrow = c(1, 1))
  dev.off()
  
  # Figure 4 - GvHD Data Positive Group CD4 v CD8B, CD4 v CD3
  
  Lo.GvHDpos.CD4.CD8B <- flowClust(GvHD.pos, varNames = c("CD4", "CD8b"), B = 100, K = 8)
  Lo.GvHDpos.CD4.CD3 <- flowClust(GvHD.pos, varNames = c("CD4", "CD3"), B = 100, K = 8)
  
  MRF.GvHDpos.CD4.CD8B <- mrf_gating(GvHD.pos[,c("CD4", "CD8b")], temperature = 4)
  MRF.GvHDpos.CD4.CD3 <- mrf_gating(GvHD.pos[,c("CD4", "CD3")], temperature = 4)
  
  png('paper/figures/GvHD_positive.png', width = 1500, height = 1000)
  par(mfrow = c(2,3))
  par(pty = "s")
  plot(Lo.GvHDpos.CD4.CD8B, data = GvHD.pos, xlab = "CD4", 
       ylab = "CD8b", las = 1, pch.outliers = 20, 
       main = "(a) t mixture with Box-Cox", show.outliers = TRUE)
  plot(MRF.GvHDpos.CD4.CD8B, main = "(b) Markov random field probabilities")
  plot(MRF.GvHDpos.CD4.CD8B$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(c) Markov random field clusters",
       col = ifelse(MRF.GvHDpos.CD4.CD8B$groups == 0, "grey", MRF.GvHDpos.CD4.CD8B$groups))
  
  plot(Lo.GvHDpos.CD4.CD3, data = GvHD.pos, xlab = "CD4", 
       ylab = "CD3", las = 1, pch.outliers = 20,
       main = "(d) t mixture with Box-Cox", show.outliers = TRUE)
  plot(MRF.GvHDpos.CD4.CD3, main = "(e) Markov random field probabilities")
  plot(MRF.GvHDpos.CD4.CD3$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(f) Markov random field clusters",
       col = ifelse(MRF.GvHDpos.CD4.CD3$groups == 0, "grey", MRF.GvHDpos.CD4.CD3$groups))
  
  par(mfrow = c(1, 1))
  dev.off()
  
  # Table 2 - GvHD Data
  
####--- End of Script ------------------------------------------------------####  