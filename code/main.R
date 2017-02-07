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

  # Required Libraries (Uncomment next three lines for installation)
  # install.packages(c('stargazer', 'SDMTools', 'tiff'))
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("flowClust") 
  library(SDMTools)
  library(stargazer)
  library(flowClust)
  library(tiff)

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
  
  # Figure 2 - Rituximab Data FSC-Height v SSC-Height & 7-AAD v Anti-BrdU
  Lo.rit <- flowClust(rituximab, varNames = c("FSC.H", "SSC.H"), K = 1, 
                      z.cutoff = 0.5)
  Lo.rit.gate <- rituximab[rituximab %in% Lo.rit, c("FL3.H", "FL1.H")]
  Lo.rit.stagetwo <- flowClust(Lo.rit.gate, K = 3)
  
  MRF.rit <- mrf_gating(rituximab[, c("FSC.H", "SSC.H")], temperature = 4)
  MRF.rit.gate <- rituximab[which(MRF.rit$groups == 1), c("FL3.H", "FL1.H")]
  MRF.rit.stagetwo <- mrf_gating(MRF.rit.gate, temperature = 4)
  
  
  tiff('paper/figures/rituximab.tiff', height = 8, width = 12, res = 600,
         units = 'in')
  par(mfrow = c(2,3))
  par(pty = "s")
  plot(Lo.rit, data = rituximab, xlab = "FSC-Height", 
       ylab = "SSC-Height", las = 1, 
       main = "(a) t mixture with Box-Cox", show.outliers = TRUE,
       pch.outliers = 20, cex = 1.5, cex.outliers = 1)
  plot(MRF.rit, xlab = "FSC-Height", ylab = "SSC-Height",
  		main = "(b) Markov random field probabilities")
  plot(MRF.rit$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(c) Markov random field clusters", las = 1,
       xlab = "FSC-Height", ylab = "SSC-Height",
       col = ifelse(MRF.rit$groups == 0, "grey", MRF.rit$groups))
  
  plot(Lo.rit.stagetwo, data = Lo.rit.gate, xlab = "7 AAD", 
       ylab = "Anti-BrdU FITC", las = 1, 
	 xlim = c(0, 1023), ylim = c(0, 1023),
       main = "(d) t mixture with Box-Cox", show.outliers = TRUE,
       pch.outliers = 20, cex = 1.5, cex.outliers = 1)
  plot(MRF.rit.stagetwo, xlab = "7 AAD", ylab = "Anti-BrdU FITC",
       main = "(e) Markov random field probabilities")
  plot(MRF.rit.stagetwo$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(f) Markov random fields cluster", las = 1,
       xlab = "7 AAD", ylab = "Anti-BrdU FITC",
       col = ifelse(MRF.rit.stagetwo$groups == 0, "grey", 
                    MRF.rit.stagetwo$groups))
  par(mfrow = c(1, 1))
  dev.off()
  
  # Table 1 - Rituximab Data
  
    # t-mixtures - (stage two isn't complete due to cluster number changing on runs)
    Lo.rit.table <- list(cluster.count = sum(!Lo.rit@flagOutliers, na.rm = TRUE), 
                         cluster.percent = round((sum(!Lo.rit@flagOutliers, na.rm = TRUE)/1545) * 100, digits = 2),
                         undefined.count = sum(Lo.rit@flagOutliers, na.rm = TRUE),
                         undefined.percent = round((sum(Lo.rit@flagOutliers, na.rm = TRUE)/1545) * 100, digits = 2),
                         doublets.count = length(which(is.na(Lo.rit@z))),
                         doublets.percent = round((length(which(is.na(Lo.rit@z)))/1545) * 100, digits = 2)
    )
    Lo.rit.table.stagetwo <- list(cluster.count = table(Lo.rit.stagetwo@label[which(!Lo.rit.stagetwo@flagOutliers)]), 
                         cluster.percent = round((table(Lo.rit.stagetwo@label[which(!Lo.rit.stagetwo@flagOutliers)])/1545) * 100, digits = 2),
                         undefined.count = sum(Lo.rit.stagetwo@flagOutliers, na.rm = TRUE),
                         undefined.percent = round((sum(Lo.rit.stagetwo@flagOutliers, na.rm = TRUE)/1545) * 100, digits = 2)
    )
    
    # Markov Random Field

    MRF.rit.table <- list(cluster.count = sum(MRF.rit$groups == 1, na.rm = TRUE), 
                         cluster.percent = round((sum(MRF.rit$groups == 1, na.rm = TRUE)/1545) * 100, digits = 2),
                         undefined.count = sum(MRF.rit$groups == 0, na.rm = TRUE),
                         undefined.percent = round((sum(MRF.rit$groups == 0, na.rm = TRUE)/1545) * 100, digits = 2),
                         doublets.count = sum(MRF.rit$removals, na.rm = TRUE),
                         doublets.percent = round((sum(MRF.rit$removals, na.rm = TRUE)/1545) * 100, digits = 2),
                         cluster.level = c(length(which(MRF.rit$probs < 0.5)),
                                           length(which(MRF.rit$probs < 0.75)) - length(which(MRF.rit$probs < 0.5)),
                                           length(which(MRF.rit$probs >= 0.75)))
    )
    
    MRF.rit.stagetwo.table <- list(cluster.count = as.vector(table(MRF.rit.stagetwo$groups))[-1], 
                          cluster.percent = as.vector(round((table(MRF.rit.stagetwo$groups)[-1]/1545) * 100, digits = 2)),
                          undefined.count = sum(MRF.rit.stagetwo$groups == 0, na.rm = TRUE),
                          undefined.percent = round((sum(MRF.rit.stagetwo$groups == 0, na.rm = TRUE)/1545) * 100, digits = 2),
                          doublets.count = sum(MRF.rit.stagetwo$removals, na.rm = TRUE),
                          doublets.percent = round((sum(MRF.rit.stagetwo$removals, na.rm = TRUE)/1545) * 100, digits = 2),
                          cluster.level = list(low = as.vector(by(MRF.rit.stagetwo$probs,
                                                                  as.factor(MRF.rit.stagetwo$groups),
                                                                  function(x) length(which(x < 0.5))))[-1],
                                               medium = as.vector(by(MRF.rit.stagetwo$probs,
                                                                     as.factor(MRF.rit.stagetwo$groups),
                                                                     function(x) length(which(x < 0.75)) - length(which(x < 0.5))))[-1],
                                               high = as.vector(by(MRF.rit.stagetwo$probs,
                                                                   as.factor(MRF.rit.stagetwo$groups),
                                                                   function(x) length(which(x >= 0.75))))[-1]
                                            )
    )
    
  
  # Figure 3 - GvHD Data Control Group CD4 v CD8B, CD4 v CD3
  
  Lo.GvHDcon.CD8B <- flowClust(GvHD.con, varNames = c("CD4", "CD8b"), K = 4)
  Lo.GvHDcon.CD3 <- flowClust(GvHD.con, varNames = c("CD4", "CD3"), K = 5)
  
  MRF.GvHDcon.CD8B <- mrf_gating(GvHD.con[,c("CD4", "CD8b")], temperature = 4, cluster.prob = 0.75)
  MRF.GvHDcon.CD3 <- mrf_gating(GvHD.con[,c("CD4", "CD3")], temperature = 4, cluster.prob = 0.75)
  

  tiff('paper/figures/GvHD_con.tiff', height = 8, width = 12, res = 600,
	 units = 'in')
  par(mfrow = c(2,3))
  par(pty = "s")
  plot(Lo.GvHDcon.CD8B, data = GvHD.con, xlab = "CD4", 
       ylab = "CD8b", las = 1, pch.outliers = 20, 
	 xlim = c(0, 1023), ylim = c(0, 1023),
	 cex = 1.5, cex.outliers = 1,
       main = "(a) t mixture with Box-Cox", show.outliers = TRUE)
  plot(MRF.GvHDcon.CD8B, main = "(b) Markov random field probabilities")
  plot(MRF.GvHDcon.CD8B$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(c) Markov random field clusters", las = 1,
       col = ifelse(MRF.GvHDcon.CD8B$groups == 0, "grey", 
                    MRF.GvHDcon.CD8B$groups))
 
  plot(Lo.GvHDcon.CD3, data = GvHD.con, xlab = "CD4", 
       ylab = "CD3", las = 1, pch.outliers = 20, 
	 xlim = c(0, 1023), ylim = c(0, 1023),
	 cex = 1.5, cex.outliers = 1,
       main = "(d) t mixture with Box-Cox", show.outliers = TRUE)
  plot(MRF.GvHDcon.CD3, main = "(e) Markov random field probabilities")
  plot(MRF.GvHDcon.CD3$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(f) Markov random field clusters", las = 1,
       col = ifelse(MRF.GvHDcon.CD3$groups == 0, "grey", 
                    MRF.GvHDcon.CD3$groups))
  par(mfrow = c(1, 1))
  dev.off()
  
  # Figure 4 - GvHD Data Positive Group CD4 v CD8B, CD4 v CD3
  
  Lo.GvHDpos.CD8B <- flowClust(GvHD.pos, varNames = c("CD4", "CD8b"), K = 6)
  Lo.GvHDpos.CD3 <- flowClust(GvHD.pos, varNames = c("CD4", "CD3"), K = 5)
  
  MRF.GvHDpos.CD8B <- mrf_gating(GvHD.pos[,c("CD4", "CD8b")], temperature = 4, cluster.prob = 0.75)
  MRF.GvHDpos.CD3 <- mrf_gating(GvHD.pos[,c("CD4", "CD3")], temperature = 4, cluster.prob = 0.75)
  
  tiff('paper/figures/GvHD_pos.tiff', height = 8, width = 12, res = 600,
	 units = 'in')
  par(mfrow = c(2,3))
  par(pty = "s")
  plot(Lo.GvHDpos.CD8B, data = GvHD.pos, xlab = "CD4", 
       ylab = "CD8b", las = 1, pch.outliers = 20, 
       xlim = c(0, 1023), ylim = c(0, 1023),
	 cex = 1.5, cex.outliers = 1,
       main = "(a) t mixture with Box-Cox", show.outliers = TRUE)
  plot(MRF.GvHDpos.CD8B, main = "(b) Markov random field probabilities")
  plot(MRF.GvHDpos.CD8B$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(c) Markov random field clusters", las = 1,
       col = ifelse(MRF.GvHDpos.CD8B$groups == 0, "grey", 
                    MRF.GvHDpos.CD8B$groups))
  
  plot(Lo.GvHDpos.CD3, data = GvHD.pos, xlab = "CD4", 
       ylab = "CD3", las = 1, pch.outliers = 20,
       xlim = c(0, 1023), ylim = c(0, 1023),
	 cex = 1.5, cex.outliers = 1,
       main = "(d) t mixture with Box-Cox", show.outliers = TRUE)
  plot(MRF.GvHDpos.CD3, main = "(e) Markov random field probabilities")
  plot(MRF.GvHDpos.CD3$x, xlim = c(0, 1023), ylim = c(0,1023),
       main = "(f) Markov random field clusters", las = 1,
       col = ifelse(MRF.GvHDpos.CD3$groups == 0, "grey", 
                    MRF.GvHDpos.CD3$groups))
  par(mfrow = c(1, 1))
  dev.off()
  
  # Table 2 - GvHD Data
  
####--- End of Script ------------------------------------------------------####  