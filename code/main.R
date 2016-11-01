#------------------------------------------------------------------------------#
#                  Automated Gating of Flow Cytometry Data via                 #
#                   Adaptive Markov Random Fields Clustering                   #
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

  # Required Functions
  functions.path <- "functions"
  sub.dir <- dir(path = functions.path, recursive = TRUE, full.names = TRUE)
  invisible(sapply(sub.dir, source))
  rm(functions.path, sub.dir)
  
  # Required Libraries
  #install.packages(c('tools', 'stargazer', 'SDMTools'))
  library(tools)
  library(SDMTools)
  library(stargazer)

####--- Rituximab Data Example ---------------------------------------------####

  # Read in Rituximab Data - online resource
  github_account <- "https://raw.githubusercontent.com/significantstats/"
  file_location <- "automatedgating/master/rituximab.csv"
  rituximab <- read.csv(file = paste0(github_account, file_location))
  
  # Exploratory Analysis
  explore_flow(rituximab)
  stargazer(rituximab, summary = TRUE, digits = 2, label = "tab:ritdatasum", 
            median = TRUE, iqr = TRUE, mean.sd = FALSE,
            title = "Summary Statistics for Rituximab Data", align = TRUE)
  
  # Gating FSC v SSC
  gate.initial <- mrf_gating(rituximab[,c('FSC.H', 'SSC.H')])
  plot(gate.initial)
  
  # Gating 7 AAD v Anti-BrdU FITC
  gate.analysis <- mrf_gating(rituximab[which(gate.initial$groups == 1),c('FL1.H', 'FL3.H')])
  plot(gate.analysis)

####--- GvHD Data Example --------------------------------------------------####

  # Read in Rituximab Data - online resource
  github_account <- "https://raw.githubusercontent.com/significantstats/"
  file_location <- "automatedgating/master/"
  gvhd.con <- read.csv(file = paste0(github_account, file_location, 'GvHD_control.csv'))
  gvhd.pos <- read.csv(file = paste0(github_account, file_location, 'GvHD_positive.csv'))
  
  # Exploratory Analysis
  explore_flow(gvhd.con)
  stargazer(gvhd.con, summary = TRUE, digits = 2, label = "tab:gvhd.con", 
            median = TRUE, iqr = TRUE, mean.sd = FALSE,
            title = "Summary Statistics for GvHD Control Data", align = TRUE)
  
  explore_flow(gvhd.pos)
  stargazer(gvhd.pos, summary = TRUE, digits = 2, label = "tab:gvhd.pos", 
            median = TRUE, iqr = TRUE, mean.sd = FALSE,
            title = "Summary Statistics for GvHD Positive Data", align = TRUE)
  
  # Gating FSC v SSC
  gate.initial <- mrf_gating(rituximab[,c('FSC.H', 'SSC.H')])
  plot(gate.initial)
  
  # Gating 7 AAD v Anti-BrdU FITC
  gate.analysis <- mrf_gating(gate.initial$population[,c('FL1.H', 'FL3.H')])
  plot(gate.analysis)


