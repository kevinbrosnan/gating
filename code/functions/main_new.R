##---- Markov Random Fields Ising model with HMRF & Simulated Annealing ----##


## Data Required

  # GitHub Access
  gh.proj <- "https://raw.githubusercontent.com/significantstats/gating/master/"
  gh.data <- paste0(gh.proj, "data/")

  # Data Loading
  GvHD.pos <- read.csv(paste0(gh.data, "GvHD/GvHD_positive.csv"))
  GvHD.con <- read.csv(paste0(gh.data, "GvHD/GvHD_control.csv"))
  rituximab <- read.csv(paste0(gh.data, "rituximab/rituximab.csv"))

## Functions Required
  
  # GitHub Access
  gh.code <- paste0(gh.proj, "code/functions/")
  
  # Function Loading
  source(paste0(gh.code, "make_grid.R"))
  source(paste0(gh.code, "grid_red.R"))
  source(paste0(gh.code, "ising_model.R"))
  source(paste0(gh.code, "energy.R"))
  source(paste0(gh.code, "energy_system.R"))
  
##---- Methodology Applied to GvHD data ----##

adc.level <- ceiling(max(log2(max(GvHD.con)), log2(max(GvHD.pos))))

CD4.v.CD8B <- make_grid(x = GvHD.con[, c(1,2)], min = 0, max = (2^adc.level - 1))
CD4.v.CD8B.levels <- grid_red(CD4.v.CD8B, red.dim = 64)
CD4.v.CD8B.prob.map <- vector('list', length = length(CD4.v.CD8B.levels))

for (i in 1:length(CD4.v.CD8B.levels)) {
  temperature <- 4
  CD4.v.CD8B.prob.map[[i]] <- ising_model(CD4.v.CD8B.levels[[i]], temp = temperature)
}
