mrf_gating <- function(x, min = 0, max = 1023, temperature) {
  
  # Put the two variables of interest in an NxN matrix (field)  
  mat.grid <- make_grid(x, min = min, max = max)
  dimension <- nrow(mat.grid)
  
  # Remove all extremity values (keep record of their positions)
  mat.grid[ , dimension] <- mat.grid[dimension, ] <- 0 
  extreme.values <- unique(c(which(x[,1] == max), which(x[,2] == max)))
  removals <- rep(0, times = nrow(x))
  removals[extreme.values] <- 1
  
  # Markov Random Field Approach
  mat.grid <- grid_red(mat.grid, red.dim = 64, dimension)
  mrf.grid <- ising_model(mat.grid[[1]], temp = temperature)

  for (i in 1:4) {
    mrf.grid <- grid_inc(x = mrf.grid$state, nrow(mrf.grid$state))
    mrf.grid <- ising_model(mrf.grid, temp = temperature)
    print(paste0(dim(mrf.grid$state)[1], " x ", dim(mrf.grid$state)[1], " grid complete"))
  }  

  # Requirement to use connected components labelling
  mrf.grid.round <- round(mrf.grid$prob)

  # Connected Components Algorithm
  groups.grid <- SDMTools::ConnCompLabel(mrf.grid.round)

  # Identify the groups and return to data frame format
  groups <- unmake_grid(x = groups.grid, original = x, min = min, max = max)

  # Print Probability of active cell for each cell in the original matrix
  probs <- rep(0, times = dim(x)[1])
  for (j in 1:dim(x)[1]) {
    probs[j] <- mrf.grid$prob[x[j,1], x[j,2]]
  }
  
  # Output to return to the user
  output <- structure(list(x = x, groups = groups, probabilities = probs,
                           removals = removals, grid.probs = mrf.grid$prob), 
                      class = "mrf_gating")
  
  return(output)
}
