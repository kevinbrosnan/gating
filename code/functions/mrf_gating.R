mrf_gating <- function(x, min = 0, max = 1023, temperature, cluster.prob = 0.25) {
  
  # Put the two variables of interest in an NxN matrix (field)  
  mat.grid <- make_grid(x, min = min, max = max)
  dimension <- nrow(mat.grid)
  
  # Remove all extremity values (keep record of their positions)
  mat.grid[ , dimension] <- mat.grid[dimension, ] <- 0 
  extreme.values <- unique(c(which(x[,1] == max), which(x[,2] == max)))
  removals <- rep(0, times = nrow(x))
  removals[extreme.values] <- 1
  
  # Markov Random Field Approach
    # 64 x 64 Grid
    mat.grid <- grid_red(mat.grid, red.dim = 64, dimension)
    mrf.grid <- ising_model(prob = mat.grid[[1]], temp = temperature)
    
    # 128 x 128 Grid
    mrf.grid.prob <- grid_inc(x = spatial_smooth(mrf.grid$prob), nrow(mrf.grid$prob))
    mrf.grid <- ising_model(prob = spatial_smooth(mrf.grid.prob), temp = temperature)

    # 256 x 256 Grid
    mrf.grid <- grid_inc(x = spatial_smooth(mrf.grid$prob), nrow(mrf.grid$prob))
    mrf.grid <- ising_model(spatial_smooth(mrf.grid), temp = temperature)

  # Requirement to use connected components labelling
  mrf.grid.round <- ifelse(spatial_smooth(mrf.grid$prob) > cluster.prob, 1, 0)

  # Connected Components Algorithm
  groups.grid <- SDMTools::ConnCompLabel(mrf.grid.round)
  
  # Scale up the grids to the correct dimension of the data
  groups.grid <- grid_inc(grid_inc(groups.grid, 256), 512)
    
  # Identify the groups and return to data frame format
  groups <- unmake_grid(x = groups.grid, original = x, min = min, max = max)
  groups[which(removals == 1)] <- NA
  
  # Output to return to the user
  output <- structure(list(x = x, groups = groups, removals = removals,
                           grid.probs = spatial_smooth(mrf.grid$prob)), 
                      class = "mrf_gating")
  
  return(output)
}

plot.mrf_gating <- function(object, xlab = NULL, ylab = NULL, main = NULL) {
  
  if (!inherits(object, "mrf_gating"))
    stop('object not of class \"mrf_gating\"')
  
  if (is.null(xlab)) {
    xlab <- names(object$x)[[1]]
  }
  
  if (is.null(ylab)) {
    ylab <- names(object$x)[[2]]
  }
  
  if (is.null(main)) {
    main <- ""
  }
  
  no.row <- nrow(object$grid.probs)
  no.col <- ncol(object$grid.probs)
  axis.pos <- c(0, 200, 400, 600, 800, 1000) * (255/1023)
  
  # Initilise Plot
  par(pty = "s")
  plot(0, type = "n", las = 1, xlab = xlab, ylab = ylab,
       main = main, xlim = c(-0.5, no.row), ylim = c(-0.5, no.col),
       xaxt = "n", yaxt = "n", bty = "o")
  axis(side = 1, at = axis.pos, labels = c(0, 200, 400, 600, 800, 1000))
  axis(side = 2, at = axis.pos, labels = c(0, 200, 400, 600, 800, 1000), las = 1)
  
  # Add probability map layer
  colours.scale <- c("yellow", "orange", "red")
  
  probs.x <- rep(1:no.row, times = no.col)
  probs.y <- sort(rep(1:no.col, times = no.row))
  probs.val <- as.vector(object$grid.probs)
  probs.col <- ifelse(probs.val < 0.5, colours.scale[1],
                      ifelse(probs.val < 0.75, colours.scale[2], colours.scale[3]))
  probs.plot <- ifelse(probs.val < 0.25, FALSE, TRUE)
  
  for (i in 1:length(probs.x)) {
    if (probs.plot[i] == TRUE) {
      polygon(x = c((probs.x[i] - 1.5), (probs.x[i] - 1.5), 
                    (probs.x[i] - 0.5), (probs.x[i] - 0.5)), 
              y = c((probs.y[i] - 1.5), (probs.y[i] - 0.5), 
                    (probs.y[i] - 0.5), (probs.y[i] - 1.5)),
              col = probs.col[i], border = NA)
    }
  }
  
  # Add points coloured by cluster
  points.characters <- ifelse(is.na(object$groups), 8, 1)
  points.colours <- ifelse(is.na(object$groups), "grey",
                         ifelse(object$groups == 0, "grey", "black"))
  points(object$x * (255/1023), col = points.colours, pch = points.characters, cex = 0.2)
  
  invisible()
}
