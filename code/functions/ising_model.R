ising_model <- function(x, temp = 4) {
  
  # Dimensions of input matrix
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  # Simulated Annealing Stopping Criterion
  SA.stop <- FALSE
  SA.update <- FALSE
  SA.updates <- 0
  
  # Defined Constant
  ln.alpha <- log(runif(1))
  energy.system <- sum(energy_system(x))
  
  # Holding matrices
  current.state <- x
  prob <- matrix(1, nrow = no.rows, ncol = no.cols)
  
  while (SA.stop != TRUE) {
    
    temp.cur <- temp * (0.95 ^ SA.updates)
    SA.updates <- SA.updates + 1
    SA.update <- FALSE
    prob.previous <- prob
    SA.update.runs <- 0
    
    while (SA.update != TRUE) {
      
      SA.update.runs <- SA.update.runs + 1
      # Pick a random point and get its linear index
      row.index <- ceiling(no.rows * runif(1))
      col.index <- ceiling(no.cols * runif(1))
      lin.index <- (col.index - 1) * no.rows + row.index 
      
      # Energy Change
      energy.cur <- energy(current.state, position = lin.index, 
                           value = current.state[[lin.index]])
      energy.swap <- energy(current.state, position = lin.index, 
                            value = 1 - current.state[[lin.index]])
      energy.change <- energy.swap - energy.cur
      
      # Metropolis Criterion                     
      if (energy.change <= 0) {
        current.state[lin.index] <- 1 - current.state[lin.index]
        
        # Update the temperature?
        if (abs(energy.change/energy.system) < 1e-5) {
          SA.update <- TRUE
        }
      } else if (ln.alpha <= -energy.change/temp.cur) {
        current.state[lin.index] <- 1 - current.state[lin.index]
        
        # Update the temperature?
        if (abs(energy.change/energy.system) < 1e-5) {
          SA.update <- TRUE
        }
      }
    }
    
    neigh.ones <- neigh_system(current.state, value = 1)
    neigh.zeros <- neigh_system(current.state, value = 0)

    prob.final.ones <- exp((1/temp) * neigh.ones)
    prob.final.zeros <- exp((1/temp) * neigh.zeros)
    prob <- prob.final.ones / (prob.final.ones + prob.final.zeros)
    
    # Are we finished annealing?
    if ((abs(prob - prob.previous) / length(prob)) < 1e-4 || SA.updates >= 43) {
      SA.stop <- TRUE
    }
  }
  
  output <- structure(list(prob = prob, temperature = temp, 
                           state = current.state, SA.update.runs))
  
  return(output)
}