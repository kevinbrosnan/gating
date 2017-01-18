energy_system <- function(x, value = NULL, dimension) {
  
  if (is.null(value)) {
  
    energy <- matrix(0, nrow = dimension, ncol = dimension)
    for (i in 1:(dimension^2)) {
      energy[i] <- energy(x, position = i, value = x[i], dimension)
    }

  } else {
    
    energy <- matrix(0, nrow = dimension, ncol = dimension)
    for (i in 1:(dimension^2)) {
      energy[i] <- energy(x, position = i, value = value)
    }
    
  }
  
  return(energy)
}