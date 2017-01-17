energy_system <- function(x, value = NULL, dimension) {
  
  if (is.null(value)) {
  
    energy <- matrix(0, nrow = dimension, ncol = dimension)
    for (i in 1:dimension) {
      for (j in 1:dimension) {
        pos <- j + (i - 1) * dimension
        energy[pos] <- energy(x, position = pos, value = x[j,i], dimension)
      }
    }

  } else {
    
    energy <- matrix(0, nrow = dimension, ncol = dimension)
    for (i in 1:dimension) {
      for (j in 1:dimension) {
        pos <- j + (i - 1) * dimension
        energy[pos] <- energy(x, position = pos, value = value)
      }
    }
    
  }
  
  return(energy)
}