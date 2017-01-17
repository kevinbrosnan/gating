neigh_system <- function(x, value = 1, dimension) {
  
  neigh <- matrix(0, nrow = dimension, ncol = dimension)
  for (i in 1:dimension) {
    for (j in 1:dimension) {
      pos <- j + (i - 1) * dimension
      neigh[pos] <- neighbours(x, pos, value, dimension)
    }
  }
    
  return(neigh)
}