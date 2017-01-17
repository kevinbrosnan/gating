energy <- function(x, position, value, dimension) {
  
  if (position %% dimension == 0) {
    neighbourhood <- position + c(-1, -dimension, dimension)
  } else if ((position - 1) %% dimension == 0) {
    neighbourhood <- position + c(1, -dimension, dimension)
  } else {
    neighbourhood <- position + c(1, -1, dimension, -dimension)
  }
  
  neighbourhood <- neighbourhood[neighbourhood > 0]
  neighbourhood <- neighbourhood[neighbourhood <= dimension ^ 2]
  
  energy <- 2 * sum(x[neighbourhood] == value) - length(neighbourhood)
  
  return(energy) 
}