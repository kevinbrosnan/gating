neighbours <- function(x, position, value, dimension) {
  
  if (position %% dimension == 0) {
    neighbourhood <- position + c(-1, -dimension, dimension)
  } else if ((position - 1) %% dimension == 0) {
    neighbourhood <- position + c(1, -dimension, dimension)
  } else {
    neighbourhood <- position + c(1, -1, dimension, -dimension)
  }
  
  neighbourhood <- neighbourhood[which(neighbourhood > 0)]
  neighbourhood <- neighbourhood[which(neighbourhood <= dimension ^ 2)]
  
  neighbours <- sum(x[neighbourhood] == value)
  
  return(neighbours) 
}