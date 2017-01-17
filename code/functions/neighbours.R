neighbours <- function(x, position, value, dimension) {
  
  top.bottom.cases <- position %% dimension
  if (top.bottom.cases == 0) {
    neighbourhood <- position + c(-1, -dimension, dimension)
  } else if (top.bottom.cases == 1) {
    neighbourhood <- position + c(1, -dimension, dimension)
  } else {
    neighbourhood <- position + c(1, -1, dimension, -dimension)
  }
  
  neighbourhood <- neighbourhood[neighbourhood > 0]
  neighbourhood <- neighbourhood[neighbourhood < (dimension ^ 2 + 1)]
  
  neighbours <- sum(x[neighbourhood] == value)
  
  return(neighbours) 
}