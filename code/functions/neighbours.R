neighbours <- function(x, position, value, dimension) {
  
  top.bottom.cases <- position %% dimension
  
  neighbourhood <- position + c(1, -1, dimension, -dimension)
  if (top.bottom.cases == 0) {
    neighbourhood <- neighbourhood[-1]
  } else if (top.bottom.cases == 1) {
    neighbourhood <- neighbourhood[-2]
  }
  
  neighbourhood <- neighbourhood[neighbourhood > 0]
  neighbourhood <- neighbourhood[neighbourhood < (dimension ^ 2 + 1)]
  
  neighbours <- sum(x[neighbourhood] == value)
  
  return(neighbours) 
}