neighbours <- function(x, position, value) {
  
  no.rows <- nrow(x)
  
  if (position %% no.rows == 0) {
    neighbourhood <- position + c(-1, -no.rows, no.rows)
  } else if ((position - 1) %% no.rows == 0) {
    neighbourhood <- position + c(1, -no.rows, no.rows)
  } else {
    neighbourhood <- position + c(1, -1, no.rows, - no.rows)
  }
  
  neighbourhood <- neighbourhood[which(neighbourhood > 0)]
  neighbourhood <- neighbourhood[which(neighbourhood <= no.rows * no.rows)]
  
  neighbours <- sum(x[neighbourhood] == value)
  
  return(neighbours) 
}