neigh_system <- function(x, value = 1) {
  
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  neigh <- matrix(0, nrow = no.rows, ncol = no.cols)
  for (i in 1:no.cols) {
    for (j in 1:no.rows) {
      pos <- j + (i - 1) * no.cols
      neigh[pos] <- neighbours(x, pos, value)
    }
  }
    
  return(neigh)
}