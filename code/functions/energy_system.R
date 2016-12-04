energy_system <- function(x, value = NULL) {
  
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  if (is.null(value)) {
  
    energy <- matrix(0, nrow = no.rows, ncol = no.cols)
    for (i in 1:no.cols) {
      for (j in 1:no.rows) {
        pos <- j + (i - 1) * no.cols
        energy[pos] <- energy(x, position = pos, value = x[j,i])
      }
    }

  } else {
    
    energy <- matrix(0, nrow = no.rows, ncol = no.cols)
    for (i in 1:no.cols) {
      for (j in 1:no.rows) {
        pos <- j + (i - 1) * no.cols
        energy[pos] <- energy(x, position = pos, value = value)
      }
    }
    
  }
  
  return(energy)
}