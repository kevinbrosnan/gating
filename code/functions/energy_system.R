energy_system <- function(x, value = NULL) {
  
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  if (is.null(value)) {
  
    energy <- 0
    for (i in 1:no.cols) {
      for (j in 1:no.rows) {
        pos <- j + (i - 1) * no.cols
        energy <- energy + energy(x, position = pos, value = x[j,i])
      }
    }

  } else {
    
    energy <- 0
    for (i in 1:no.cols) {
      for (j in 1:no.rows) {
        pos <- j + (i - 1) * no.cols
        energy <- energy + energy(x, position = pos, value = value)
      }
    }
    
  }
  
  return(energy)
}