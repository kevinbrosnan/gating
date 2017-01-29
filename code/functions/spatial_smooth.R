spatial_smooth <- function(x) {
  
  no.row <- nrow(x)
  no.col <- ncol(x)

  x <- x
  
  for (i in 2:(no.row - 1)) {
    for (j in 2:(no.col - 1)) {
      x[i,j] <- x[i,j] + 0.5 * (x[(i - 1), j] + x[(i + 1), j] + x[i, (j - 1)] + x[i, (j + 1)]) + (1/(2*sqrt(2))) * (x[(i - 1), (j - 1)] + x[(i - 1), (j + 1)] + x[(i + 1), (j + 1)] + x[(i + 1), (j - 1)])
      x[i,j] <- x[i,j]/(3 + sqrt(2))
    }
  }  
  
  
  return(x)
}