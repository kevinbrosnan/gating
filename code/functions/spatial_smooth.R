spatial_smooth <- function(x, lambda = 0.5) {
  
  no.row <- nrow(x)
  no.col <- ncol(x)

  x <- x
  
  # Corners
  x[1,1] <- x[1,1] + lambda * (x[1,2] + x[2,1]) + lambda * (1/sqrt(2)) * x[2,2]
  x[1,1] <- x[1,1] / (1 + 2 * lambda + lambda / sqrt(2))
  
  x[no.row,1] <- x[no.row,1] + lambda * (x[no.row,2] + x[(no.row - 1),1]) + lambda * (1/sqrt(2)) * x[(no.row - 1),2]
  x[no.row,1] <- x[no.row,1] / (1 + 2 * lambda + lambda / sqrt(2))
  
  x[1,no.col] <- x[1,no.col] + lambda * (x[2,no.col] + x[1,(no.col - 1)]) + lambda * (1/sqrt(2)) * x[2,(no.col - 1)]
  x[1,no.col] <- x[no.row,1] / (1 + 2 * lambda + lambda / sqrt(2))
  
  x[no.row,no.col] <- x[no.row,no.col] + lambda * (x[no.row,(no.col - 1)] + x[(no.row - 1),no.col]) + lambda * (1/sqrt(2)) * x[(no.row - 1),(no.col - 1)]
  x[no.row,no.col] <- x[no.row,no.col] / (1 + 2 * lambda + lambda / sqrt(2))
  
  for (i in 2:(no.row - 1)) {
    x[i, 1] <- (x[i, 1] + lambda * (x[(i-1),1] + x[(i+1),1] + x[i,2]) + lambda * (1/sqrt(2)) * (x[(i-1), 2] + x[(i+1),2])) / (1 + 3 * lambda + sqrt(2) * lambda)
    x[i,no.col] <- (x[i,no.col] + lambda * (x[(i-1),no.col] + x[(i+1),no.col] + x[i,(no.col - 1)]) + lambda * (1/sqrt(2)) * (x[(i-1), (no.col - 1)] + x[(i+1),(no.col - 1)])) / (1 + 3 * lambda + sqrt(2) * lambda)  
  }
  
  for (i in 2:(no.col - 1)) {
    x[1, i] <- (x[1, i] + lambda * (x[1, (i - 1)] + x[1, (i + 1)] + x[2, i]) + lambda * (1/sqrt(2)) * (x[2, (i - 1)] + x[2, (i + 1)])) / (1 + 3 * lambda + sqrt(2) * lambda)
    x[no.row,i] <- (x[no.row,i] + lambda * (x[no.row,(i-1)] + x[no.row,(i+1)] + x[(no.row - 1),i]) + lambda * (1/sqrt(2)) * (x[(no.row - 1),(i-1)] + x[(no.row - 1),(i+1)])) / (1 + 3 * lambda + sqrt(2) * lambda)
  }
  
  for (i in 2:(no.row - 1)) {
    for (j in 2:(no.col - 1)) {
      x[i,j] <- x[i,j] + lambda * (x[(i - 1), j] + x[(i + 1), j] + x[i, (j - 1)] + x[i, (j + 1)]) + lambda * (1/(sqrt(2))) * (x[(i - 1), (j - 1)] + x[(i - 1), (j + 1)] + x[(i + 1), (j + 1)] + x[(i + 1), (j - 1)])
      x[i,j] <- x[i,j]/(1 + 4 * lambda + 2 * sqrt(2) * lambda)
    }
  }  
  
  for (i in (no.row - 1):2) {
    for (j in (no.col - 1):2) {
      x[i,j] <- x[i,j] + lambda * (x[(i - 1), j] + x[(i + 1), j] + x[i, (j - 1)] + x[i, (j + 1)]) + lambda * (1/(sqrt(2))) * (x[(i - 1), (j - 1)] + x[(i - 1), (j + 1)] + x[(i + 1), (j + 1)] + x[(i + 1), (j - 1)])
      x[i,j] <- x[i,j]/(1 + 4 * lambda + 2 * sqrt(2) * lambda)
    }
  }
  
  return(x)
}