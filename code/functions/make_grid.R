make_grid <- function(x, y = NULL, min = 0, max = 1023) {
  
  if (is.null(y)) {
    if (is.null(dim(x)) || dim(x)[2] < 2) {
      stop("If y is not supplied x needs to be a matrix of at least two columns")
    }
    if (dim(x)[2] > 2) {
      warning("x has more than 2 dimensions. Using the first and second columns")
    }
    y <- x[,2]
    x <- x[,1]
  }
  
  unique.recorded <- max - min + 1
  mat.grid <- matrix(0, nrow = unique.recorded, ncol = unique.recorded)
  
  index <- (x + 1) + y * unique.recorded
  index <- unique(index)
  
  mat.grid[index] <- 1
  class(mat.grid) <- "flow_grid"
  return(mat.grid)  
}