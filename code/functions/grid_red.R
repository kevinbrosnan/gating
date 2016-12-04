grid_red <- function(x, red.dim = 64) {
  
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  output <- vector("list", length = length(c(log2(red.dim):log2(no.rows))))
  
  if (no.rows != no.cols) {
    stop("This function currently only works on square matrices")
  }
  
  if (log2(no.rows) != round(log2(no.rows), digits = 0)) {
    stop("The number of rows and columns should be a power of 2 for fcs data")
  }
  
  if (!(red.dim %in% c(2^c(0:log2(no.rows))))) {
    stop("The reduced matrix should be a power of 2 for fcs data")
  }
  
  output[[length(output)]] <- x
  
  for (i in (length(output)-1):1) {
    
    tmp <- output[[i + 1]]
    evens <- seq(from = 1, to = nrow(tmp), by = 2)
    odds <- seq(from = 2, to = nrow(tmp), by = 2)
    tmp <- tmp[evens,] + tmp[odds,]
    tmp <- tmp[,evens] + tmp[,odds]
    tmp[tmp > 0] <- 1
    output[[i]] <- tmp 
    
  }
  
  
  return(output)
}