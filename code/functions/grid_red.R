grid_red <- function(x, red.dim = 64, dimension) {
  
  output <- vector("list", length = length(c(log2(red.dim):log2(dimension))))
  
  if (dimension != nrow(x)) {
    stop("Dimension parameter passed is incorrect")
  }
  
  if (log2(dimension) != round(log2(dimension), digits = 0)) {
    stop("The number of rows and columns should be a power of 2 for fcs data")
  }
  
  if (!(red.dim %in% c(2^c(0:log2(dimension))))) {
    stop("The reduced matrix should be a power of 2 for fcs data")
  }
  
  output[[length(output)]] <- x
  
  for (i in (length(output) - 1):1) {
    tmp <- output[[i + 1]]
    tmp.dim <- nrow(tmp)
    evens <- seq(from = 1, to = tmp.dim, by = 2)
    odds <- seq(from = 2, to = tmp.dim, by = 2)
    tmp <- tmp[evens,] + tmp[odds,]
    tmp <- tmp[,evens] + tmp[,odds]
    tmp[tmp > 0] <- 1
    output[[i]] <- tmp 
  }
  
  return(output)
}