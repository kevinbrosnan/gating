grid_inc <- function(x) {
  
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  y <- matrix(0, nrow = 2 * no.rows, ncol = 2 * no.cols)
  
  odd <- seq(1, nrow(y), by = 2)
  even <- seq(2, nrow(y), by = 2)
  order <- rep(0, times = nrow(y))
  
  order[odd] <- seq(1, nrow(y)/2, by = 1)
  order[even] <- seq((nrow(y)/2 + 1), nrow(y), by = 1)
  
  for (i in 1:length(odd)) {
    y[,c(odd[i], even[i])] <- x[, i]
  }
  
  y <- y[order,]
  
  return(y)
}