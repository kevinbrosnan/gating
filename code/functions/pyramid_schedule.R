pyramid_schedule <- function(x) {
  
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  if (no.rows != no.cols || no.rows %% 2 != 0) {
    stop('Number of rows should equal number of columns and should be divisible by 2')
  }
  
  res.levels <- 2^c(6:log2(no.rows))
  resolutions <- vector('list', length(res.levels))
  
  for (i in 1:(length(res.levels)-1)) {
    tmp <- matrix(0, nrow = res.levels[i], ncol = res.levels[i])
    reduction <- no.rows / res.levels[i]
    for (j in 1:res.levels[i]) {
      for (k in 1:res.levels[i]) {
        tmp[j,k] <- sum(x[((j-1) * reduction + 1):(j * reduction),((k-1) * reduction + 1):(k * reduction)])
      }
    }
    tmp[tmp > 0] <- 1
    resolutions[[i]] <- tmp
  }
  resolutions[[length(res.levels)]] <- x
  
  return(resolutions)
}