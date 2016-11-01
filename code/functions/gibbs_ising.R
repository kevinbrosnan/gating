gibbs_ising <- function(x) {
  
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  x[x == 0] <- -1
  
  iterations <- 0
  
  for (jj in 1:30) {
    for (i in 1:10000) {
      row.index <- ceiling(no.rows * runif(1))
      col.index <- ceiling(no.cols * runif(1))
      cur.x <- x[row.index, col.index]
      cur.pos <- (row.index + 1) + no.rows * col.index
      neighbourhood <- cur.pos + c(-1, 1, -no.rows, +no.rows)
      if (row.index == 1) {
        if (col.index == 1) {
          neighbourhood <- neighbourhood[-c(1, 2)]
        } else if (col.index == no.cols) {
          neighbourhood <- neighbourhood[-c(2, 4)]
        } else {
          neighbourhood <- neighbourhood[-c(2)]
        }
      } else if (row.index == no.rows) {
        if (col.index == 1) {
          neighbourhood <- neighbourhood[-c(1, 3)]
        } else if (col.index == no.cols) {
          neighbourhood <- neighbourhood[-c(3, 4)]
        } else {
          neighbourhood <- neighbourhood[-c(3)]
        }
      } else if (col.index == 1) {
        neighbourhood <- neighbourhood[-c(1)]
      } else if (col.index == no.cols) {
        neighbourhood <- neighbourhood[-c(4)]
      }
      
      potential <- sum(x[neighbourhood])
      if (runif(1) < exp(-beta * potential) / (exp(-beta * potential) + exp(beta * potential))) {
        x[row.index, col.index] <- -1
      } else {
        x[row.index, col.index] <- 1
      }
      iterations <- iterations + 1
    }
  }
}