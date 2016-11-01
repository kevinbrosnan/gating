ising_denosing <- function(x, J = 0.5) {
  
  x[which(x == 0)] <- -1
  
  no.rows <- nrow(x)
  no.cols <- ncol(x)
  
  theta <- matrix(1, nrow = no.rows, ncol = no.cols)
  mf <- matrix(0, nrow = no.rows, ncol = no.cols)
  misfit <- vector(mode = "numeric")
  iter <- 0
  prob.array <- vector('list', 100)
  while (iter < 1000001) {
    # Pick a random point and then using its linear index
    row.index <- ceiling(no.rows * runif(1))
    col.index <- ceiling(no.cols * runif(1))
    lin.index <- (col.index - 1) * no.rows + row.index 
    
    # These won't be needed once I apply changes below!
    thetap <- -theta[lin.index]
    lik.ratio <- exp(x[lin.index] * (thetap - theta[lin.index]))
    
    # Neighbourhood System for 2-D Ising
    neighbourhood <- lin.index + c(-no.rows, -1, +1, no.rows)
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
    
    # I just need to add in the energy function here!
    disagree <- sum(theta[neighbourhood] != theta[lin.index])
    disagreep <- sum(theta[neighbourhood] != thetap)
    DelLogPr <- 2 * J * (disagree - disagreep)
    
    alpha <- exp(DelLogPr) * lik.ratio
    if (runif(1) < alpha) {
      theta[lin.index] <- thetap
    }
    
    # Leave this as is - It makes the output more useable!
    iter <- iter + 1
    if (iter %% 10000 == 0) {
      mf <- mf + theta
      prob <- exp(mf) / (1 + exp(mf))
      prob.array[[iter/10000]] <- prob 
      if ((iter/10000) > length(misfit)) {
        misfit <- c(misfit, rep(0, times = 100))
        misfit[iter/10000] <- sum(sum((x - theta) ^ 2))
      }
    }
  }
  
  return(prob)
}