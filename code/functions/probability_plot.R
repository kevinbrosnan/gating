probability_map <- function(x, data = NULL) {
  
  no.row <- nrow(x)
  no.col <- ncol(x)
  
  colours.scale <- c("darkgreen", "green", "yellow", "orange", "red")
  
  x.pos <- rep(1:no.row, times = no.col)
  y.pos <- sort(rep(1:no.col, times = no.row))
  val <- as.vector(x)
  val.col <- as.factor(val)
  
  par(pty = "s")
  plot(x.pos, y.pos, type = "n", ylim = c(-0.5, no.col), xlim = c(-0.5, no.col), las = 1)
  
  for (i in 1:length(x.pos)) {
    polygon(x = c((x.pos[i] - 1.5), (x.pos[i] - 1.5), (x.pos[i] - 0.5), (x.pos[i] - 0.5)), 
            y = c((y.pos[i] - 1.5), (y.pos[i] - 0.5), (y.pos[i] - 0.5), (y.pos[i] - 1.5)),
            col = colours.scale[val.col[i]], border = NA)
  }
  
  if (!is.null(data)) {
    points(data)
  }
  
  return(invisible(NULL))
}