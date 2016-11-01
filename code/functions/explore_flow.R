
explore_flow <- function(flow_data) {
  
  panel.hist.density <- function(x) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h <- hist(x, plot = FALSE, breaks = 128)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan")
    d <- density(x, na.rm = TRUE, bw = "nrd0", adjust = 1.2)
    d$y <- d$y/max(d$y)
    lines(d, col = "red", lwd = 2)
    rug(x)
  }
  
  panel.cor <- function(x, y, digits = 2, prefix = "") {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    paste0("$", txt, "$")
    text(0.5, 0.5, txt, cex = 3)
  }
  
  pairs(flow_data,  upper.panel = panel.cor, diag.panel = panel.hist.density)
  return(invisible(NULL))
}