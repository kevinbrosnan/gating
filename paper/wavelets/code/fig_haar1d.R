
library(tikzDevice)
tikz(width = 7, height = 3.5, file = "fig_haar1d.tex")
par(mfrow = c(1,2))

# Father Wavelet
x.father <- seq(0, 1, length.out = 100)
y.father <- rep(1, times = length(x.father))
plot(x = x.father, y = y.father, type = "l", xlim = c(0, 1), ylim = c(-1,1),
     xlab = "$t$", ylab = "", bty = "l", xaxp = c(0, 1, 2), yaxp = c(-1,1,2))
segments(x0 = c(0,1), y0 = c(0,0), y1 = c(1,1), lty = 2)
abline(h = 0, col = "gray", lty = 2)
mtext(text = "Father Wavelet - $\\phi(t)$")

# Mother Wavelet
plot(x = 1, type = "n", xlim = c(0,1), ylim = c(-1,1), xlab = "$t$", 
     ylab = "", bty = "l", yaxp = c(-1,1,2), xaxp = c(0,1,2))
mtext(text = "Mother Wavelet - $\\psi(t)$")
segments(x0 = c(0, 0, 0.5, 0.5, 1), y0 = c(0, 1, 1, -1, -1), 
         x1 = c(0, 0.5, 0.5, 1, 1), y1 = c(1, 1, -1, -1, 0), 
         lty = c(2, 1, 2, 1, 2))
abline(h = 0, col = "gray", lty = 2)
dev.off()
par(mfrow = c(1,1))
  