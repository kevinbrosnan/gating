#### MACSI SFI Review Poster ####

source('code/startup.R')

## Existing Approach ##

par(pty = "s")
plot(x = rit_10bit[,1], y = rit_10bit[,2], xlab = "FSC.H", ylab = "SSC.H", las = 1, col = "gray")
polygon(x = c(50,50,440,380), y = c(0,200,580,80), border = "red", lwd = 2, lty = 1)
polygon(x = c(50,50,420,600), y = c(0,200,610,110), border = "blue", lwd = 2, lty = 2)
legend('topright', legend = c('Expert 1', 'Expert 2'), col = c('red', 'blue'), lty = c(1,2), lwd = c(2,2))

## Results ##
rit_6bit_image <- as_image(rit_6bit[,1:2], cur.dim = 6)
par(pty = "s")
image(-rit_6bit_image, xaxt = 'n', yaxt = 'n', ann = FALSE, col = grey.colors(2, start = 0, end = 1))
denois_rit_6bit <- ising_denosing(rit_6bit_image)
par(pty = "s")
image(-denois_rit_6bit$mf, yaxt = 'n', xaxt = 'n', ann = FALSE, col = grey.colors(10, start = 0, end = 1))
rit_6bit_clust <- ConnCompLabel(denois_rit_6bit$prob)
rit_6bit_clust[rit_6bit_clust > 0] <- 1
par(pty = "s")
image(-rit_6bit_clust, xaxt = 'n', yaxt = 'n', ann = FALSE, col = grey.colors(10, start = 0, end = 1))
