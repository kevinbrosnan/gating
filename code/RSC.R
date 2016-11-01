####---- Code for use in RSC 2016 ----####

source('code/startup.R')

## Flow cytometry data (Slide 5)
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

panel.cor <- function(x, y, digits = 2, prefix = "")
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  paste0("$", txt, "$")
  text(0.5, 0.5, txt, cex = 3)
}

pairs(rit_10bit[,1:3],  upper.panel = panel.cor, diag.panel = panel.hist.density)

stargazer(rit_10bit, summary = TRUE, digits = 2, label = "tab:ritdatasum", 
          title = "Summary Statistics for Rituximab Data", align = TRUE)

## Sparsity (Slide 6)
pairs(rit_8bit[,1:3],  upper.panel = panel.cor, diag.panel = panel.hist.density)

stargazer(rit_8bit, summary = TRUE, digits = 2, label = "tab:ritdatasum", 
          title = "Summary Statistics for Rituximab Data", align = TRUE)

## Industry Standard (Slide 8)
par(pty="s")
plot(x = rit_10bit[,1], y = rit_10bit[,2], xlab = "FSC.H", ylab = "SSC.H", las = 1, col = "gray")
polygon(x = c(50,50,440,380), y = c(0,200,580,80), border = "red", lwd = 2, lty = 1)
#plot(x = rit_10bit[,1], y = rit_10bit[,2], xlab = "FSC.H", ylab = "SSC.H", las = 1, col = "gray")
polygon(x = c(50,50,420,600), y = c(0,200,610,110), border = "blue", lwd = 2, lty = 2)
legend('topright', legend = c('Expert 1', 'Expert 2'), col = c('red', 'blue'), lty = c(1,2), lwd = c(2,2))

## Model Based Clustering (Slide 9)
mbc_4 <- Mclust(data = rit_10bit[,1:2], G = 4)
mbc_6 <- Mclust(data = rit_10bit[,1:2], G = 6)
par(mfrow = c(1,2))
plot(mbc_4, what = "classification", main = "4 Groups by mclust")
title(main = "4 groups by mclust")
plot(mbc_6, what = "classification", main = "6 Groups by mclust")
title(main = "6 groups by mclust")
par(mfrow = c(1,1))

## t-mixtures (Slide 10)
model <- flowClust(rit_10bit, varNames = c("FSC.H", "SSC.H"), K = 1:6)
best <- which.max(criterion(object = model, "BIC"))
plot(model[[best]], data = rit_10bit, las = 1, xlab = "FSC.H", ylab = "SSC.H")

## Markov Random Fields (Slide 13)
X <- read.bmp("data/letters/lettera.bmp")
X <- t(X)[,nrow(X):1] 
d <- 2 * ((X - mean(X)) > 0) - 1
Y <- d + rnorm(length(d))
png(paste0(getwd(), '/images/letterA/noisy.png'))
image(Y, yaxt = 'n', xaxt = 'n', ann = FALSE)
title(main = paste0("Noisy Image"))
dev.off()
denois_letterA <- ising_denosing(Y)
for (i in 1:length(denois_letterA)) {
  iter <- i * 10000 + 1
  png(paste0(getwd(), "/images/letterA/", "fig_", i, ".png"))
  image(denois_letterA[[i]], main = paste0("Iteration ", iter), yaxt = 'n', xaxt = 'n', ann = FALSE)
  title(main = paste0("Iteration ", iter))
  dev.off()
}

## Results (Slide 17)
rit_8bit_image <- as_image(rit_8bit[,1:2], cur.dim = 8)
png(paste0(getwd(), '/images/rit_8bit/original.png'))
par(pty = "s")
image(-rit_8bit_image, xaxt = 'n', yaxt = 'n', ann = FALSE, col = grey.colors(n = 2, start = 0, end = 1))
title(main = "Original 8-bit")
dev.off()
denois_rit_8bit <- ising_denosing(rit_8bit_image)
for (i in 1:length(denois_rit_8bit)) {
  iter <- i * 10000 + 1
  png(paste0(getwd(), "/images/rit_8bit/", "fig_", i, ".png"))
  image(denois_rit_8bit[[i]], main = paste0("Iteration ", iter), yaxt = 'n', xaxt = 'n', ann = FALSE, col = grey.colors(n = 10, start = 0, end = 1))
  title(main = paste0("Iteration ", iter))
  dev.off()
}
rit_8bit_clust <- ConnCompLabel(denois_rit_8bit)
png(paste0(getwd(), '/images/rit_8bit/clustered.png'))
image(-rit_8bit_clust, xaxt = 'n', yaxt = 'n', ann = FALSE)
title(main = "Clustered 8-bit")
dev.off()

rit_7bit_image <- as_image(rit_7bit[,1:2], cur.dim = 7)
png(paste0(getwd(), '/images/rit_7bit/original.png'))
image(rit_7bit_image, xaxt = 'n', yaxt = 'n', ann = FALSE)
title(main = "Original 7-bit")
dev.off()
denois_rit_7bit <- ising_denosing(rit_7bit_image)
for (i in 1:length(denois_rit_7bit)) {
  iter <- i * 10000 + 1
  png(paste0(getwd(), "/images/rit_7bit/", "fig_", i, ".png"))
  image(denois_rit_7bit[[i]], main = paste0("Iteration ", iter), yaxt = 'n', xaxt = 'n', ann = FALSE)
  title(main = paste0("Iteration ", iter))
  dev.off()
}
rit_7bit_clust <- ConnCompLabel(denois_rit_7bit[[40]])
png(paste0(getwd(), '/images/rit_7bit/clustered.png'))
image(rit_7bit_clust, xaxt = 'n', yaxt = 'n', ann = FALSE)
title(main = "Clustered 7-bit")
dev.off()

rit_6bit_image <- as_image(rit_6bit[,1:2], cur.dim = 6)
png(paste0(getwd(), '/images/rit_6bit/original.png'))
image(-rit_6bit_image, xaxt = 'n', yaxt = 'n', ann = FALSE, col = grey.colors(2, start = 0, end = 1))
title(main = "Original 6-bit")
dev.off()
denois_rit_6bit <- ising_denosing(rit_6bit_image)
for (i in 1:length(denois_rit_6bit)) {
  iter <- i * 10000 + 1
  png(paste0(getwd(), "/images/rit_6bit/", "fig_", i, ".png"))
  image(denois_rit_6bit[[i]], main = paste0("Iteration ", iter), yaxt = 'n', xaxt = 'n', ann = FALSE)
  title(main = paste0("Iteration ", iter))
  dev.off()
}
rit_6bit_clust <- ConnCompLabel(denois_rit_6bit)
png(paste0(getwd(), '/images/rit_6bit/clustered.png'))
image(rit_6bit_clust, xaxt = 'n', yaxt = 'n', ann = FALSE)
title(main = "Clustered 6-bit")
dev.off()