rm(list = ls())

n = 473
pi = 0.5
Vraw = 0.099

htau = 0.211
aa = Vraw^2
bb = 2 * sqrt(8) * Vraw / sqrt(n) * sqrt(pi * (1 - pi))
cc = 8 / n * pi * (1 - pi)

alphas = 1/(2:1000)

bounds = sapply(alphas, function(alpha){
  zalpha = -qnorm(alpha/2)
  high = uniroot(function(tt) ((htau - tt)^2 - zalpha^2 * (aa + bb * abs(tt) + cc * tt^2)),
                 c(htau, 10))$root
  low = uniroot(function(tt)((htau - tt)^2 - zalpha^2 * (aa + bb * abs(tt) + cc * tt^2)),
                c(htau, -10))$root
  c(low=low, high=high, low0=htau-zalpha*sqrt(aa), high0=htau+zalpha*sqrt(aa))
})

pdf("sensitivity.pdf")
pardef = par(mar = c(5, 4.5, 2, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
cols = RColorBrewer::brewer.pal(3, "Set1")
plot(NA, NA, xlim = range(1/alphas), ylim = range(bounds), log = "x",
     xlab = "1/alpha", ylab = "tau")
lines(1/alphas, rep(htau, length(alphas)), lwd = 2, col = cols[3])
abline(h=0, lty = 3)
#abline(v=20, lty = 3)
lines(1/alphas, bounds["low",], lwd = 2, col = cols[1])
lines(1/alphas, bounds["high",], lwd = 2, col = cols[1])
lines(1/alphas, bounds["low0",], lwd = 2, col = cols[2], lty = 2)
lines(1/alphas, bounds["high0",], lwd = 2, col = cols[2], lty = 2)
par=pardef
dev.off()

alphas[19]
round(bounds[,19], 3)
