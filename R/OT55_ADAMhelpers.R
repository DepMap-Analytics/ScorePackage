ADAM.tradeoffEO_TPR<-function (EO, TPR, point = NULL, filename)
{
  if (length(point) == 0) {
    CCOL <- "red"
  }
  else {
    CCOL <- rgb(255, 0, 0, alpha = 100, maxColorValue = 255)
  }
  x <- EO
  x[x == Inf] <- max(x[x < Inf])
  x[is.na(x)] <- 0
  x <- (x - min(x))/(max(x) - min(x))
  y <- TPR
  y <- (y - min(y))/(max(y) - min(y))
  orEO <- EO
  orEO[orEO == Inf] <- max(orEO[orEO < Inf])
  orEO[is.na(orEO)] <- 0
  orTPR <- TPR
  EO <- x
  TPR <- y
  #pdf(paste0(filename, "_Tradeoff.pdf"))
  par(mar = c(4, 4, 4, 4))
  MAIN <- c("log10 (obs/Expct) n.genes [red, left]", paste("% covered ",
                                                           filename, " [blue, right]", sep = ""))
  plot(EO, type = "l", xlab = "genes depleted in >= # cell lines",
       ylab = "", axes = FALSE, lwd = 4, main = MAIN, col = CCOL,
       cex.main = 0.8, xlim = c(0, length(EO)))
  axis(2, at = seq(0, 1, 0.2), format(seq(min(orEO), max(orEO),
                                          (max(orEO) - min(orEO))/5), digits = 2))
  axis(1)
  par(new = TRUE)
  plot(TPR, type = "l", xlab = "", ylab = "", axes = FALSE,
       lwd = 4, col = "blue", ylim = c(0, 1), xlim = c(0, length(EO)))
  axis(4, at = seq(0, 1, 0.2), format(seq(min(orTPR), max(orTPR),
                                          (max(orTPR) - min(orTPR))/5), digits = 2))
  if (length(point) == 0) {
    point <- min(which(!y > x))
    abline(v = point)
    abline(h = y[point], lty = 2)
    points(point, y[point], pch = 16, cex = 2)
  }
  else {
    abline(v = point)
    abline(h = y[point], lty = 2)
    points(point, y[point], pch = 16, cex = 2)
  }
  legend("top", paste(format(100 * orTPR[point], digits = 2),
                      "% covered", sep = ""), bg = NULL, bty = "n")
  #dev.off()
  #print(point)
  return(point)
}
