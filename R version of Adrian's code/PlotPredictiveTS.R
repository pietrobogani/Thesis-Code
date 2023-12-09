PlotPredictiveTS <- function(YQ, QQ, Time, yh) {

  # Ensure YQ, QQ, and Time are of the same length
  stopifnot(NCOL(YQ) == length(QQ), NROW(YQ) == length(Time))

  # Select a few important quantiles: 5%, 10%, 25%, 50%, 75%, 90%, 95%
  Qsel <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
  PosSel <- sapply(Qsel, function(q) which.min(abs(QQ - q)))
  mat_quant <- YQ[, PosSel]

  # Check for non-finite values and omit them
  finite_inds <- (rowSums(is.finite(mat_quant)) == length(Qsel))
  Time <- Time[finite_inds]
  yh <- yh[finite_inds]
  mat_quant <- mat_quant[finite_inds, ]

  # Establish plotting area
  plot(Time, yh, type = "n", xlim = range(Time), ylim = range(c(mat_quant, yh), na.rm = TRUE), xlab = "", ylab = "", xaxt = "n")

  # Fill areas between quantiles
  polygon(c(Time, rev(Time)), c(mat_quant[,1], rev(mat_quant[,2])), col = gray(0.95), border = NA)
  polygon(c(Time, rev(Time)), c(mat_quant[,2], rev(mat_quant[,3])), col = gray(0.90), border = NA)
  polygon(c(Time, rev(Time)), c(mat_quant[,3], rev(mat_quant[,5])), col = gray(0.75), border = NA)
  polygon(c(Time, rev(Time)), c(mat_quant[,5], rev(mat_quant[,6])), col = gray(0.90), border = NA)
  polygon(c(Time, rev(Time)), c(mat_quant[,6], rev(mat_quant[,7])), col = gray(0.95), border = NA)

  # Overlay the median and realized values
  lines(Time, yh, col = "blue", lty = 2)
  lines(Time, mat_quant[,4], col = "black")

  # X-axis date labels
  years <- as.POSIXlt(Time)$year + 1900
  XTickTimes <- Time[years %% 5 == 0]
  axis(1, at = XTickTimes, labels = years[years %% 5 == 0])

  # Legend
  legend("topright", legend = c("Realized", "Median"), lty = c(2, 1), col = c("blue", "black"))
}


