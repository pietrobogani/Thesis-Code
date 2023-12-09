PlotQRbands <- function(BQ, bBQ, B2, bB2, QQ, ylims) {


  qq <- c(0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975)
  Nboot <- dim(bB2)[3]
  
  cat("Dimensions of bBQ:", dim(bBQ), "\n")  
  
  qBQ <- t(apply(bBQ, 1, function(x) quantile(x, probs = qq)))
  qB2 <- quantile(bB2, probs = qq)
  
  #print(head(qBQ))
  #print(head(bBQ))
  #print(head(qB2))
  
  
  # For area plot
  mat_quant <- qBQ
  matm <- mat_quant
  
  for (i in 2:ncol(matm)) {
    matm[, i] <- matm[, i] - mat_quant[, i - 1]
    if (is.na(matm[nrow(matm), i])) {
      matm[nrow(matm), i] <- matm[nrow(matm) - 1, i]
    }
  }
  if (is.na(matm[nrow(matm), 1])) {
    matm[nrow(matm), 1] <- matm[nrow(matm) - 1, 1]
  }
  
  colors <- colorRampPalette(c(rgb(1,1,1), rgb(0.75,0.75,0.75)))(7)
  
  # Plotting

  plot(QQ, BQ, type = "l", xlab = expression(tau), ylab = expression(beta(tau)), 
       xlim = c(QQ[1], QQ[length(QQ)]), ylim = ylims)

  
  for (i in 7:1) {
    if (i == 1) {
      y_values <- c(mat_quant[,i], rep(0, length(QQ)))
    } else {
      y_values <- c(mat_quant[,i], rev(mat_quant[,i-1]))
    }
    

    
    polygon(c(QQ, rev(QQ)), y_values, col = colors[i], border = NA)
  }
  
  
  lines(QQ, BQ, col = "red", lwd = 2)
  lines(QQ, mat_quant[,4], col = "black", lty = 2, lwd = 2)
  lines(QQ, rep(B2, length(QQ)), col = "blue", lty = 2, lwd = 2)
  abline(h = 0, col = "black", lwd = 1)
  
  # legend(legendLocation, legend = c("In-sample fit", "Median", "OLS"),
  #        col = c("red", "black", "blue"), lty = c(1, 2, 2), lwd = c(2, 2, 2))
  # 
  
}
