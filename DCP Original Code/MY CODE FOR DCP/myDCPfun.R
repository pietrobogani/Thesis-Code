library(quantreg)

fullDCP <- function(Y, X, X_new, alpha, ytrial) {

  cumulative_prob <- function(X, quantiles, values) {
    if (length(values[values <= X]) == 0) {
      return(0) # If X is less than the smallest value, the cumulative probability is 0
    }

    max_value_not_exceeding_X <- max(values[values <= X], na.rm = TRUE)

    quantile_index <- which(values == max_value_not_exceeding_X)
    if (length(quantile_index) == 0) {
      return(1) # If X is greater than all values, return the maximum cumulative probability
    }

    return(quantiles[min(quantile_index)])
  }



  QQ <- seq(0.05, 0.95, by = 0.05)
  Xadd <- as.matrix(c(X, X_new))
  b <- array(NA, dim = c(length(ytrial), 2, length(QQ)))
  U <- array(NA, dim = c(length(ytrial), length(Xadd)))
  qqTarg <- array(NA, dim = c(length(Xadd),length(QQ)))
  pval <- array(NA, dim = c(length(ytrial)))

  # Compute F^

    for (y in 1:length(ytrial)){

      Yadd <- c(Y, ytrial[y])

      # Compute F^
      for (jq in 1:length(QQ)) {
        beta <- rq(Yadd ~ Xadd, tau = QQ[jq])
        b[y, , jq] <- coef(beta)
        qqTarg[,jq] <-  as.vector(cbind(1, Xadd) %*% b[y, , jq])
      }

      # Compute ranks
      for (x in 1:length(Xadd)) {
        qqTarg[x, ] <- sort(qqTarg[x, ])
         # if (is.unsorted(qqTarg[x, ])) {
         #  cat("Error: 'qqTarg[x, ]' should be sorted.", x)  #check if we're doing good.
         # }
         #
        U[y, x] <-  cumulative_prob(Yadd[x], QQ, qqTarg[x, ] )

      }

      # Compute p-value
      pval[y] <- mean(U[y, ] >= U[y, length(Xadd)]) #I select the last element of U (for the given y)
    }

return(pval)
}



# Test the function

n <- 204
jtFirstOOS <- 80 #First index for out-of-sample computations



#----Generate AR(2)
phi_ar2 <- c(0.5, -0.2)  # AR coefficients for AR(2)
Y_ar2 <- numeric(2*n)
Y_ar2[1] <- 1
Y_ar2[2] <- -2
for (i in 3:(2*n)) {
  Y_ar2[i] <- phi_ar2[1] * Y_ar2[i-1] + phi_ar2[2] * Y_ar2[i-2] + rt(n = 1, df = 2) #heavy tails
}
Y_ar2 <- Y_ar2[n: (2*n - 1)] #I implemented a burn-in

# Create lagged variables
Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])


# Prepare the dataset for AR(3) (excluding the first three NA values)
data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)], Y_lag2 = Y_lag2_ar2[-c(1:3)])
alpha <- 0.1

in_or_out<- numeric(length(jtFirstOOS:(n-1)))

#for (i in jtFirstOOS:(n-1)) {
i <- jtFirstOOS+10
Y <- data_ar2$Y[1:i]
X <- data_ar2$Y_lag1[1:i]
X_new <- data_ar2$Y_lag1[1+i]
Y_new <- data_ar2$Y[1+i]
ytrial <- seq(-max(abs(Y)), max(abs(Y)), 0.1)

pvalues <- fullDCP(Y, X, X_new, alpha, ytrial)

indices <- which(pvalues > alpha)
CI <- c(ytrial[indices[1]], ytrial[indices[length(indices)]])
# if (is.na(CI[1])) {
#   in_or_out[i - jtFirstOOS + 1] = 0
# } else if (Y_new >= CI[1] && Y_new <= CI[2]) {
#   in_or_out[i - jtFirstOOS + 1] = 1
# } else {
#   in_or_out[i - jtFirstOOS + 1] = 0
# }

# }
#
# length(which(in_or_out==1))
# length(which(in_or_out==0))












library(ggplot2)
data("diamonds")


n <-500 #length(diamonds[,1])
jtFirstOOS <- 390
alpha <- 0.1
in_or_out<- numeric(length(jtFirstOOS:(n-1)))

for (i in jtFirstOOS:(n-1)) {

X <- diamonds$x[1:i]
Y <- diamonds$carat[1:i]
X_new <- diamonds$x[i+1]
Y_new <- diamonds$carat[i+1]
ytrial <- seq(-max(abs(Y)), max(abs(Y)), 0.01)
pvalues <- fullDCP(Y, X, X_new, alpha, ytrial)

indices <- which(pvalues > alpha)
CI <- c(ytrial[indices[1]], ytrial[indices[length(indices)]])
print(CI)

if (is.na(CI[1])) {
  in_or_out[i - jtFirstOOS + 1] = 0
} else if (Y_new >= CI[1] && Y_new <= CI[2]) {
  in_or_out[i - jtFirstOOS + 1] = 1
} else {
  in_or_out[i - jtFirstOOS + 1] = 0
}


print(in_or_out[i-jtFirstOOS+1])
}
print(sum(in_or_out)/length(in_or_out))









n=100
x=runif(n)
u=rnorm(n)
y=x+x*u
ytrial <- seq(-max(abs(y)), max(abs(y)), 0.01)
alpha<-0.1
x0 <- 0.8
pvalues <- fullDCP(y, x, x0, alpha, ytrial)
indices <- which(pvalues > alpha)
CI <- c(ytrial[indices[1]], ytrial[indices[length(indices)]])
print(CI)


ynew <- 0.8 + 0.8*u
mean(ynew <=CI[2])
