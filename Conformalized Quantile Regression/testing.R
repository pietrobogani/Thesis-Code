  library(ggplot2)
  library(MASS)
  library(quantregForest) 


# Set seed for reproducibility
set.seed(1)

# Desired miscoverage error
alpha <- 0.1
# Low and high target quantiles
quantiles <- c(5, 95)

# Maximal number of test points to plot
max_show <- 1000

# Save figures?
save_figures <- FALSE

# Parameters of random forests
n_estimators <- 100
min_samples_leaf <- 40
max_features <- 1 # 1D signal
random_state = 0

plot_func <- function(x, y, y_u=NULL, y_l=NULL, pred=NULL, shade_color="",
                      method_name="", title="", filename=NULL,
                      save_figures=FALSE) {
  
  x_ <- x[1:max_show]
  y_ <- y[1:max_show]
  if (!is.null(y_u)) {
    y_u_ <- y_u[1:max_show]
  }
  if (!is.null(y_l)) {
    y_l_ <- y_l[1:max_show]
  }
  if (!is.null(pred)) {
    pred_ <- pred[1:max_show]
  }
  
  p <- ggplot() +
    geom_point(data=data.frame(x=x_, y=y_), aes(x=x, y=y), color="black", alpha=0.2, size=3) +
    xlim(range(x_)) +
    ylim(c(-2.5, 7)) +
    xlab("$X$") +
    ylab("$Y$") +
    theme_minimal()
  
  if (!is.null(y_u) && !is.null(y_l)) {
    p <- p + geom_ribbon(data=data.frame(x=x_, y_u=y_u_, y_l=y_l_), aes(x=x, ymin=y_l, ymax=y_u),
                         fill=shade_color, alpha=0.3) +
      labs(title=title, subtitle=method_name + " prediction interval") +
      theme_minimal()
  }
  
  if (!is.null(pred)) {
    if (length(dim(pred_)) == 2) {
      p <- p + geom_line(data=data.frame(x=x_, y_l=pred_[,1], y_u=pred_[,2]), aes(x=x, y=y_l), color="black", size=1) +
        labs(title=title, subtitle="Predicted low and high quantiles") +
        theme_minimal()
    } else {
      p <- p + geom_line(data=data.frame(x=x_, y=pred_), aes(x=x, y=y), linetype="dashed", color="black", size=1) +
        labs(title=title, subtitle="Predicted value") +
        theme_minimal()
    }
  }
  
  print(p)
  
  if (save_figures && !is.null(filename)) {
    ggsave(filename, p, width=6, height=4, dpi=300)
  }
}




# Number of training examples
n_train <- 2000
# Number of test examples (to evaluate average coverage and length)
n_test <- 5000

# Function to construct data (1D example)
f <- function(x) {
  ax <- rep(0, length(x))
  for (i in seq_along(x)) {
    ax[i] <- rpois(1, lambda = sin(x[i])^2 + 0.1) + 0.03 * x[i] * rnorm(1)
    ax[i] <- ax[i] + 25 * (runif(1) < 0.01) * rnorm(1)
  }
  return(as.numeric(ax))
}

# Training features
x_train <- runif(n_train, 0, 5.0)

# Test features
x_test <- runif(n_test, 0, 5.0)

# Generate labels
y_train <- f(x_train)
y_test <- f(x_test)

# Reshape the features
x_train <- matrix(x_train, ncol = 1)
x_test <- matrix(x_test, ncol = 1)

# Display the test data in full range (including the outliers)
library(ggplot2)
fig <- ggplot(data.frame(x = x_test, y = y_test)) +
  geom_point(aes(x = x, y = y), color = "black", alpha = 0.3, size = 3, fill = NA) +
  labs(x = "$X$", y = "$Y$", title = "Test data (visualize outliers)") +
  theme_minimal()

print(fig)

# Save figure if specified
if (save_figures) {
  ggsave("illustration_test_data.png", fig, width = 6, height = 4, dpi = 300)
}

# Display the test data without outliers (zoom in)
plot_func(x_test, y_test, title = "Test data (zoom in)")


#---------------------- fin qui ok



# Set seed for reproducibility
set.seed(1)

# Divide the data into proper training set and calibration set
idx <- sample(seq_len(n_train))
n_half <- floor(n_train / 2)
idx_train <- idx[1:n_half]
idx_cal <- idx[(n_half + 1):(2 * n_half)]

# Create training and calibration sets
x_cal <- x_train[idx_cal, , drop = FALSE]
y_cal <- y_train[idx_cal]


# Define quantile random forest (QRF) parameters
params_qforest <- list(
  ntree = n_estimators,
  minsize = min_samples_leaf,
  mtry = max_features
)

# Train the quantile random forest model
qrf_model <- quantregForest(x = x_train, y = y_train, control = params_qforest)

# Predict quantiles on the calibration set
quantile_pred <- predict(qrf_model, x_cal, what = c(0.05,0.95))

# Initialize a vector for errors
E_i <- rep(NA, length(idx_cal))

# Calculate errors for each point in the test set I2
for (i in 1:length(E_i)) {
  E_i[i] <- max(quantile_pred[i,2] - y_cal[i], y_cal[i] - quantile_pred[i,3])

}

# Compute Q(1-??)(E, I2)
quantile_E <- quantile(E_i, pmin(1, pmax(0, (1 - QQ[jq]) * (1 + 1/length(Yh2)))))

























#bunc_low <- rq(Yh1[(h + 1):length(Yh)] ~ 1, tau=0.001)
YQunclow_IS[(h + 1):length(Yh), jq] <- - Inf #rep(coef(bunc_low), length(Time) - h)

bunc_high <- rq(Yh1[(h + 1):length(Yh)] ~ 1, tau=QQ[jq])
YQunchigh_IS[(h + 1):length(Yh), jq] <- rep(coef(bunc_high), length(Time) - h)

# Initialize a vector for errors
E_i <- rep(NA, length(Yh2))

# Calculate errors for each point in the test set I2
for (i in 1:length(E_i)) {
  E_i[i] <- max(YQunclow_IS[h + i, jq] - Yh2[i], Yh2[i] - YQunchigh_IS[h + i, jq])
}

# Compute Q(1-??)(E, I2)
quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))

YQunc_IS[(h + 1):length(Yh), jq] <- rep(coef(bunc_high), length(Time) - h) + quantile_E

#da aggiungere
YQGDPonly_high_OOS
YQGDPonly_low_OOS
YQunclow_OOS
YQunchigh_OOS

YQunchigh_OOS
YQunclow_OOS