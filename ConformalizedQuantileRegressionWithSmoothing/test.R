# Install and load necessary packages
library(quantreg)
library(lubridate)
library(pracma)
library(readxl)
library(sn)

# Clear workspace 
rm(list = ls())

# Set forecast horizon (run script separately for h = 1 and h = 4)
h <- 4

loadsavedresults = FALSE; # If I ran code already results are stored in ResOOS_H11 and ResOOS_H14, in that case, set = TRUE


# Graphics settings - R's graphics system differs from MATLAB's, so some modifications are necessary
par(mfrow = c(1, 1))  # Reset plot window to single pane


# Load data 
file_path <- "DataVulnerabilityAppendix.xls"

# Read the file
data <- read_excel(file_path)
data<-data[,1:3]

# Filter data for 1973Q1-2015Q4
colnames(data)[1] <- "Time"
data$Time <- as.Date(data$Time)


# Subset the data
data <- data[data$Time >= as.Date("1973-01-01") & data$Time <= as.Date("2015-10-01"), ]
X <- data[,2:3]
Time <- data$Time


# Set forecast settings
QQ <- seq(0.05, 0.95, by = 0.05)
deltaYY <- 0.1
YY <- seq(-20, 20, by = deltaYY)
jtFirstOOS <- which(year(data$Time) == 1993 & month(data$Time) == 1)
indices <- which(QQ %in% c(0.05, 0.25, 0.5, 0.75, 0.95))
jq05 <- indices[1]
jq25 <- indices[2]
jq50 <- indices[3]
jq75 <- 15 #couldn't automatically translate from MATLAB, I set it manually
jq95 <- 19 #couldn't automatically translate from MATLAB, I set it manually

# Construct average growth rates
y <- X$A191RL1Q225SBEA
Yh <- matrix(0, nrow=length(y), ncol=4)


Yh <- filter(y, rep(1/h, h), sides=1) #If h = 1, y = Yh
if (h>1){
  Yh[1:(h-1)] <- NA
}

#Construct matrices of regressors
Z <- cbind(1, X[,2], y)
ZGDPonly <- cbind(1, y)
Z <-as.matrix(Z)


# Get length of Time and QQ/YY
len_time <- length(data$Time)
len_qq <- length(QQ)
len_yy <- length(YY)

# Initialize matrices to store forecasts
{
  # Raw quantiles
  YQ_NaNs <- matrix(NA, len_time, len_qq)
  #YQ_low_adj_IS <- YQ_NaNs #IN the original script I've just YQ_IS, because I do a point estimate of the quantile. Here I do a confidence interval estimation and I need 2
  #YQ_high_adj_IS <- YQ_NaNs
  YQ_low_IS <- YQ_NaNs
  YQ_high_IS <- YQ_NaNs
  YQ_IS <- YQ_NaNs
  #YQ_low_adj_OOS <- YQ_NaNs 
  #YQ_high_adj_OOS <- YQ_NaNs
  YQ_low_OOS <- YQ_NaNs
  YQ_high_OOS <- YQ_NaNs
  YQ_OOS <- YQ_NaNs
  YQ_lowGDPonly_IS <- YQ_NaNs
  YQ_highGDPonly_IS <- YQ_NaNs
  YQGDPonly_IS <- YQ_NaNs
  YQGDPonly_OOS <- YQ_NaNs
  YQunc_IS <- YQ_NaNs
  YQunc_OOS <- YQ_NaNs
  YQunclow_IS <- YQ_NaNs
  YQunchigh_IS <- YQ_NaNs
  YQGDPonly_high_OOS <- YQ_NaNs
  YQGDPonly_low_OOS <- YQ_NaNs
  YQunclow_OOS <- YQ_NaNs
  YQunchigh_OOS <- YQ_NaNs
  
  YQunchigh_OOS <- YQ_NaNs
  YQunclow_OOS <- YQ_NaNs
  
  
  # PDFs (evaluated over grid)
  P_NaNs <- matrix(NA, len_time, len_yy)
  PST_IS <- P_NaNs
  PST_OOS <- P_NaNs
  PSTGDPonly_IS <- P_NaNs
  PSTGDPonly_OOS <- P_NaNs
  PSTunc_IS <- P_NaNs
  PSTunc_OOS <- P_NaNs
  
  # Smoothed quantiles
  Q_NaNs <- matrix(NA, len_time, len_qq)
  QST_IS <- Q_NaNs
  QST_OOS <- Q_NaNs
  QSTGDPonly_IS <- Q_NaNs
  QSTGDPonly_OOS <- Q_NaNs
  QSTunc_IS <- Q_NaNs
  QSTunc_OOS <- Q_NaNs
  
  # CDFs (evaluated over grid)
  C_NaNs <- matrix(NA, len_time, len_yy)
  CST_IS <- C_NaNs
  CST_OOS <- C_NaNs
  CSTGDPonly_IS <- C_NaNs
  CSTGDPonly_OOS <- C_NaNs
  CSTunc_IS <- C_NaNs
  CSTunc_OOS <- C_NaNs
  
  # Skewed t-distribution parameters
  STpar_NaNs <- matrix(NA, len_time, 4)
  STpar_IS <- STpar_NaNs
  STpar_OOS <- STpar_NaNs
  STparGDPonly_IS <- STpar_NaNs
  STparGDPonly_OOS <- STpar_NaNs
  STparunc_IS <- STpar_NaNs
  STparunc_OOS <- STpar_NaNs
  
  # Predictive scores
  Score_NaNs <- rep(NA, len_time)
  ScoreST_IS <- Score_NaNs
  ScoreST_OOS <- Score_NaNs
  ScoreSTGDPonly_IS <- Score_NaNs
  ScoreSTGDPonly_OOS <- Score_NaNs
  ScoreSTunc_IS <- Score_NaNs
  ScoreSTunc_OOS <- Score_NaNs
  
  # Probability integral transforms
  Pit_NaNs <- rep(NA, len_time)
  PitST_IS <- Pit_NaNs
  PitST_OOS <- Pit_NaNs
  PitSTGDPonly_IS <- Pit_NaNs
  PitSTGDPonly_OOS <- Pit_NaNs
  PitSTunc_IS <- Pit_NaNs
  PitSTunc_OOS <- Pit_NaNs
  
  # Left entropy
  Entropy_NaNs <- rep(NA, len_time)
  LeftEntropy_IS <- Entropy_NaNs
  LeftEntropy_OOS <- Entropy_NaNs
  
  #Split creating I1 and I2
  full_length <- length(Yh[(h + 1):length(Yh)])
  
  #permuted_indices <- c(1,2,3,4,sample((h+1):length(Yh))) #if you want to permute
  permuted_indices <- c(1:length(Yh)) # if you don't want to permute
  
  Yh_perm <- Yh[permuted_indices]
  Z_perm <- Z[permuted_indices,]
  ZGDPonly_perm <- ZGDPonly[permuted_indices,]
  Yh1 <- Yh_perm[(h+1):(h+full_length/2)]
  Yh2 <- Yh_perm[(h+1+full_length/2):length(Yh)]
  Z1 <- Z_perm[1:(full_length/2),]
  Z2 <- Z_perm[(full_length/2+1):(length(Yh) - h),]
  ZGDPonly1 <- ZGDPonly_perm[1:(full_length/2),]
  ZGDPonly2 <- ZGDPonly_perm[(full_length/2+1):(length(Yh) - h),]
}


  #-------------------    %% In-sample estimation of conditional quantiles
  
  for (jq in 1:length(QQ)) {
    
    #CQR is built for intervals, not quantiles. Thus I create a prediction interval of the form [-Inf, quantile X] creating a prediction interval with
    #confidence X %
    
    
    #------ Conformalized Quantile regression with both NFCI and GDP
    
    #b_low <- rq(Yh1 ~ Z1[,-1], tau= 0.01) #Train on I1
    #YQ_low_IS[(h + 1):(h + full_length/2), jq] <-as.vector(Z2 %*% coef(b_low)) #Evaluate on I2
    YQ_low_IS[(h + 1):(h + full_length/2), jq] <- -Inf 
    
    
    b_high <- rq(Yh1 ~ Z1[,-1], tau= QQ[jq])
    YQ_high_IS[(h + 1):(h + full_length/2), jq] <- as.vector(Z2 %*% coef(b_high)) 
    
    # Initialize a vector for errors
    E_i <- rep(NA, length(Yh2))
    
    # Calculate errors for each point in the test set I2
    for (i in 1:length(E_i)) {
      E_i[i] <- max(YQ_low_IS[h + i, jq] - Yh2[i], Yh2[i] - YQ_high_IS[h + i, jq])
    }
    
    # Compute Q(QQ[jq])(E, I2) N.B 1 - ?? = QQ[jq]
    quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))
    
    YQ_IS[(h + 1):length(Yh), jq] <- as.vector(Z[1:(length(Yh) - h),] %*% coef(b_high)) + quantile_E
    
    
    
    
    
    #------ Conformalized Quantile regression with GDP only
    
    #b_lowGDPonly <- rq(Yh1 ~ ZGDPonly1[,-1], tau= 0.001) #Train on I1
    #YQ_lowGDPonly_IS[(h + 1):(h + full_length/2), jq] <- as.vector(ZGDPonly2 %*% coef(b_lowGDPonly)) #Evaluate on I2
    
    YQ_lowGDPonly_IS[(h + 1):(h + full_length/2), jq] <- -Inf 
    
    b_highGDPonly <- rq(Yh1 ~ ZGDPonly1[,-1], tau= QQ[jq])
    YQ_highGDPonly_IS[(h + 1):(h + full_length/2), jq] <- as.vector(ZGDPonly2 %*% coef(b_highGDPonly)) 
    
    # Initialize a vector for errors
    E_i <- rep(NA, length(Yh2))
    
    # Calculate errors for each point in the test set I2
    for (i in 1:length(E_i)) {
      E_i[i] <- max(YQ_lowGDPonly_IS[h + i, jq] - Yh2[i], Yh2[i] - YQ_highGDPonly_IS[h + i, jq])
    }
    
    # Compute Q(QQ[jq])(E, I2) N.B 1 - ?? = QQ[jq]
    quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))
    
    YQGDPonly_IS[(h + 1):length(Yh), jq] <- as.vector(ZGDPonly[1:(length(Yh) - h),] %*% coef(b_highGDPonly)) + quantile_E
    
    
    
    
    #------ Unconditional quantiles (conformalized quantile regression on constant)
    
    #bunc_low <- rq(Yh1 ~ 1, tau=0.001)
    #YQunclow_IS[(h + 1):(h + full_length/2), jq] <- rep(coef(bunc_low), length(Time) - h)
    YQunclow_IS[(h + 1):(h + full_length/2), jq] <- - Inf 
    
    bunc_high <- rq(Yh1 ~ 1, tau=QQ[jq])
    YQunchigh_IS[(h + 1):(h + full_length/2), jq] <- rep(coef(bunc_high), length((h + 1):(h + full_length/2)))
    
    # Initialize a vector for errors
    E_i <- rep(NA, length(Yh2))
    
    # Calculate errors for each point in the test set I2
    for (i in 1:length(E_i)) {
      E_i[i] <- max(YQunclow_IS[h + i, jq] - Yh2[i], Yh2[i] - YQunchigh_IS[h + i, jq])
    }
    
    # Compute Q(QQ[jq])(E, I2) N.B 1 - ?? = QQ[jq]
    quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))
    
    YQunc_IS[(h + 1):length(Yh), jq] <- rep(coef(bunc_high), length(Time) - h) + quantile_E
  }
  
  
  #---------    %% Fit skewed-t distribution for in-sample unconditional quantiles
  
  
    # This part is identical to the original paper
    
    QuantilesInterpolation_env <- new.env()
    source("QuantilesInterpolation.r",local = QuantilesInterpolation_env)
    

    #METODO BASE
    qqTarg <- YQunc_IS[nrow(YQunc_IS), ]
    list_params <- QuantilesInterpolation_env$QuantilesInterpolation(qqTarg, QQ) 
    lc <- list_params$lc #4.77
    sc <- list_params$sc #2.5
    sh <- list_params$sh #-8.62
    df <- list_params$df #30
    
    
    
    QuantilesInterpolation_env1 <- new.env()
    source("QuantilesInterpolationfaster1.r",local = QuantilesInterpolation_env1)
    
    #PARALLELIZZO IL CICLO SUI 30 GRADI DI LIBERTA': 3 VOLTE PIU' VELOCE
    qqTarg <- YQunc_IS[nrow(YQunc_IS), ]
    list_params1 <- QuantilesInterpolation_env1$QuantilesInterpolationfaster1(qqTarg, QQ) 
    lc <- list_params1$lc #4.77
    sc <- list_params1$sc #2.5
    sh <- list_params1$sh #-8.6
    df <- list_params1$df #30
    
    
    
    
    
    
    
    
    QuantilesInterpolation_env3 <- new.env()
    source("QuantilesInterpolationfaster3.r",local = QuantilesInterpolation_env3)
    
    qqTarg <- YQunc_IS[nrow(YQunc_IS), ]
    list_params3 <- QuantilesInterpolation_env3$QuantilesInterpolation(qqTarg, QQ) 
    lc <- list_params3$lc #4.77
    sc <- list_params3$sc #2.5
    sh <- list_params3$sh #-8.6
    df <- list_params3$df #30
    
    
    
    library(microbenchmark)
    results <- microbenchmark(
      QuantilesInterpolation_env$QuantilesInterpolation(qqTarg, QQ),
      QuantilesInterpolation_env3$QuantilesInterpolation(qqTarg, QQ),
      times = 1 # Number of times each function is executed
    )    
    summary(results)
    