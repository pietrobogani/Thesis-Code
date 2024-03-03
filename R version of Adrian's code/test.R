
# Create a cluster
cl <- makeCluster(detectCores() - 1) # Use one less than the total number of cores

# Export necessary objects and functions to the cluster if needed
clusterExport(cl, c("QuantilesInterpolation_env", "YQ_OOS", "YQGDPonly_OOS", "YQunc_OOS", "QQ"))

# List of tasks as strings to be evaluated in the global environment
tasks <- c("QuantilesInterpolation_env$QuantilesInterpolation(YQ_OOS, QQ)",
           "QuantilesInterpolation_env$QuantilesInterpolation(YQGDPonly_OOS, QQ)",
           "QuantilesInterpolation_env$QuantilesInterpolation(YQunc_OOS, QQ)")

# Run tasks in parallel using clusterEvalQ to evaluate expressions
results <- parLapply(cl, tasks, function(task) eval(parse(text=task)))

# Stop the cluster
stopCluster(cl)

# Unpack results
params <- results[[1]]
params_GDPonly <- results[[2]]
params_unc <- results[[3]]

PST_OOS[jt + h, ] <- dst(YY, params$lc, params$sc, params$sh, params$df)
QST_OOS[jt + h, ] <- qst(QQ, params$lc, params$sc, params$sh, params$df)
CST_OOS[jt + h, ] <- pst(YY, params$lc, params$sc, params$sh, params$df)
STpar_OOS[jt + h, ] <- c(params$lc, params$sc, params$sh, params$df)
ScoreST_OOS[jt + h] <- dst(YhRealized, params$lc, params$sc, params$sh, params$df)
PitST_OOS[jt + h] <- pst(YhRealized, params$lc, params$sc, params$sh, params$df) # is the probability to observe a value < of YhRealized in this distribution 

# Fit skewed-t distribution for quantile regression with GDP only, out-of-sample
PSTGDPonly_OOS[jt + h, ] <- dst(YY, params_GDPonly$lc, params_GDPonly$sc, params_GDPonly$sh, params_GDPonly$df)
QSTGDPonly_OOS[jt + h, ] <- qst(QQ, params_GDPonly$lc, params_GDPonly$sc, params_GDPonly$sh, params_GDPonly$df)
CSTGDPonly_OOS[jt + h, ] <- pst(YY, params_GDPonly$lc, params_GDPonly$sc, params_GDPonly$sh, params_GDPonly$df)
STparGDPonly_OOS[jt + h, ] <- c(params_GDPonly$lc, params_GDPonly$sc, params_GDPonly$sh, params_GDPonly$df)
ScoreSTGDPonly_OOS[jt + h] <- dst(YhRealized, params_GDPonly$lc, params_GDPonly$sc, params_GDPonly$sh, params_GDPonly$df)
PitSTGDPonly_OOS[jt + h] <- pst(YhRealized, params_GDPonly$lc, params_GDPonly$sc, params_GDPonly$sh, params_GDPonly$df) # is the probability to observe a value < of YhRealized in this distribution 

# Fit skewed t-distribution for unconditional quantiles, out-of-sample
PSTunc_OOS[jt + h, ] <- dst(YY, params_unc$lc, params_unc$sc, params_unc$sh, params_unc$df)
QSTunc_OOS[jt + h, ] <- qst(QQ, params_unc$lc, params_unc$sc, params_unc$sh, params_unc$df)
CSTunc_OOS[jt + h, ] <- pst(YY, params_unc$lc, params_unc$sc, params_unc$sh, params_unc$df)
STparunc_OOS[jt + h, ] <- c(params_unc$lc, params_unc$sc, params_unc$sh, params_unc$df)
ScoreSTunc_OOS[jt + h] <- dst(YhRealized, params_unc$lc, params_unc$sc, params_unc$sh, params_unc$df)
PitSTunc_OOS[jt + h] <- pst(YhRealized, params_unc$lc, params_unc$sc, params_unc$sh, params_unc$df) # is the probability to observe a value < of YhRealized in this distribution 
