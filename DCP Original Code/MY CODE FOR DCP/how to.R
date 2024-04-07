#After training model1 and 2 on train data, I generate prediction intervals at level 95% for y_test using x_test.
#The model is perfectly calibrated if 95%   of values of y_test actually fall inside their own prediction interval.


#1) Train the 2 models In-Sample for lots of coverages = c(0.01,0.02,...,0.98,0.99)
#   For each coverage, calculate 2 regressions for each model: 0.01 => QQ = c(0.495,0.505)

#2) With OOS data, use regressors values (like NFCI) to produce confidence intervals of their respective OOS target variable.
#3) If the model is perfect calibrated, a coverage=0.95 interval will contain 95% of all the OOS values
