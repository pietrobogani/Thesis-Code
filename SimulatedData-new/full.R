#' Distributional conformal prediction intervals.
#'
#' Compute prediction intervals using distributional conformal inference.
#'
#' @param x Matrix of features, of dimension (say) n x p.
#' @param y Vector of responses, of length (say) n.
#' @param x0 Matrix of features, each row being a point at which we want to
#'   form a prediction interval, of dimension (say) n0 x p.
#' @param train.fun A function to perform model training, i.e., to produce an
#'   estimator of Q(\tau|x), the conditional quantiles of the response variable
#'   Y given features X. Its input arguments should be x: matrix of features,
#'   y: vector of responses, \tau list of orders of the quantile to be
#'   estimated.
#' @param predict.fun A function to perform prediction for the conditional
#'   quantiles of the responses at new feature values. Its input arguments
#'   should be out: output produced by train.fun, and newx: feature values at
#'   which we want to make predictions.
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param num.grid.pts Number of grid points used when forming the conformal
#'   intervals (each grid point is a trial value for the interval). Default is
#'   100.
#'@param num.tau.pts Number of grid points for the discretisation of \tau to be
#'   used in the estimation of the conditional quantiles
#'@param grid.factor Expansion factor used to define the grid for the conformal
#'   intervals, i.e., the grid points are taken to be equally spaced in between
#'   -grid.factor*max(abs(y)) and grid.factor*max(abs(y)). Default is 1.25. In
#'   this case (and with exchangeable data, thus unity weights) the restriction
#'   of the trial values to this range costs at most 1/(n+1) in coverage. See
#'   details below.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#'
#' @return A list with the following components: pred, lo, up, fit. The first
#'   three are matrices of dimension n0 x m, and the last is a matrix of
#'   dimension n x m. Recall that n0 is the number of rows of x0, and m is the
#'   number of tuning parameter values internal to predict.fun. In a sense, each
#'   of the m columns really corresponds to a different prediction function;
#'   see details below. Hence, the rows of the matrices pred, lo, up give
#'   the predicted value, and lower and upper confidence limits (from conformal
#'   inference), respectively, for the response at the n0 points given in
#'   x0. The rows of fit give the fitted values for the n points given in x.
#'
#'
#'
#'#TO BE EDITED!
#' @details For concreteness, suppose that we want to use the predictions from
#'   forward stepwise regression at steps 1 through 5 in the path. In this case,
#'   there are m = 5 internal tuning parameter values to predict.fun, in the
#'   notation used above, and each of the returned matrices pred, lo, and up will
#'   have 5 rows (one corresponding to each step of the forward stepwise path).
#'   The code is structured in this way so that we may defined a single pair of
#'   functions train.fun and predict.fun, over a set of m = 5 tuning parameter
#'   values, instead of calling the conformal function separately m = 5 times.
#'
#'   The third arugment to train.fun, as explained above, is the output produced
#'   by a previous call to train.fun, but importantly, at the \emph{same}
#'   features x. The function train.fun may (optionally) leverage this returned
#'   output for efficiency purposes. Here is how it is used, in this function:
#'   when successively querying several trial values for the prediction
#'   interval, denoted, say, y0[j], j = 1,2,3,..., we set the out argument in
#'   train.fun to be the result of calling train.fun at the previous trial value
#'   y0[j-1]. Note of course that the function train.fun can also choose to
#'   ignore the returned output out, and most default training functions made
#'   available in this package will do so, for simplicity. An exception is the
#'   training function produced by \code{\link{lm.funs}}. This will use this
#'   output efficiently: the first time it trains by regressing onto a matrix of
#'   features x, it computes an appropriate Cholesky factorization; each
#'   successive time it is asked to train by regressing onto the same matrix x,
#'   it simply uses this Cholesky factorization (rather than recomputing one).
##   TODO implement and mention this for other functions too?
#'   The analogous explanation and discussion here applies to the out argument
#'   used by mad.train.fun.
#'
#'   If the data (training and test) are assumed to be exchangeable, the basic
#'   assumption underlying conformal prediction, then the probability that a new
#'   response value will lie outside of [-max(abs(y)), max(abs(y))], where y is
#'   the vector of training responses, is 1/(n+1).  Thus the restriction of the
#'   trials values to [-grid.factor*max(abs(y)), grid.factor*max(abs(y))], for
#'   all choices grid.factor >= 1, will lead to a loss in coverage of at most
#'   1/(n+1). This was also noted in "Trimmed Conformal Prediction for
#'   High-Dimensional Models" by Chen, Wang, Ha, Barber (2016) (who use this
#'   basic fact as motivation for proposing more refined trimming methods).
#'
#' @seealso \code{\link{dist.conformal.pred.split}},
#' EDIT!!!!
#' @references See "Algorithmic Learning in a Random World" by Vovk, Gammerman,
#'   Shafer (2005) as the definitive reference for conformal prediction; see
#'   also "" by Lei,
#'   G'Sell, Rinaldo, Tibshirani, Wasserman (2018) for another description; and
#'   "Conformal Prediction Under Covariate Shift" by Barber, Candes, Ramdas,
#'   Tibshirani (2019) for the weighted extension.
#' @example examples/ex.conformal.pred.R EDITTTT!!!!
#' @export dist.conformal.pred

dist.conformal.pred = function(x, y, x0, train.fun, predict.fun, alpha=0.1,
  num.grid.pts=500, num.grid.tau=100, grid.factor=1.25, verbose=FALSE) {

  # Set up data
  x = cbind(1,as.matrix(x))
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = cbind(1,matrix(x0,ncol=p-1))
  n0 = nrow(x0)
  paral = requireNamespace("future.apply", quietly=TRUE)

  # Check input arguments EDIT!
  #check.args(x=x,y=y,x0=x0,alpha=alpha,train.fun=train.fun,
  #          predict.fun=predict.fun,mad.train.fun=mad.train.fun,
  #           mad.predict.fun=mad.predict.fun)

  if (length(num.grid.pts) != 1 || !is.numeric(num.grid.pts)
      || num.grid.pts <= 1 || num.grid.pts >= 1000
      || round(num.grid.pts) != num.grid.pts) {
    stop("num.grid.pts must be an integer between 1 and 1000")
  }
  #check.pos.num(grid.factor)


  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }

  if (verbose) cat(sprintf("%sInitial training on full data set ...\n",txt))

  # Train, fit, and predict on full data set
  out = train.fun(x,y,0.5) #QR di livello 0.5. In output ho i parametri
  fit = matrix(predict.fun(out,x),nrow=n)
  pred = matrix(predict.fun(out,x0),nrow=n0) #Predizione del quantile 0.5 per i valori di x0 (nel mio caso forse i valori OOS....)


  # Trial values for y, empty lo, up matrices to fill
  ymax = max(abs(y))
  yvals = seq(-grid.factor*ymax, grid.factor*ymax,length=num.grid.pts)
  pvals = matrix(0,num.grid.pts)
  xx = rbind(x,rep(0,p))
  lo = up = matrix(0,n0)

  # Discretisation of conditional CDF
  tauvals = seq(0.001,0.999,length=num.grid.tau)



  for (i in 1:n0) { #cycle on all the points I want to get predictions for
    if (verbose) {
      cat(sprintf("\r%sProcessing prediction point %i (of %i) ...",txt,i,n0))
      flush.console()
    }

    xx[n+1,] = x0[i,]

    # Refit for each point in yvals, compute conformal p-value
    pvals = future.apply::future_sapply (1:num.grid.pts, function (j) {
      if (verbose) {
        cat(sprintf("\r%sProcessing grid point %i (of %i) ...",txt,j,num.grid.pts))
        flush.console()
      }
      yy = c(y,yvals[j])

      out = train.fun(xx,yy,tauvals)
      Q.yx = t(apply(predict.fun(out,xx),1,sort))
      u.hat = rowMeans(Q.yx <= matrix(yy,length(yy),length(tauvals)))
      #ncm=abs(u.hat-0.5)
      ncm <- u.hat             #CAMBIATO IO!!!!!!!!!!
      return(sum(ncm>=ncm[n+1])/(n+1))
    })

    grid_ok = yvals[pvals>alpha]
    lo[i,] = min(grid_ok)
    up[i,] = max(grid_ok)

}
  if (verbose) cat("\n")

  # Remove parallel structure
  if(paral){
    ## To avoid CRAN check errors
    ## R CMD check: make sure any open connections are closed afterward
    future::plan(future::sequential)
  }

  return(list(pred=pred,lo=lo,up=up,fit=fit))
}
