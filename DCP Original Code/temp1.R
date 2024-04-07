#rq.fit() does not exist, I try to use rq from quantreg package

library(quantreg)

rq.funs=function(){   #I fit several QR at levels inside tau_list. I return a matrix with all the coefficients and as columns different
                      #levels of alpha. max is (length(x)+1)*length(tau_list) Matrix
train.fun=function(x,y,tau_list){

  p=ncol(x)
  n_tau=length(tau_list)
  mat=sapply(tau_list,function(tau_element){rq.fit(x,y,tau=tau_element)$coefficients})
  colnames(mat)=tau_list
  return(mat)
}

predict.fun=function(out,x){

    x %*% out

  }

return(list(train.fun=train.fun, predict.fun=predict.fun))

}





