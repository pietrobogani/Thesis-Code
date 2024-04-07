
rq.funs=function(){
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





