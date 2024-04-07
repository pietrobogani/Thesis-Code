n=100
x=runif(n)
u=rnorm(n)
y=x+x*u
plot(x,y)



funs=rq.funs()
predict.fun=funs$predict.fun
train.fun=funs$train.fun



x0=c(0.1,0.8)


test=dist.conformal.pred(x,y,x0,train.fun = train.fun,predict.fun = predict.fun, verbose = T)
