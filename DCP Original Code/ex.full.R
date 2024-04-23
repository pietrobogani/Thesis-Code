source("temp1.R")
source("full.R")


#  PER FARLO FUNZIONARE HO DOVUTO USARE SOURCE DUE VOLTE E AGGIUNGERE LA LIBRARY quantreg in temp1

set.seed(1)
n=100
x=runif(n)
u=rnorm(n)
y=x+x*u
plot(x,y)



funs=rq.funs()
predict.fun=funs$predict.fun
train.fun=funs$train.fun



x0=x


test=dist.conformal.pred(x,y,x0,train.fun = train.fun, predict.fun = predict.fun, verbose = T)
points(c(0.1,0.8),test$lo[,1],col='red')
points(c(0.1,0.8),test$up[,1],col='green')







ytrial <- seq(-max(abs(y)), max(abs(y)), 0.01)
alpha<-0.1
x0 <- 0.8
pvalues <- fullDCP(y, x, x0, alpha, ytrial)
indices <- which(pvalues > alpha)
CI <- c(ytrial[indices[1]], ytrial[indices[length(indices)]])
print(CI)



ynew <- 0.8 + 0.8*u
mean(ynew <=CI[2])
mean(ynew <= test$up[2,1])
