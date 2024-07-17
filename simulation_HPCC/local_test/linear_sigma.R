rm(list = ls())
library(LaplacesDemon)
library(Rmpfr)
options(width=200)
args <- commandArgs(TRUE)
parameters = as.numeric(args)
Sys.sleep(0.1)
print(parameters)

exp1 <- mpfr(exp(1), precBits=1024)


n=100#parameters[1]
seed=16#parameters[2]
set.seed(seed)

p = 3
b0 = c(-1.5, 3, 1)
a0 = -2


x <- runif(n * p, -3, 4)
x <- matrix(x, ncol = p)
y <- a0 + x %*% b0 + rnorm(n, 0, 2)
d <- cbind(y, x)




li_reg2 <- function(pars) #loglkh
{
  data = d
  a <- pars[1] #intercept
  b <- pars[-1] #slope
  log_likelihood <- sum((data[, 1] - a - data[, -1] %*% b) ^ 2) * -1
  
  return(log_likelihood)
}




srange=5
rho=function(theta,lambda=0){
  s=dnorm(theta[1],lambda,sigma)
  s=s*dnorm(theta[2],lambda,sigma)
  s=s*dnorm(theta[3],lambda,sigma)
  s=s*dnorm(theta[4],lambda,sigma)
  
  #s=s/dnorm(theta[1],0,sigma)
  #s=s/dnorm(theta[2],0,sigma)
  #s=s/dnorm(theta[3],0,sigma)
  #s=s/dnorm(theta[4],0,sigma)
  
  
  s=s/dunif(theta[1],-srange,srange)
  s=s/dunif(theta[2],-srange,srange)
  s=s/dunif(theta[3],-srange,srange)
  s=s/dunif(theta[4],-srange,srange)
  
  return(s)
}

srange=5
rho=function(theta,lambda=1){
  s=dnorm(theta[1],0,lambda)
  s=s*dnorm(theta[2],0,lambda)
  s=s*dnorm(theta[3],0,lambda)
  s=s*dnorm(theta[4],0,lambda)
  
  #s=s/dnorm(theta[1],0,sigma)
  #s=s/dnorm(theta[2],0,sigma)
  #s=s/dnorm(theta[3],0,sigma)
  #s=s/dnorm(theta[4],0,sigma)
  
  
  s=s/dunif(theta[1],-srange,srange)
  s=s/dunif(theta[2],-srange,srange)
  s=s/dunif(theta[3],-srange,srange)
  s=s/dunif(theta[4],-srange,srange)
  
  return(s)
}

r.sim=50000
#mu.sim=rnorm(r.sim*4,0,100)
mu.sim=runif(r.sim*4,-srange,srange)
mu.sim=matrix(mu.sim,ncol=4)
index=apply(mu.sim,1,li_reg2)
exp1power=exp1^index


lambda.sim = seq(1, 10, by =0.5)

Iv.sim=c()
result=data.frame()
Iv.max=0
for (i in 1:length(lambda.sim)){
  lambda=lambda.sim[i]
  impf=apply(mu.sim,1,rho,lambda=lambda)
  Iv=mean(exp1power*impf)
  Iv=mpfr(Iv,precBits=1024)
  temp=Iv
  if (length(Iv.sim)==0){
    Iv.sim=temp
  }else{
    Iv.sim=cbind(Iv.sim,temp)
  }
}
Iv.sim=Iv.sim/max(Iv.sim)
Iv.sim=as.numeric(Iv.sim)
print(Iv.sim)


omega.sim=c(1, 0.95, 0.9,0)
lambda.sim.new=c()
Iv.sim.new=c()
for (omega in omega.sim){
  lambda.index=which(Iv.sim>=omega)
  lambda.index=lambda.index[which.min(Iv.sim[lambda.index])]
  lambda.sim[lambda.index]
  lambda.sim.new=c(lambda.sim.new, lambda.sim[lambda.index])
  Iv.sim.new=c(Iv.sim.new, Iv.sim[lambda.index])
}




J=p+1
esti=c()
for (lambda in lambda.sim.new){
  mon.names <- "LP"
  parm.names <- as.parm.names(list(beta=rep(0,J)))
  pos.beta <- grep("beta", parm.names)
  #pos.sigma <- grep("sigma", parm.names)
  
  PGF <- function(Data) { 
    beta <- rnorm(Data$J) 
    return(beta) 
  }
  
  MyData <- list(J=p+1, PGF=PGF, X=x, mon.names=mon.names,  parm.names=parm.names, pos.beta=pos.beta, y=y)
  
  Model <- function(parm, Data){
    
    beta <- parm[Data$pos.beta]
    
    
    
    beta.prior <- dnormv(beta, 0, lambda, log=TRUE)
    
    data=d
    a<-beta[1] #intercept
    b<-beta[-1] #slope
    fvalue=plogis(a+data[,-1]%*%b)
    
    LL=sum((data[, 1] - a - data[, -1] %*% b) ^ 2) * -1
    LP <- LL + sum(beta.prior)
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(Data$y), Data$y, 1), parm=parm)
    return(Modelout)
  }
  
  Initial.Values <- c(rep(0,J))
  
  
  fit <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
                       Covar=NULL, Iterations=3000, Status=100, Thinning=1, 
                       Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))
  
  print(fit)
  
  chain=fit$Posterior2
  
  thetahat=colMeans(chain)
  omegahat=t(chain)%*%chain/dim(chain)[1]-thetahat%*%t(thetahat)
  tempx=cbind(1,x)
  for (i in 1:(dim(tempx)[2])){
    tempx[,i]=tempx[,i]*(y-cbind(1,x)%*%thetahat)*2
  }
  vhat=t(tempx)%*%tempx
  
  std=sqrt(diag(n*omegahat%*%vhat%*%omegahat))
  
  
  esti=rbind(esti,c(colMeans(chain),std))
  
  
}

result = data.frame(cbind(omega.sim, lambda.sim.new, Iv.sim.new, esti))
colnames(result) = c("omega","lambda", "Iv", "a", "b1", "b2", "b3", "sda", "sd1", "sd2", "sd3")
result$n = n
result$seed = seed


