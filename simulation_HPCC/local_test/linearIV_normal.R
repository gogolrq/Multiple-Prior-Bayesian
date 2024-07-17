rm(list = ls())
library(LaplacesDemon)
library(Rmpfr)
options(width=200)
args <- commandArgs(TRUE)
parameters = as.numeric(args)
Sys.sleep(0.1)
print(parameters)

exp1 <- mpfr(exp(1), precBits=1024)


n=4000#parameters[1]
sigma=1 #parameters[2]
seed=16#parameters[3]


p=3
set.seed(seed)

b0=c(-2.5,1.5,-1,2)
x<-runif(n*p,-2,4)
x=matrix(x,ncol=p)
#e=rlaplace(n, 0, 1)
e=rnorm(n, 0, 1)
z=runif(n, -2, 6)
x[,1]=z+2*e
error=e
y<-b0[1]+x%*%b0[-1]+error

v=cbind(1,z,x[,-1])
#v=cbind(1,x[,2]>0,x[,2],x[,3])
d<-cbind(y,x)

J=p+1


li_reg<-function(pars,data)
{
  # data=d
  a<-pars[1] #intercept
  b<-pars[-1] #slope
  g0=data[,1]-a - data[,-1]%*%b
  g=v
  for (i in 1:(dim(v)[2])){
    g[,i]=g0*g[,i]
  }
  #w=wm
  w=t(g)%*%g/n
  log_likelihood<- colMeans(g)%*%solve(w)%*%colMeans(g)*(-n)
  prior<- prior_reg(pars)
  return(log_likelihood + prior)
}




li_reg2<-function(pars)
{
  data=d
  a<-pars[1] #intercept
  b<-pars[-1] #slope
  g0=data[,1]-a - data[,-1]%*%b
  g=v
  for (i in 1:(dim(v)[2])){
    g[,i]=g0*g[,i]
  }
  w=t(g)%*%g/n
  
  log_likelihood<- colMeans(g)%*%solve(w)%*%colMeans(g)*(-n)
  # prior<- prior_reg(pars)
  return(log_likelihood)
}

rho=function(theta,lambda=0){
  s=dnorm(theta[1],lambda,sigma)
  s=s*dnorm(theta[2],lambda,sigma)
  s=s*dnorm(theta[3],lambda,sigma)
  s=s*dnorm(theta[4],lambda,sigma)
  
  #s=s/dnorm(theta[1],0,sigma)
  #s=s/dnorm(theta[2],0,sigma)
  #s=s/dnorm(theta[3],0,sigma)
  #s=s/dnorm(theta[4],0,sigma)
  
  
  s=s/dunif(theta[1],-100,100)
  s=s/dunif(theta[2],-100,100)
  s=s/dunif(theta[3],-100,100)
  s=s/dunif(theta[4],-100,100)
  
  return(s)
}

r.sim=50000
#mu.sim=rnorm(r.sim*4,0,100)
mu.sim=runif(r.sim*4,-100,100)
mu.sim=matrix(mu.sim,ncol=4)
index=apply(mu.sim,1,li_reg2)
exp1power=exp1^index


lambda.sim = seq(-10, 10, length=50)

Iv.sim=c()
result=data.frame()
Iv.max=0
for (i in 1:length(lambda.sim)){
  lambda=lambda.sim[i]
  impf=apply(mu.sim,1,rho,lambda=lambda)
  Iv=mean(exp1power*impf)
  Iv=mpfr(Iv,8)
  temp=Iv
  if (length(Iv.sim)==0){
    Iv.sim=temp
  }else{
    Iv.sim=cbind(Iv.sim,temp)
  }
}
Iv.sim=Iv.sim/max(Iv.sim)
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
#lambda.sim.new=c(lambda.sim.new, 0)
#Iv.sim.new=c(Iv.sim.new, Iv.sim[which.min(Iv.sim)])



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
  
  MyData <- list(J=p+1, PGF=PGF, X=x, mon.names=mon.names,  parm.names=parm.names, pos.beta=pos.beta, y=y, v=v)
  
  Model <- function(parm, Data){
    
    beta <- parm[Data$pos.beta]
    
    
    
    beta.prior <- dnormv(beta, lambda, sigma, log=TRUE)
    
    data=d
    a<-beta[1] #intercept
    b<-beta[-1] #slope
    g0=Data$y-a - Data$X%*%b
    g=Data$v
    for (i in 1:(dim(v)[2])){
      g[,i]=g0*g[,i]
    }
    w=t(g)%*%g/n
    LL <- colMeans(g)%*%solve(w)%*%colMeans(g)*(-n)
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
  
  std=sqrt(diag(n*omegahat))
  
  esti=rbind(esti,c(colMeans(chain),std))
  

}

result = data.frame(cbind(omega.sim,lambda.sim.new, Iv.sim.new, esti))
colnames(result) = c("omega","lambda", "Iv", "a", "b1", "b2", "b3", "sda", "sd1", "sd2", "sd3")
result$n = n
result$seed = seed
result$sigma = sigma


x1hat=lm(x[,1]~z)$fit
summary(lm(y~z+x[,2]+x[,3]))