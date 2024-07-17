rm(list = ls())
setwd("/home/ruiqliu/MPB/v2/linear_normal")
library(LaplacesDemon)
library(Rmpfr)
options(width=200)
args <- commandArgs(TRUE)
parameters = as.numeric(args)
Sys.sleep(0.1)
print(parameters)

exp1 <- mpfr(exp(1), precBits=1024)


n=parameters[1]
sigma=parameters[2]
seed=parameters[3]


set.seed(seed)

p = 3
b0 = c(-1.5, 3, 1)
a0 = -2


x <- runif(n * p, -3, 4)
x <- matrix(x, ncol = p)
y <- a0 + x %*% b0 + rnorm(n, 0, 2)
d <- cbind(y, x)


li_reg <- function(pars) #loglkh+prior
{
  data = d
  a <- pars[1] #intercept
  b <- pars[-1] #slope
  log_likelihood <- sum((data[, 1] - a - data[, -1] %*% b) ^ 2) * -1
  prior <- prior_reg(pars)
  return(log_likelihood + prior)
}





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

prior_reg <- function(pars)
{
  s = 0
  a <- pars[1] #intercept
  prior_a <-
    dunif(a, lambda, sigma) ## non-informative (flat) priors on all
  if (prior_a == 0) {
    prior_a = -1000000000
  } else{
    prior_a = log(prior_a)
  }
  s = s + prior_a
  for (i in 2:(p + 1)) {
    b <- pars[i] #slope
    prior_b <- dunif(b, lambda, sigma) ## parameters.
    if (prior_b == 0) {
      prior_b = -1000000000
    } else{
      prior_b = log(prior_b)
    }
    s = s + prior_b
  }
  #print(s)
  #print(lambda)
  return(s)
}



r.sim = 50000
mu.sim=runif(r.sim*4,-srange,srange)
mu.sim = matrix(mu.sim, ncol = 4)
index = apply(mu.sim, 1, li_reg2)
exp1power = exp1 ^ index


lambda.sim = seq(-10, 10, by=1)

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
    
    
    
    beta.prior <- dnormv(beta, lambda, sigma, log=TRUE)
    
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

result = data.frame(cbind(omega.sim,lambda.sim.new, Iv.sim.new, esti))
colnames(result) = c("omega","lambda", "Iv", "a", "b1", "b2", "b3", "sda", "sd1", "sd2", "sd3")
result$n = n
result$seed = seed
result$sigma = sigma



print(result)
filename="data/result.txt"
write.table(result, file = filename, append = TRUE, quote = TRUE, sep = " ", col.names=FALSE)