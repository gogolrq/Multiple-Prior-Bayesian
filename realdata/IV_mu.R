rm(list = ls())
library(LaplacesDemon)
library(ivpack)
library(Rmpfr)
set.seed(1)

setwd("E:/Dropbox/Multiple Prior Bayesian/code/realdata")
df=read.csv("mroz.csv")
index=which(df$lwage!=".")
data=df[index,]
data$lwage=as.numeric(data$lwage)
exp1 <- mpfr(exp(1), precBits=1024)
y=data$lwage
z=data$fatheduc
#x=cbind(data$educ, data$age, data$kidsge6, data$kidslt6, data$exper, data$expersq)
x=cbind(data$educ, data$kidsge6, data$kidslt6, data$exper, data$expersq)

v=cbind(1,z, x[,-1])
d<-cbind(y,x)

n=dim(data)[1]
p=5
J=p+1
li_reg2<-function(pars)
{
  data=d
  a<-pars[1] #intercept
  b<-pars[-1] #slope
  #print(a)
  #print(b)
  g0=data[,1]-a - data[,-1]%*%b

  g=v
  for (i in 1:(dim(v)[2])){
    g[,i]=g0*g[,i]
  }
  w=t(g)%*%g/n
  
  log_likelihood<- 0.5*colMeans(g)%*%solve(w)%*%colMeans(g)*(-n)
  # prior<- prior_reg(pars)
  return(log_likelihood)
}

srange=5
sigma=2
rho=function(theta,lambda=0){
  s=dnorm(theta[1],lambda, sigma)
  for (i in 2:6){
    s=s*dnorm(theta[i],lambda, sigma)
  }
  
  
  for (i in 1:6){
    s=s/dunif(theta[i],-srange,srange)
  }

  return(s)
}

r.sim=50000
mu.sim=runif(r.sim*(p+1),-srange,srange)
mu.sim=matrix(mu.sim,ncol=(p+1))
index=apply(mu.sim,1,li_reg2)
exp1power=exp1^index

lambda.sim = seq(-10, 10, by=1)


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
Iv.sim=as.numeric(Iv.sim)
print(Iv.sim)



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
Iv.sim=as.numeric(Iv.sim)
print(Iv.sim)




index=order(Iv.sim,decreasing = TRUE)
Iv.sim=Iv.sim[index]
lambda.sim=lambda.sim[index]


J=p+1
esti=c()

lambda.sim.new=lambda.sim
Iv.sim.new=Iv.sim
for (i in 1:length(lambda.sim.new)){
  lambda=lambda.sim.new[i]
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
    LL <- 0.5*colMeans(g)%*%solve(w)%*%colMeans(g)*(-n)
    LP <- LL + sum(beta.prior)
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(Data$y), Data$y, 1), parm=parm)
    return(Modelout)
  }
  
  Initial.Values <- rep(0,J)
  
  
  fit <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
                       Covar=NULL, Iterations=3500, Status=1000, Thinning=1, 
                       Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))
  
  print(fit)
  
  chain=fit$Posterior2
  thetahat=colMeans(chain)
  
  omegahat=t(chain)%*%chain/dim(chain)[1]-thetahat%*%t(thetahat)
  
  std=sqrt(diag(n*omegahat))
  
  std=std/sqrt(n)
  
  esti=rbind(esti,c(colMeans(chain),std))
  
  posterior=cbind(Iv.sim.new[i],fit$Posterior2)
  filename="posterior_mu.txt"
  write.table(posterior, file = filename, append = TRUE, quote = TRUE, sep = " ", col.names=FALSE)
}

result = data.frame(cbind(lambda.sim.new, Iv.sim.new, esti))
col.names=c("lambda", "Iv", paste("b", seq(0,p), sep=""), paste("sd", seq(0,p), sep=""))
colnames(result) = col.names



