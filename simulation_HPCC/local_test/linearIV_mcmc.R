rm(list = ls())
library(LaplacesDemon)
library(Rmpfr)
options(width=200)
args <- commandArgs(TRUE)
parameters = as.numeric(args)
Sys.sleep(0.1)
print(parameters)

exp1 <- mpfr(exp(1), precBits=1024)


n=1000#parameters[1]
sigma=20 #parameters[2]
seed=10#parameters[3]


p=3
set.seed(seed)

b0=c(-2.5,1.5,-1,2)
x<-runif(n*p,-2,4)
x=matrix(x,ncol=p)
#e=rlaplace(n, 0, 1)
e=rnorm(n, 0, 1)
z=runif(n, -2, 6)
z=x[,1]
#x[,1]=z+2*e
error=e
y<-b0[1]+x%*%b0[-1]+error

v=cbind(1,z,x[,-1])
#v=cbind(1,x[,2]>0,x[,2],x[,3])
d<-cbind(y,x)

J=p+1

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
  
  
  
  beta.prior <- dnormv(beta, 0, sigma, log=TRUE)
  
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

Initial.Values <-b0  #c(rep(0,J))


fit <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
                     Covar=NULL, Iterations=3000, Status=100, Thinning=1, 
                     Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))

print(fit)

chain=fit$Posterior2
thetahat=colMeans(chain)

omegahat=t(chain)%*%chain/dim(chain)[1]-thetahat%*%t(thetahat)

std=sqrt(diag(n*omegahat))


summary(lm(y~x))
#x1hat=lm(x[,1]~z)$fit
#summary(lm(y~z+x[,2]+x[,3]))