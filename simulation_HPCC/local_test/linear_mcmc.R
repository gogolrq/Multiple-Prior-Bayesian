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
sigma=1 #parameters[2]
seed=16#parameters[3]
set.seed(seed)

p = 3
b0 = c(-1.5, 3, 1)
a0 = -2


x <- runif(n * p, -3, 4)
x <- matrix(x, ncol = p)
y <- a0 + x %*% b0 + rnorm(n, 0, 1)
d <- cbind(y, x)

J=p+1
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
  
  
  
  beta.prior <- dnormv(beta, 0, sigma, log=TRUE)
  
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
                     Covar=NULL, Iterations=3000, Status=100, Thinning=2, 
                     Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))

chain=fit$Posterior2

thetahat=colMeans(chain)
omegahat=t(chain)%*%chain/dim(chain)[1]-thetahat%*%t(thetahat)
tempx=cbind(1,x)
for (i in 1:(dim(tempx)[2])){
  tempx[,i]=tempx[,i]*(y-cbind(1,x)%*%thetahat)*2
}
vhat=t(tempx)%*%tempx

std=sqrt(diag(n*omegahat%*%vhat%*%omegahat))