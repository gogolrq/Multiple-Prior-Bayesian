rm(list = ls())
library(ivpack)
library(ggplot2)
library(latex2exp)


setwd("E:/Dropbox/Multiple Prior Bayesian/code/realdata")
df=read.table("posterior.txt")
df=df[,-1]
colnames(df)=c("omega", "b0", "b1", "b2", "b3", "b4", "b5")
omega.sim=unique(df$omega)

data=read.csv("mroz.csv")
index=which(data$lwage!=".")
data=data[index,]
data$lwage=as.numeric(data$lwage)
y=data$lwage
z=data$fatheduc
x=cbind(data$educ, data$age, data$kidsge6, data$kidslt6)
v=cbind(1,z, x[,-1])
d=cbind(y,x)
n=dim(data)[1]



lmod1=lm(lwage~educ+kidsge6+kidslt6+exper+expersq,data=data)
lmod2=ivreg(lwage~educ+kidsge6+kidslt6+exper+expersq|fatheduc+kidsge6+kidslt6+exper+expersq,data=data)
summary(lmod1)
summary(lmod2)



vnames=c("educ","kidsge6", "kidslt6","exper","expersq")
beta=1 #change beta from 1-4


subdf=df[which(df$omega==omega.sim[1]),] # 
subdf=subdf[-c(1:500),-1] 
chain=as.matrix(subdf)
thetahat=colMeans(chain)
omegahat=t(chain)%*%chain/dim(chain)[1]-thetahat%*%t(thetahat)
std=sqrt(diag(n*omegahat))
std=std/sqrt(n)

cbind(thetahat,std, pnorm(abs(thetahat/std),lower.tail = FALSE)*2)

esti=c()

for (omega in omega.sim){
  subdf=df[which(df$omega==omega),]
  subdf=subdf[-c(1:500),-1]
  chain=as.matrix(subdf)
  thetahat=colMeans(chain)
  esti=rbind(esti,c(omega,thetahat))
}
colnames(esti)[1]="omega"
esti=data.frame(esti)