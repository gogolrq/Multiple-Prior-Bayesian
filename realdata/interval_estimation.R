rm(list = ls())
library(ivpack)
library(ggplot2)
library(latex2exp)


setwd("E:/Dropbox/Multiple Prior Bayesian/code/realdata")
df=read.table("posterior.txt")
df=df[,-1]
colnames(df)=c("omega", "b0", "b1", "b2", "b3", "b4", "b5")
omega.sim=unique(df$omega)
vnames=c("educ","kidsge6", "kidslt6","exper","expersq")




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



num=dim(esti)[1]
esti=esti[num:1,]

result=c()

for (j in 1:5){
  esti$b=esti[,j+2]
  for (i in 1:num){
    temp=esti[i:num,]
    print(temp)
    lb=min(temp$b)
    ub=max(temp$b)
    result=rbind(result,c(j,esti[i,1],lb,ub))
  }

}
result=data.frame(result)
colnames(result)=c("j","w","lb","ub")
result$j=factor(result$j)

for (j in 1:5){
  subdf=result[which(result$j==j),]
  subdf=subdf[-c(1:6),]
  filename=paste("../../Figures/IE",j,".eps",sep="")
  ps.options(horizontal = F)
  ps.options(height=3, width=3)
  postscript(filename)
  par(mar=c(2,2,3,2))
  p1=ggplot(subdf, aes(x=w))+
    geom_line(aes(y=ub),linewidth=0.6)+
    geom_line(aes(y=lb),linewidth=0.6)+
    geom_point(aes(y=ub), size = 1.5, color="black")+
    geom_point(aes(y=lb), size = 1.5, color="black")+
    theme(text = element_text(size=15))+
    xlab(TeX("$\\omega$"))+
    ylab("")+
    ggtitle(vnames[j])+
    theme(plot.title = element_text(hjust = 0.5,face="bold",size=10),axis.text = element_text(face="bold",size=9),legend.text =element_text(size=15) )
  print(p1)
  dev.off()
}

