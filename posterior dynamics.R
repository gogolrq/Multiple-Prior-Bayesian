rm(list=ls())
library(ggplot2)
library(latex2exp)
set.seed(1)
setwd("E:/Dropbox/Multiple Prior Bayesian/code/realdata")
k=80
z=rnorm(k)
w=0.95
w.sim=c(1, 0.90, 0)
result=c()
for (n in 15:k){
  for (w in w.sim){
    zbar=mean(z[1:n])
    if (w>0){
      ub=zbar+sqrt(2*(n+1)*log(1/w)/n)
      lb=zbar-sqrt(2*(n+1)*log(1/w)/n)      
    }else{
      ub=3
      lb=-3
    }

    
    
    
    ub=(ub+n*zbar)/(n+1)
    lb=(lb+n*zbar)/(n+1)
    temp=c(n,w, ub,lb)
    result=rbind(result,temp)
  }
}
result=data.frame(result)
colnames(result)=c("n","w","ub","lb")

#plot(result$n,result$lb,type="l")
#lines(result$n,result$ub)

data=result
data$w=factor(data$w, levels=w.sim)
filename="../../Figures/pdynamics.eps"
ps.options(horizontal = F)
ps.options(height=5, width=6)
postscript(filename)
p1=ggplot(data, aes(x=n))+
  geom_line(aes(y=ub,color=w, linetype=w),linewidth=1)+
  geom_line(aes(y=lb,color=w, linetype=w),linewidth=1)+
  labs("test")+
  xlab("n")+
  ylab("")+
 # ylab(TeX("$\\theta$"))+
  scale_color_manual(values=c("black", "red", "purple"))+
  scale_linetype_manual(values=c("solid","longdash","twodash"))+
  labs(color = TeX("$\\omega$    "))+
  labs(linetype = TeX("$\\omega$    ") )+
  #theme(legend.position = "bottom")+
  theme(legend.title.align=0.5,legend.position = c(0.9, 0.8),legend.background = element_rect(fill="transparent",
                                                                       size=0.3, linetype="solid", 
                                                                       colour="black"))+
  theme(text = element_text(size=15))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=15),axis.text = element_text(face="bold",size=13),legend.text =element_text(size=15) )
  #theme(plot.title = element_text(hjust = 0.5,face="bold",size=20),axis.text = element_text(face="bold",size=15),legend.text =element_text(size=15) )
  
print(p1)
dev.off()


n.sim=c(80, 50 , 30)
result=c()
for (n in n.sim){
  w.sim=seq(0.01, 1, by=0.01)
  ub=mean(z[1:n])+sqrt(2*(n+1)*log(1/w.sim)/n)
  lb=mean(z[1:n])-sqrt(2*(n+1)*log(1/w.sim)/n)
  ub=(ub+sum(z[1:n]))/(n+1)
  lb=(lb+sum(z[1:n]))/(n+1)
  temp=cbind(n,w.sim, ub,lb)
  result=rbind(result,temp)
}

result=data.frame(result)
colnames(result)=c("n","w","ub","lb")
data2=result
data2$n=factor(data2$n, levels=n.sim)
filename="../../Figures/pdynamics2.eps"
ps.options(horizontal = F)
ps.options(height=5, width=6)
postscript(filename)
p2=ggplot(data2, aes(x=w))+
  geom_line(aes(y=ub,color=n, linetype=n),linewidth=1)+
  geom_line(aes(y=lb,color=n, linetype=n),linewidth=1)+
  labs("test")+
  xlab(TeX("$\\omega$"))+
  ylab("")+
  scale_color_manual(values=c("black", "red", "purple"))+
  scale_linetype_manual(values=c("solid","longdash","twodash"))+
  labs(color = "n")+
  labs(linetype = "n" )+
  theme(legend.title.align=0.5,legend.position = c(0.9, 0.8),legend.background = element_rect(fill="transparent",
                                                                       size=0.3, linetype="solid", 
                                                                       colour="black"))+
  theme(text = element_text(size=15))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=15),axis.text = element_text(face="bold",size=13),legend.text =element_text(size=15) )

print(p2)
dev.off()