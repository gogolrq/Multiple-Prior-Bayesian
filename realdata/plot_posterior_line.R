rm(list = ls())
library(ivpack)
library(ggplot2)
library(latex2exp)
library(reshape2)

setwd("E:/Dropbox/Multiple Prior Bayesian/code/realdata")
df=read.table("posterior.txt")
df=df[,-1]
colnames(df)=c("omega", "b0", "b1", "b2", "b3", "b4", "b5")
omega.sim=unique(df$omega)





vnames=c("educ","kidsge6", "kidslt6","exper","expersq")
#w=1 #change from 1-3
#beta=1#change from 1-5

for (w in c(1,2,3)){
  subdf=df[which(df$omega==omega.sim[w]),]
  subdf=subdf[-c(1:500),]
  for (beta in c(1,2,3,4,5)){
    
    
    colMeans(subdf)
    subdf$index=1:dim(subdf)[1]
    #subdf$b=subdf[,beta+2]
    #subdf[,beta+2]=(subdf[,beta+2]-min(subdf[,beta+2]))/diff(range(subdf[,beta+2]))
    
    subdf[,beta+2]=(subdf[,beta+2]-mean(subdf[,beta+2]))/sd(subdf[,beta+2])/4+(5-beta)*2
    print(mean(subdf[,beta+2]))
  }
  subdf=subdf[,-2]
  subdf=melt(subdf, id.vars=c("index", "omega"))
  filename=paste("../../Figures/line_beta",w,".eps",sep="")
  ps.options(horizontal = F)
  ps.options(height=5, width=5.5)
  postscript(filename)
  par(mar=c(2,2,3,2))
  p1=ggplot(data=subdf, aes(x=index,y=value, group=variable))+
    geom_line(aes(color=variable), linewidth=1.2)+
    scale_color_manual(values=c("black", "red", "purple","cyan","chartreuse"), labels = vnames)+
    theme(text = element_text(size=15))+
    xlab("")+
    ylab("")+
    labs(color = "")+
    theme(axis.text.y=element_blank())+
    #ggtitle(vnames[beta])+ legend.position = "bottom",
    theme(legend.position = "top",plot.title = element_text(hjust = 0.5,face="bold",size=15),axis.text = element_text(face="bold",size=15),legend.text =element_text(size=15) )
  print(p1)
  dev.off()
}
