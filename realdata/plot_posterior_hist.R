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
#w=1 #change from 1-3
#beta=1#change from 1-5

for (w in c(1,2,3)){
  for (beta in c(1,2,3,4,5)){
    subdf=df[which(df$omega==omega.sim[w]),]
    subdf=subdf[-c(1:500),]
    colMeans(subdf)
    subdf$index=1:dim(subdf)[1]
    subdf$b=subdf[,beta+2]
    
    # filename=paste("../../Figures/line_beta",w,beta,".eps",sep="")
    # ps.options(horizontal = F)
    # ps.options(height=4, width=6)
    # postscript(filename)
    # par(mar=c(2,2,3,2))
    # p1=ggplot(data=subdf, aes(x=index,y=b))+
    #   geom_line()+
    #   theme(text = element_text(size=15))+
    #   xlab("")+
    #   ylab("")+
    #   ggtitle(vnames[beta])+
    #   theme(plot.title = element_text(hjust = 0.5,face="bold",size=10),axis.text = element_text(face="bold",size=9),legend.text =element_text(size=15) )
    # print(p1)
    # dev.off()
    
    filename=paste("../../Figures/hist_beta",w,beta,".eps",sep="")
    ps.options(horizontal = F)
    ps.options(height=5, width=5.5)
    postscript(filename)
    par(mar=c(2,2,3,2))
    p2=ggplot(subdf, aes(x=b, y=after_stat(density)))+
      geom_histogram(color="darkblue", fill="lightblue")+
      geom_density(alpha=.2, fill="#FF6666")+
      theme(text = element_text(size=15))+
      #xlab(TeX(paste("$\\beta_",1,"$",sep="")))+
      xlab("")+
      ylab("")+
      ggtitle(vnames[beta])+
      theme(plot.title = element_text(hjust = 0.5,face="bold",size=15),axis.text = element_text(face="bold",size=15),legend.text =element_text(size=15) )
    print(p2)
    dev.off()
  }
}
