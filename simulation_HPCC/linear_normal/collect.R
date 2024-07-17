rm(list = ls())
setwd("/home/ruiqliu/MPB/v2/linear_normal")
options(width = 200)
n.sim = c(100,300,500, 800, 1000)
output = c()

result = c() 

filename = paste("data/result.txt")
df = read.table(filename, sep = " ")
df=df[,-1]

#df=result
colnames(df)=c("omega","lambda","Iv","a","b1","b2","b3","sda", "sd1", "sd2", "sd3","n","seed","sigma")


b0=c(-2,-1.5, 3, 1)
error=(df$a-b0[1])^2+(df$b1-b0[2])^2+(df$b2-b0[3])^2+(df$b3-b0[4])^2
error=(df$b1-b0[2])^2+(df$b2-b0[3])^2+(df$b3-b0[4])^2
error=sqrt(error)


check1=b0[2]<df$b1+df$sd1*1.96/sqrt(df$n)
check2=b0[2]>df$b1-df$sd1*1.96/sqrt(df$n)
checkb1=check1*check2

check1=b0[3]<df$b2+df$sd2*1.96/sqrt(df$n)
check2=b0[3]>df$b2-df$sd2*1.96/sqrt(df$n)
checkb2=check1*check2

check1=b0[4]<df$b3+df$sd3*1.96/sqrt(df$n)
check2=b0[4]>df$b3-df$sd3*1.96/sqrt(df$n)
checkb3=check1*check2

df$error=error
df$Ivmax=(df$Iv==1)
df$checkb1=checkb1
df$checkb2=checkb2
df$checkb3=checkb3
#df=df[which(df$lambda %in% c(-4,-2,0,2,4)),]
best=df[which(df$Iv==1),]
best$lambda=100
df=rbind(df,best)
re=c()
for (n in n.sim){
  subdf=df[which(df$n==n & df$omega!=0 & df$sigma==5),]
  #avg=aggregate(cbind(subdf$error,subdf$Ivmax,subdf$Iv), list(subdf$lambda),FUN=mean)
  avg=aggregate(cbind(error,checkb1, checkb2, checkb3, Iv)~n+omega+sigma,data=subdf, FUN="mean")
  print(paste("n=",n, " count=",dim(subdf)[1],sep=""))
  print(avg)
  re=rbind(re,cbind(n,avg))
}


re=re[order(re$n, re$sigma, -re$omega),]
re=re[,c(2,3,4,5,6,7,8)]