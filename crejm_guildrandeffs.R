
library(Rfast)
library(igraph)
library(ggb)
library(foreach)
library(doParallel)

setwd('.../crejm-code')
load(".../crejm-code/selection.RData")
rm(list=setdiff(ls(), c("output_ccrem")))
source('.../crejm-code/library/rcodelib_stableplayer_postselection.R')
sourceCpp('.../crejm-code/library/codelib_postselection.cpp')

idx.1<-which(output_ccrem[[4]]$Gamma1[1:26]!=0)
idx.2<-which(output_ccrem[[4]]$Gamma2[1:27]!=0)
idx.3<-which(output_ccrem[[4]]$Gamma3[1:26]!=0)

idz.1<-which(diag(output_ccrem[[4]]$Sigma1[1:18,1:18])>0.015)
idz.2<-which(diag(output_ccrem[[4]]$Sigma1[19:36,19:36])>0.015)
idz.3<-which(diag(output_ccrem[[4]]$Sigma1[37:54,37:54])>0.015)

load(".../crejm-code/postselection.est.RData")

#1. Load data---------------------------------------
load(".../crejm-code/out.RData")
y<-out$y.pred
evec<-out$e.pred
x<-cbind(rep(1,dim(out$x.pred)[1]),out$x.pred)
z<-cbind(rep(1,dim(out$z.pred)[1]),out$z.pred)
gg<-out$g.pred
guild_mem<-out$gm.pred

m<-16
n<-dim(x)[1]/(m-1)
pc<-dim(z)[2]
pfpc<-dim(x)[2]
pg<-dim(gg)[2]
K<-dim(guild_mem)[2]
t.cut<-3
D<-1000
rm('out')

#2. Create G matrix first n(m-1) by pg
G<- matrix(0,n*(m-1),pg)
for(i in 1:n){
  
  s<-1+(i-1)*(m-1)
  e<-(m-1)*i
  temp1<-matrix(rep(c(guild_mem[s:e,]),pg),ncol=pg,byrow = F)
  temp2<-temp1*gg
  idx<-which(temp2!=0,arr.ind = T)#(temp2!=0)
  G[s:e,]<-matrix(gg[unique(idx[,1]),],ncol=pg,byrow = F)
}
xx<-cbind(x,G)

#3. Create vector of parameter estimates
Gamma1.final<-matrix(0,dim(xx)[2],1)
Gamma1.final[idx.1]<-postselection.est$Gamma1.est
Gamma2.final<-matrix(0,dim(xx)[2]+1,1)
Gamma2.final[idx.2]<-postselection.est$Gamma2.est
Gamma3.final<-matrix(0,dim(xx)[2],1)
Gamma3.final[idx.3]<-postselection.est$Gamma3.est
Sigma1.final<-postselection.est$Sigma1.est
Sigma2.final<-postselection.est$Sigma2.est

#4. Begin prediction

Xb1 = as.matrix(xx)%*%Gamma1.final
Xb2 = as.matrix(xx)%*%(Gamma2.final[1:26]/Gamma2.final[27])
Xb3 = as.matrix(xx)%*%Gamma3.final

set.seed(42000) 
C = sample_re(K*D,Sigma2.final)
C0 = process_guildre(C,guild_mem,n,K,D)
ii<-1:(n*(m-1))
C1 = C0[ii,]
C2 = C0[ii+(n*(m-1)),]
C3 = C0[ii+(2*n*(m-1)),]

Sigma1.inv<-solve(Sigma1.final)
Sigma2.inv<-solve(Sigma2.final)

cl <- makeCluster(10)
registerDoParallel(cl)

result<-foreach(i = 1:n,.packages=c("Rfast"))%dopar%{

  s = (i-1)*(m-1)+2;
  e = i*(m-1);

  B<-matrix(0,D,3*pc)
  set.seed(i)
  b<-rmvnorm(D,rep(0,dim(Sigma1.final)[1]),Sigma1.final)
  B[,c(idz.1,18+idz.2,36+idz.3)]<-b

  bchat<-empbayes_ccrejm(B,C,C1,C2,C3,y,evec,xx,z,Gamma1.final,Gamma2.final[1:26],
                      Gamma3.final,Sigma1.inv,Sigma2.inv,
                      1/Gamma2.final[27],s,e,idz.1,idz.2,idz.3)
  bhat<-bchat$B
  C1hat<-bchat$C1
  C2hat<-bchat$C2
  C3hat<-bchat$C3
  idx<- which(guild_mem[s:e,]>0,arr.ind = T)[,2]

  return(list("C1hat"=C1hat,"C2hat"=C2hat,"C3hat"=C3hat,"mem" = idx))

}

stopCluster(cl)
registerDoSEQ()

rm(list=c('x','z','C0','C1','C2','C3','C','G','xx','y','Xb1','Xb2','Xb3'))

#----- Get plots ---------------------------------
library(ggplot2)
library(ggpubr)

idx<-c(4,13,42)
nn<-n
C1hat<- sapply(1:nn,function(i) result[[i]]$C1hat)
C2hat<- sapply(1:nn,function(i) result[[i]]$C2hat)
C3hat<- sapply(1:nn,function(i) result[[i]]$C3hat)
memb<- c(sapply(1:nn,function(i) result[[i]]$mem))

C1_cum<- c(apply(C1hat, 2,cumsum))
C2_cum<- c(apply(C2hat, 2,cumsum))
C3_cum<- c(apply(C3hat, 2,cumsum))

C1<- c(C1hat)
C2<- c(C2hat)
C3<- c(C3hat)

plotdata1<- as.data.frame(cbind(C1_cum,C1,memb,rep(1:14,nn)))
ii<-(plotdata1$memb==27)
plotdata1<- plotdata1[!ii,]
g1<-ggplot()+
  geom_smooth(data=plotdata1,aes(x=V4,y=C1,group=memb),color="gray",
              size=0.2,alpha = 0.1,span=0.25,se=FALSE,method = "loess")+
  geom_smooth(data=plotdata1[plotdata1$memb%in%idx,],aes(x=V4,y=C1,color=factor(memb)),
              size=2,span=0.25,se=FALSE,method = "loess")+
  theme_bw()+labs(color = "Guilds")+
  theme(legend.position = "top",
        legend.text=element_text(size=15),legend.title = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+coord_cartesian(ylim = c(-0.35,0.2))+
  ylab("Guild Random Effects - Login Ind.")+
  xlab("Days")+scale_x_continuous(breaks=1:14,labels = 17:30)+
  scale_y_continuous(breaks = c(-0.35,-0.25,-0.15,-0.05,0.05,0.15))

plotdata2<- as.data.frame(cbind(C2_cum,C2,memb,rep(1:14,nn)))
ii<-(plotdata2$memb==27)
plotdata2<- plotdata2[!ii,]
g2<-ggplot()+
  geom_smooth(data=plotdata2,aes(x=V4,y=C2,group=memb),color="gray",
              size=0.2,alpha = 0.1,span=0.25,se=FALSE,method = "loess")+
  geom_smooth(data=plotdata2[plotdata2$memb%in%idx,],aes(x=V4,y=C2,color=factor(memb)),
              size=2,alpha = 0.2,span=0.25,se=FALSE,method = "loess")+
  theme_bw()+labs(color = "Guilds")+
  theme(legend.position = "top",
        legend.text=element_text(size=15),legend.title = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+coord_cartesian(ylim = c(-0.6,0.1))+
  ylab("Guild Random Effects - Duration")+
  xlab("Days")+scale_x_continuous(breaks=1:14,labels = 17:30)

plotdata3<- as.data.frame(cbind(C3_cum,C3,memb,rep(1:14,nn)))
ii<-(plotdata3$memb==27)
plotdata3<- plotdata3[!ii,]
g3<-ggplot()+
  geom_smooth(data=plotdata3,aes(x=V4,y=C3,group=memb),color="gray",
              size=0.2,alpha = 0.1,span=0.25,se=FALSE,method = "loess")+
  geom_smooth(data=plotdata3[plotdata3$memb%in%idx,],aes(x=V4,y=C3,color=factor(memb)),
              size=2,alpha = 0.2,span=0.25,se=FALSE,method = "loess")+
  theme_bw()+labs(color = "Guilds")+
  theme(legend.position = "top",legend.title = element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+coord_cartesian(ylim = c(-0.3,0.2))+
  ylab("Guild Random Effects - Purch Prop. ")+
  xlab("Days")+scale_x_continuous(breaks=1:14,labels = 17:30)

ggarrange(g1,g2,g3,ncol=3,common.legend = TRUE,legend="top")

save.image("crejm_guildrandeffs.RData")


