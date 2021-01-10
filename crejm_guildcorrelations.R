
library(Rfast)
library(foreach)
library(doParallel)
library(ggplot2)
library(gridExtra)

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
D<-500
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

Xb1 = as.matrix(xx)%*%Gamma1.final
Xb2 = as.matrix(xx)%*%(Gamma2.final[1:26]/Gamma2.final[27])
Xb3 = as.matrix(xx)%*%Gamma3.final

#-------------------------------------------------------------------------
#---- First get the predicted means and second moments for all players ---
#-------------------------------------------------------------------------

set.seed(42000)
C = sample_re(K*D,Sigma2.final)
C0 = process_guildre(C,guild_mem,n,K,D)
ii<-1:(n*(m-1))
C1 = C0[ii,]
C2 = C0[ii+(n*(m-1)),]
C3 = C0[ii+(2*n*(m-1)),]

rm('C0')

Sigma1.inv<-solve(Sigma1.final)
Sigma2.inv<-solve(Sigma2.final)

cl <- makeCluster(20)
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

  #Login Indicator
  temp = Xb1[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,1:18])+C1hat
  pii = (1+exp(-temp))^(-1);
  pii_var = pii*(1-pii)

  #Activity Time
  temp = Xb2[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,19:36])+C2hat
  A = exp(temp+0.5*(1/Gamma2.final[27]))*(1*(pii>0.5))
  A_var = (1*(pii>0.5))*exp(2*temp+2*(1/Gamma2.final[27]))-(A^2)

  #Purchase Indicator
  temp = Xb3[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,37:54])+C3hat
  qii = ((1+exp(-temp))^(-1))*(1*(pii>0.5))
  qii_var = qii*(1-qii)


  return(list("pii"=pii,"pii_var" = pii_var,
              "A"=A,"A_var"=A_var,
              "qii"=qii,"qii_var"=qii_var))

}

stopCluster(cl)
registerDoSEQ()

pii_mean<-t(sapply(1:n,function(i) result[[i]]$pii))
pii_var<-t(sapply(1:n,function(i) result[[i]]$pii_var))
A_mean<-t(sapply(1:n,function(i) result[[i]]$A))
A_var<-t(sapply(1:n,function(i) result[[i]]$A_var))
qii_mean<-t(sapply(1:n,function(i) result[[i]]$qii))
qii_var<-t(sapply(1:n,function(i) result[[i]]$qii_var))

save.image("crejm_guildcorrs_part1.RData")

#-------------------------------------------------------------------------
#--- Estimate player Covariances by each guild ---------------------------
#-------------------------------------------------------------------------

idd.time<-matrix(1:(n*(m-1)),ncol=m-1,byrow = T)
#-------------------------------------------------------------
#1. Start with login indicator
#----------------------------------------------------
guild_corr_login<- list()
guild_corr_act<- list()
guild_corr_purch<- list()

# Guilds 1:K
cl <- makeCluster(20)
registerDoParallel(cl)

result_guildcorr<-foreach(g=1:K,.packages = c("Rfast"))%dopar%{
  
  guild<- matrix(guild_mem[,g],ncol=m-1,byrow=T)
  
  # Times 2:15
  l1<-list()
  l2<-list()
  l3<-list()
  for(tt in 2:(m-1)){
    idd.player<- which(guild[,tt]>0)
    if(length(idd.player)>1){
      ij_login<- ij_act<- ij_purch<- matrix(0,length(idd.player),length(idd.player))
      for(i in 1:length(idd.player)){
        xx1<- xx[idd.time[idd.player[i],1:(tt-1)],]
        y1<- y[idd.time[idd.player[i],1:(tt-1)]]
        evec1<- evec[idd.time[idd.player[i],1:(tt-1)]]
        z1<- z[idd.time[idd.player[i],1:(tt-1)],]
        B1<-matrix(0,D,3*pc)
        set.seed(i)
        b<-rmvnorm(D,rep(0,dim(Sigma1.final)[1]),Sigma1.final)
        B1[,c(idz.1,18+idz.2,36+idz.3)]<-b
        if(i<length(idd.player)){
        for(v in (i+1):length(idd.player)){
          
            xx2<- xx[idd.time[idd.player[v],1:(tt-1)],]
            y2<- y[idd.time[idd.player[v],1:(tt-1)]]
            evec2<- evec[idd.time[idd.player[v],1:(tt-1)]]
            z2<- z[idd.time[idd.player[v],1:(tt-1)],]
            B2<-matrix(0,D,3*pc)
            set.seed(v)
            b<-rmvnorm(D,rep(0,dim(Sigma1.final)[1]),Sigma1.final)
            B2[,c(idz.1,18+idz.2,36+idz.3)]<-b
            
          bchat<-empbayes_ccrejm_guildcorr(B1,B2,C,
                                 matrix(C1[idd.time[idd.player[i],1:(tt-1)],],nrow=tt-1),
                                 matrix(C2[idd.time[idd.player[i],1:(tt-1)],],nrow=tt-1),
                                 matrix(C3[idd.time[idd.player[i],1:(tt-1)],],nrow=tt-1),
                                  y1,evec1,xx1,z1,y2,evec2,xx2,z2,
                                  Gamma1.final,Gamma2.final[1:26],
                                 Gamma3.final,Sigma1.inv,Sigma2.inv,
                                 1/Gamma2.final[27],idz.1,idz.2,idz.3,tt)
          bhat1<-bchat$B1
          bhat2<-bchat$B2
          C1hat<-bchat$C1
          C2hat<-bchat$C2
          C3hat<-bchat$C3
          
          idx1<- idd.time[idd.player[i],tt]
          idx2<- idd.time[idd.player[v],tt]
          
          #Login Indicator
          temp1 = Xb1[idx1]+sum(as.matrix(z[idx1,])*bhat1[1:18])+C1hat[tt-1]
          temp2 = Xb1[idx2]+sum(as.matrix(z[idx2,])*bhat2[1:18])+C1hat[tt-1]
          ij_login[i,v]<- ((1+exp(-temp1))^(-1))*((1+exp(-temp2))^(-1))
          
          #Activity Time
          temp1 = Xb2[idx1]+sum(as.matrix(z[idx1,])*bhat1[19:36])+C2hat[tt-1]
          temp2 = Xb2[idx2]+sum(as.matrix(z[idx2,])*bhat2[19:36])+C2hat[tt-1]
          ij_act[i,v]<- exp(temp1+0.5*(1/Gamma2.final[27]))*exp(temp2+0.5*(1/Gamma2.final[27]))
          
          #Purchase Indicator
          temp1 = Xb3[idx1]+sum(as.matrix(z[idx1,])*bhat1[37:54])+C3hat[tt-1]
          temp2 = Xb3[idx2]+sum(as.matrix(z[idx2,])*bhat2[37:54])+C3hat[tt-1]
          ij_purch[i,v]<- ((1+exp(-temp1))^(-1))*((1+exp(-temp2))^(-1))
          #print(c(i,v))
        }
        }
      }
    }
    cov_login<- ij_login-(pii_mean[idd.player,tt-1]%*%t(pii_mean[idd.player,tt-1]))
    corr_login<- cov_login/sqrt(pii_var[idd.player,tt-1]%*%t(pii_var[idd.player,tt-1]))
    corr_login<- corr_login[upper.tri(corr_login)]
    corr_login[corr_login< -1]<- -1
    corr_login[corr_login> 1]<- 1
    l1[[tt-1]]<- corr_login
    
    cov_act<- ij_act-(A_mean[idd.player,tt-1]%*%t(A_mean[idd.player,tt-1]))
    corr_act<- cov_act/sqrt(A_var[idd.player,tt-1]%*%t(A_var[idd.player,tt-1]))
    corr_act<- corr_act[upper.tri(corr_act)]
    corr_act<- corr_act[corr_act<Inf]
    corr_act[corr_act< -1]<- -1
    corr_act[corr_act> 1]<- 1
    l2[[tt-1]]<- corr_act
    
    cov_purch<- ij_purch-(qii_mean[idd.player,tt-1]%*%t(qii_mean[idd.player,tt-1]))
    corr_purch<- cov_purch/sqrt(qii_var[idd.player,tt-1]%*%t(qii_var[idd.player,tt-1]))
    corr_purch<- corr_purch[upper.tri(corr_purch)]
    corr_purch<- corr_purch[corr_purch<Inf]
    corr_purch[corr_purch< -1]<- -1
    corr_purch[corr_purch> 1]<- 1
    l3[[tt-1]]<- corr_purch
  }

  return(list("guild_corr_login"=l1,"guild_corr_act" = l2,
              "guild_corr_purch"=l3))
  
}

stopCluster(cl)
registerDoSEQ()

save.image("crejm_guildcorrs_part2.RData")

#---- get the plots -----
library(lattice)
library(viridisLite)
library(RColorBrewer)
library(ggpubr)

# Activity

plotdata_act<-matrix(NA,ncol=3)

for(k in 1:49){
  for(t in 1:14){
    
    tempdata<- result_guildcorr[[k]]$guild_corr_act[[t]]
    plotdata_act<-rbind(plotdata_act,c(mean(tempdata),t,k))
  }
}
plotdata_act<- plotdata_act[-1,]
plotmat2<-matrix(plotdata_act[,1],ncol=14,byrow = T)

coul <- colorRampPalette(brewer.pal(11, "Spectral"))
h2<-levelplot(t(plotmat2),col.regions=coul,at=seq(-0.1,0.3,length.out = 17),
              xlab="days",ylab="Guilds",main=list('Activity',side=1,line=0.5))
# Login

plotdata_login<-matrix(NA,ncol=3)

for(k in 1:49){
  for(t in 1:14){
    
     tempdata<- result_guildcorr[[k]]$guild_corr_login[[t]]
     plotdata_login<-rbind(plotdata_login,c(mean(tempdata),t,k))
  }
}
plotdata_login<- plotdata_login[-1,]
plotmat1<-matrix(plotdata_login[,1],ncol=14,byrow = T)


coul <- colorRampPalette(brewer.pal(11, "Spectral"))
h1<-levelplot(t(plotmat1),col.regions=coul, at=seq(-0.1,0.3,length.out = 17),
              xlab="days",ylab="Guilds",main=list('Login Ind.',side=1,line=0.5))


# Purchase Propensity

plotdata_purch<-matrix(NA,ncol=3)

for(k in 1:49){
  for(t in 1:14){
    
    tempdata<- result_guildcorr[[k]]$guild_corr_purch[[t]]
    plotdata_purch<-rbind(plotdata_purch,c(mean(tempdata),t,k))
  }
}
plotdata_purch<- plotdata_purch[-1,]
plotmat3<-matrix(plotdata_purch[,1],ncol=14,byrow = T)

coul <- colorRampPalette(brewer.pal(11, "Spectral"))
h3<-levelplot(t(plotmat3),col.regions=coul,at=seq(-0.1,0.3,length.out = 17),
              xlab="days",ylab="Guilds",main=list('Purchase Prop.',side=1,line=0.5))

ggarrange(h1,h2,h3,nrow=1,ncol=3)

#----------------- Functional Cluster Analysis -----------
#--------------------------------------------------------
library(fda.usc)

do1<- plotmat1
do1.smth<-matrix(0,K-1,m-2)

for(i in 1:(K-1)){
  
  do1.smth[i,]<- smooth.spline(1:(m-2),do1[i,],df=7)$y
  
  
}
do1<-do1.smth
n<- nrow(do1.smth)

set.seed(10)
out<- kmeans.fd(fdata(do1),ncl=2,method="exact",max.iter=50000)
plot(out$centers)
table(out$cluster)

centers<-as.data.frame(c(t(out$centers$data)))
names(centers)<-c('centers')
centers$clust<- factor(rep(1:2,each=14))
centers$t = rep(17:30,2)

g1<- ggplot() + geom_line(data=centers,aes(x=t,y=centers,color=clust),size=2)+
  xlab('Days')+ylab("Correlation (Login Ind.)")+scale_x_continuous(breaks = 17:30)+
  theme_bw()+labs(color = "Cluster Centroids")+
  theme(legend.position = "top",legend.title = element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plot(g1)
#----------------------------------------------------------
do2<- plotmat2
do2.smth<-matrix(0,K-1,m-2)

for(i in 1:(K-1)){
  
  do2.smth[i,]<- smooth.spline(1:(m-2),do2[i,],df=7)$y
  
  
}
do2<-do2.smth
n<- nrow(do2.smth)

set.seed(10)
out<- kmeans.fd(fdata(do2),ncl=3,method="exact",max.iter=50000)
plot(out$centers)
table(out$cluster)

centers<-as.data.frame(c(t(out$centers$data)))
names(centers)<-c('centers')
centers$clust<- factor(rep(1:3,each=14))
centers$t = rep(17:30,3)

g2<- ggplot() + geom_line(data=centers,aes(x=t,y=centers,color=clust),size=2)+
  xlab('Days')+ylab("Correlation (Duration)")+scale_x_continuous(breaks = 17:30)+
  theme_bw()+labs(color = "Cluster Centroids")+
  theme(legend.position = "top",legend.title = element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
plot(g2)
#-------------------------------------------------------
do3<- plotmat3
do3.smth<-matrix(0,K-1,m-2)

for(i in 1:(K-1)){
  
  do3.smth[i,]<- smooth.spline(1:(m-2),do3[i,],df=7)$y
  
  
}
do3<-do3.smth
n<- nrow(do3.smth)

set.seed(10)
out<- kmeans.fd(fdata(do3),ncl=3,method="exact",max.iter=50000)
plot(out$centers)
table(out$cluster)

centers<-as.data.frame(c(t(out$centers$data)))
names(centers)<-c('centers')
centers$clust<- factor(rep(1:3,each=14))
centers$t = rep(17:30,3)

g3<- ggplot() + geom_line(data=centers,aes(x=t,y=centers,color=clust),size=2)+
  xlab('Days')+ylab("Correlation (Purchase Prop.)")+scale_x_continuous(breaks = 17:30)+
  theme_bw()+labs(color = "Cluster Centroids")+theme(legend.position = "top",legend.title = element_text(size=15),
                                                     legend.text=element_text(size=15),
                                                     axis.text.x = element_text(size=15),
                                                     axis.text.y = element_text(size=15),
                                                     axis.title.x = element_text(size=15),
                                                     axis.title.y = element_text(size=15))
plot(g3)
#--------------------------------------------------------
library(ggpubr)
ggarrange(g1,g2,g3,ncol=3)