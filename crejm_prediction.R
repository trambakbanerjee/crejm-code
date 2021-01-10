
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
pii = matrix(0,n,m-2)
pii_polished = pii
qii = matrix(0,n,m-2)
qii_polished = qii
A = matrix(0,n,m-2)

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

  #Login Indicator
  temp = Xb1[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,1:18])+C1hat
  pii = (1+exp(-temp))^(-1);
  pii_polished = 1*(pii>0.5);

  #Activity Time
  temp = Xb2[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,19:36])+C2hat
  A = exp(temp+0.5*(1/Gamma2.final[27]))*(pii_polished)

  #Purchase Indicator
  temp = Xb3[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,37:54])+C3hat
  qii = ((1+exp(-temp))^(-1))*(pii_polished)
  qii_polished = 1*(qii>0.5)

  return(list("pii_polished"=pii_polished,
              "A"=A,
              "qii_polished"=qii_polished))

}

stopCluster(cl)
registerDoSEQ()

pii_polished<-t(sapply(1:n,function(i) result[[i]]$pii_polished))
A<-t(sapply(1:n,function(i) result[[i]]$A))
qii_polished<-t(sapply(1:n,function(i) result[[i]]$qii_polished))

#----- Compute metrices---------------------------------
yy<-y
yy[y>0]<-1
deltamat<-matrix(yy,ncol=m-1,byrow = T)
deltamat<-deltamat[1:n,2:(m-1)]
C11 = pii_polished-deltamat
FP_login = 100*mean(colSums(C11>0)/n)
FN_login = 100*mean(colSums(C11<0)/n)

Alog <-log(A)
Alog[Alog==-Inf]<- 0
ymat <- matrix(y,ncol=m-1,byrow = T)
ymat = log(ymat)
ymat[ymat==-Inf]<- 0
C22 = abs(ymat[,2:(m-1)]-Alog)
err_activity<- as.matrix(colSums(C22)/n)

alphamat<-matrix(evec,ncol=m-1,byrow = T)
alphamat<-alphamat[,2:(m-1)]
C33 = qii_polished-alphamat
FP_purchase = 100*mean(colSums(C33>0)/n)
FN_purchase = 100*mean(colSums(C33<0)/n)


save.image("ccrejm_prediction.RData")

