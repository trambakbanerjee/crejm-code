
library(Rfast)
library(rpql)

setwd('.../crejm-code')

# Load full training data-------------------------------
load(".../crejm-code/out.RData")

y<-out$y.train
evec<-out$e.train
x<-cbind(rep(1,dim(out$x.train)[1]),out$x.train)
z<-cbind(rep(1,dim(out$z.train)[1]),out$z.train)
gg<-out$g.train
guild_mem<-out$gm.train

m<-16
n<-dim(x)[1]/(m-1)
pc<-dim(z)[2]
pfpc<-dim(x)[2]
pg<-dim(gg)[2]
K<-dim(guild_mem)[2]

# Create G matrix first n(m-1) by pg
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
grp<-rep(1:n,each=m-1)

#---- Prepare data for rpql ------------------------
set.seed(10)
player_select<-sample(1:n,250,replace = F)

#1. Login Indicator
yy<-matrix(0,n*(m-1),1)
yy[y>0]<-1

m1.data<- as.data.frame(cbind(yy,xx))
m1.data$grp<-grp
idx<-which(m1.data$grp %in% player_select,arr.ind = T)
m1.data<-m1.data[idx,]
zz<-list(cluster=m1.data[,c(2:15,17:19,22)])
idd<-list(cluster=as.numeric(m1.data$grp))
lambda.grid<- c(0.05,0.1,1,10,100)
m1<-rpqlseq(y=m1.data$V1,X=m1.data[,2:27],Z=zz,id=idd,
            family=binomial(), 1,lambda.grid,
            hybrid.est=TRUE)

#2. Activity
y.selected<-y[idx]
m2.data<-m1.data
m2.data$V1<-log(y.selected)
m2.data$V1[m2.data$V1==-Inf]<-0.001
zz<-list(cluster=m2.data[,c(2:15,17:19,22)])
idd<-list(cluster=as.numeric(m2.data$grp))
lambda.grid<- c(0.0001,0.001,0.01,0.1,1,10,100)
m2<-rpqlseq(y=m2.data$V1,X=m2.data[,2:27],Z=zz,id=idd,
            family=gaussian(),lambda=lambda.grid,
            hybrid.est=TRUE)


#3. Purchase Indicator
evec.selected<-evec[idx]
m3.data<-m1.data
m3.data$V1<-evec.selected
zz<-list(cluster=m3.data[,c(2:15,17:19,22)])
idd<-list(cluster=as.numeric(m3.data$grp))
lambda.grid<- c(0.0001,0.001,0.01,0.1,1,10,100)
m3<-rpqlseq(y=m3.data$V1,X=m3.data[,2:27],Z=zz,id=idd,
            family=binomial(), 1,lambda.grid,
            hybrid.est=TRUE)


#---------------------- Prediction -------------------
source('..../crejm-code/library/rcodelib_stableplayer_postselection.R')
D<-100
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

# Create G matrix first n(m-1) by pg
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
grp<-rep(1:n,each=m-1)
rm(list=c('out'))

m1fitout<-m1$best.fits$`BIC: log(nrow(X))*sum(beta!=0) + log(nrow(X))*sum(ranef_covmat!=0)`
m2fitout<-m2$best.fits$`BIC: log(nrow(X))*sum(beta!=0) + log(nrow(X))*sum(ranef_covmat!=0)`
m3fitout<-m3$best.fits$`BIC: log(nrow(X))*sum(beta!=0) + log(nrow(X))*sum(ranef_covmat!=0)`

q1<-m1fitout$nonzero.ranef$cluster # no random effects selected for model 1
q2<-m2fitout$nonzero.ranef$cluster #
q3<-m3fitout$nonzero.ranef$cluster # 
Sigma2=m2fitout$ran.cov$cluster[q2,q2]
Sigma3=m3fitout$ran.cov$cluster[q3,q3]

pii = matrix(0,n,m-2)
pii_polished = pii
qii = matrix(0,n,m-2)
qii_polished = qii
A = matrix(0,n,m-2)

Xb1 = as.matrix(xx)%*%as.vector(m1fitout$fixef)
Xb2 = as.matrix(xx)%*%as.vector(m2fitout$fixef)
Xb3 = as.matrix(xx)%*%as.vector(m3fitout$fixef)

for(i in 1:n){
  
  s = (i-1)*(m-1)+2;
  e = i*(m-1);
  
  B<-matrix(0,D,3*pc)
  set.seed(i)
  b3<-rnorm(D,0,sqrt(Sigma3))
  set.seed(10*i)
  b2<-rmvnorm(D,rep(0,length(q2)),Sigma2)
  B[,32+q3]<-b3#scale(b3,center = TRUE,scale=TRUE)#
  B[,16+q2]<-b2#scale(b2,center = TRUE,scale=TRUE)#b2
  
  

  bhat<-empbayes_rpql_0720(B,y,evec,xx,z,m1fitout$fixef,m2fitout$fixef,
                     m3fitout$fixef,solve(Sigma2),solve(Sigma3),
                     m2fitout$phi,s,e,q1,q2,q3) 
  
  #Login Indicator
  temp = Xb1[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,1:18])
  pii[i,] = (1+exp(-temp))^(-1);
  pii_polished[i,] = 1*(pii[i,]>0.5);
  
  #Activity Time
  temp = Xb2[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,19:36])
  A[i,] = exp(temp+0.5*ceiling(m2fitout$phi))*(t(pii_polished[i,]))#temp*(t(pii_polished[i,]))#
  
  #Purchase Indicator
  temp = Xb3[s:e]+rowSums(as.matrix(z[s:e,])*bhat[,37:54])
  qii[i,] = ((1+exp(-temp))^(-1))*(t(pii_polished[i,]))
  qii_polished[i,] = 1*(qii[i,]>0.5)
 
  print(i)
  
}
#----- Compute metrices---------------------------------
yy<-y
yy[y>0]<-1
deltamat<-matrix(yy,ncol=m-1,byrow = T)
deltamat<-deltamat[,2:(m-1)]
C1 = pii_polished-deltamat
FP_login = 100*mean(colSums(C1>0)/n)
FN_login = 100*mean(colSums(C1<0)/n)

Alog <-log(A)
Alog[Alog==-Inf]<- 0
ymat <- matrix(y,ncol=m-1,byrow = T)
ymat = log(ymat)
ymat[ymat==-Inf]<- 0
C2 = abs(ymat[,2:(m-1)]-Alog)
err_activity<- as.matrix(colSums(C2)/n)

alphamat<-matrix(evec,ncol=m-1,byrow = T)
alphamat<-alphamat[,2:(m-1)]
C3 = qii_polished-alphamat
FP_purchase = 100*mean(colSums(C3>0)/n)
FN_purchase = 100*mean(colSums(C3<0)/n)

save.image("rpql_prediction.RData")

