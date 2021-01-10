
library(Rfast)
library(glmmLasso)

setwd('.../crejm-code')

# Load full training data-------------------------------
load(".../crejm-code/data/out.RData")

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

#---- Prepare data for GLMM Lasso ------------------------
set.seed(1)
player_select<-sample(1:n,1000,replace = F)

xx<-cbind(x,G)

#1. Login Indicator
yy<-matrix(0,n*(m-1),1)
yy[y>0]<-1

m1.data<- as.data.frame(cbind(yy,xx))
m1.data$grp<-as.factor(rep(1:n,each=m-1))
idx<-which(m1.data$grp %in% player_select,arr.ind = T)
m1.data<-m1.data[idx,]
vars <- names(m1.data[c(-1,-2,-28)])
fla <- as.formula(paste("V1~", paste(vars, collapse="+")))
lambda.grid<- c(0.1,1,10,20,50,100,500,1000)
gc()
m1.bic<-matrix(0,length(lambda.grid),1)
for (lam in 1:length(lambda.grid)){
  l<-lambda.grid[lam]
  m1.bic[lam]<- glmmLasso(fla,rnd=NULL,
                          data=m1.data, lambda=l, 
                          family = binomial(link="logit"),
                          switch.NR=FALSE, final.re=FALSE, control = list())$bic
  gc()
  print(lam)
}
l<- lambda.grid[which(m1.bic==min(m1.bic))]
m1<- glmmLasso(fla,rnd=NULL, data=m1.data, lambda=l, 
               family = binomial(link="logit"),
               switch.NR=FALSE, final.re=TRUE, control = list())
gc()

#2. Activity
yy.selected<-y[idx]
xx.selected<-xx[idx,]
idy<-(yy.selected>0)
yy.selected<-log(yy.selected[idy])
xx.selected<-xx.selected[idy,]

m2.data<- as.data.frame(cbind(yy.selected,xx.selected))
vars <- names(m2.data[c(-1,-2)])
fla <- as.formula(paste("yy.selected~", paste(vars, collapse="+")))
lambda.grid<- c(0.1,1,10,20,50,100,500,1000)
gc()
m2.bic<-matrix(0,length(lambda.grid),1)
for (lam in 1:length(lambda.grid)){
  l<-lambda.grid[lam]
  m2.bic[lam]<- glmmLasso(fla,rnd=NULL,
                          data=m2.data, lambda=l, 
                          family = gaussian(link="identity"),
                          switch.NR=FALSE, final.re=FALSE, control = list())$bic
  gc()
  print(lam)
}
l<- lambda.grid[which(m2.bic==min(m2.bic))]
m2<- glmmLasso(fla,rnd=NULL, data=m2.data, lambda=l, 
               family = gaussian(link="identity"),
               switch.NR=FALSE, final.re=TRUE, control = list())
gc()

#3. Purchase Indicator
yy.selected<-y[idx]
xx.selected<-xx[idx,]
evec.selected<-evec[idx]
idy<-(yy.selected>0)
evec.selected<-evec.selected[idy]
xx.selected<-xx.selected[idy,]

m3.data<- as.data.frame(cbind(evec.selected,xx.selected))
vars <- names(m3.data[c(-1,-2)])
fla <- as.formula(paste("evec.selected~", paste(vars, collapse="+")))
lambda.grid<- c(0.1,1,10,20,50,100,500,1000)
gc()
m3.bic<-matrix(0,length(lambda.grid),1)
for (lam in 1:length(lambda.grid)){
  l<-lambda.grid[lam]
  m3.bic[lam]<- glmmLasso(fla,rnd=NULL,
                          data=m3.data, lambda=l, 
                          family = binomial(link="logit"),
                          switch.NR=FALSE, final.re=FALSE, control = list())$bic
  gc()
  print(lam)
}
l<- lambda.grid[which(m3.bic==min(m3.bic))]
m3<- glmmLasso(fla,rnd=NULL, data=m3.data, lambda=l, 
               family = binomial(link="logit"),
               switch.NR=FALSE, final.re=TRUE, control = list())
gc()



#----------------- Now Prediction -------------------------
y.pred<-out$y.pred
evec.pred<-out$e.pred
x.pred<-cbind(rep(1,dim(out$x.pred)[1]),out$x.pred)
z.pred<-cbind(rep(1,dim(out$z.pred)[1]),out$z.pred)
gg.pred<-out$g.pred
guild_mem.pred<-out$gm.pred

m<-16
n<-dim(x.pred)[1]/(m-1)
pc<-dim(z.pred)[2]
pfpc<-dim(x.pred)[2]
pg<-dim(gg.pred)[2]
K<-dim(guild_mem.pred)[2]

# Create G matrix first n(m-1) by pg
G.pred<- matrix(0,n*(m-1),pg)
for(i in 1:n){
  
  s<-1+(i-1)*(m-1)
  e<-(m-1)*i
  temp1<-matrix(rep(c(guild_mem.pred[s:e,]),pg),ncol=pg,byrow = F)
  temp2<-temp1*gg.pred
  idx<-which(temp2!=0,arr.ind = T)#(temp2!=0)
  G.pred[s:e,]<-matrix(gg.pred[unique(idx[,1]),],ncol=pg,byrow = F)
}

xx.pred<-cbind(x.pred,G.pred)

#1. Login Indicator
yy.pred<-matrix(0,n*(m-1),1)
yy.pred[y.pred>0]<-1
m1.coeff<-as.vector(m1$coefficients)
temp<- as.matrix(xx.pred)%*%m1.coeff
yy_predicted<-1/(1+exp(-temp))
yy_predicted[yy_predicted>0.5]<- 1
yy_predicted[yy_predicted<=0.5]<- 0
yy_predicted.mat<-matrix(yy_predicted,ncol=m-1,byrow = T)
yy.pred.mat<-matrix(yy.pred,ncol=m-1,byrow = T)

C1 = yy_predicted.mat-yy.pred.mat
FP_login_pred = 100*mean(colSums(C1>0)/n)
FN_login_pred = 100*mean(colSums(C1<0)/n)

#2. Activity
activity.predicted<-matrix(0,n*(m-1),1)
yy.activity<-log(y.pred)
yy.activity[yy.activity==-Inf]<-0
m2.coeff<-as.vector(m2$coefficients)
activity.predicted<- exp(as.matrix(xx.pred)%*%m2.coeff+0.5*m2$phi)
activity.predicted<-log(activity.predicted)
activity.predicted[yy_predicted==0]<-0
C2 = abs(matrix(yy.activity,ncol=m-1,byrow = T)-
           matrix(activity.predicted,ncol=m-1,byrow = T))
err_activity_pred<-colSums(C2)/n


#3. Purchase Indicator
purchase_predicted<-matrix(0,n*(m-1),1)
m3.coeff<-as.vector(m3$coefficients)
temp<- as.matrix(xx.pred)%*%m3.coeff
purchase_predicted<-1/(1+exp(-temp))
purchase_predicted[purchase_predicted>0.5]<- 1
purchase_predicted[purchase_predicted<=0.5]<- 0
purchase_predicted[yy_predicted==0]<-0

C3= matrix(purchase_predicted,ncol = m-1,byrow = T)-
  matrix(evec.pred,ncol=m-1,byrow = T)
FP_purchase_pred = 100*mean(colSums(C3>0)/n)
FN_purchase_pred = 100*mean(colSums(C3<0)/n)

rm(list=c('guild_mem.pred','G.pred'))

save.image("glmmlasso_prediction.RData")
