
library(Rfast)
library(igraph)
library(ggb)

setwd('.../crejm-code/spcov')
files.sources = list.files()
sapply(files.sources, source)
setwd('.../crejm-code/library')
source('rcodelib_stableplayer.R')
sourceCpp('codelib.cpp')
setwd('.../crejm-code')
# ---- Load full data-------------------------------
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
t.cut<-3
D<-500
itermax<-5

# Load Initial Estimates -------------------------
load(".../crejm-code/init.est.RData")

Sigma1.0<-init.est$Sigma1.est
Sigma2.0<-init.est$Sigma2.est
Gamma1.0<-init.est$Gamma1.est
Gamma1.0[1]<-sign(Gamma1.0[1])
Gamma2.0<-init.est$Gamma2.est
Gamma3.0<-init.est$Gamma3.est
Gamma3.0[1]<-sign(Gamma3.0[1])
A.1<-0*diag(length(Gamma1.0))
a.1<-A.1%*%Gamma1.0
A.2<-matrix(c(rep(0,length(Gamma2.0)-1),-1),1,length(Gamma2.0))
a.2<- -0.1
A.3<-0*diag(length(Gamma3.0))
a.3<-A.3%*%Gamma3.0

lambda_mat<-matrix(c(0.001,0.001,0.5,
                     0.005,0.005,0.5,
                     0.01,0.01,0.5,
                     0.05,0.05,0.5,
                     0.1,0.1,0.5),ncol=3,byrow=T)
e_ind<-c(1,2,3,4)

bic.mat<-matrix(0,dim(lambda_mat)[1],1)
output_ccrem<-vector("list",dim(lambda_mat)[1])

for(i1 in 1:dim(lambda_mat)[1]){
  lam1<-lambda_mat[i1,1]
  lam2<-lambda_mat[i1,2]
  lam3<-lambda_mat[i1,3]
  
  out_ccrem<-dccjm_select(y,evec,x,z,gg,Gamma1.0,Gamma2.0,Gamma3.0,
                          Sigma1.0,Sigma2.0,guild_mem,
                          A.1,a.1,A.2,a.2,A.3,a.3,
                          n,t.cut,lam1,lam2,lam3,D,itermax)
  
  # prepare for BIC calculation
  Gamma1.new<-out_ccrem$Gamma1.est
  Gamma2.new<-out_ccrem$Gamma2.est
  Gamma3.new<-out_ccrem$Gamma3.est
  Sigma1.new<-out_ccrem$Sigma1.est
  Sigma2.new<-out_ccrem$Sigma2.est
  k<-sum(1*(Gamma1.new!=0))+sum(1*(Gamma2.new!=0))+sum(1*(Gamma3.new!=0))+
    sum(1*(abs(Sigma1.new)>0.025))+
    sum(1*(abs(Sigma2.new)>0.01))
  
  bic.mat[i1]<- -2*out_ccrem$lQ.val+log(n)*k
  
  #store
  output_ccrem[[i1]]<-out_ccrem
  
  # prepare for warm start for the next lambda sequence
  if(i1<(dim(lambda_mat)[1])){
    Gamma1.0<-output_ccrem[[e_ind[i1]]]$Gamma1.est
    Gamma2.0<-output_ccrem[[e_ind[i1]]]$Gamma2.est
    Gamma3.0<-output_ccrem[[e_ind[i1]]]$Gamma3.est
    Sigma1.0<-output_ccrem[[e_ind[i1]]]$Sigma1.est
    Sigma2.0<-output_ccrem[[e_ind[i1]]]$Sigma2.est
  }
}

save.image("selection.RData")
