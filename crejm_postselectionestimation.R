
library(Rfast)
library(igraph)
library(ggb)



#----------- Load Selection Results --------------
load(".../crejm-code/selection.RData")
rm(list=setdiff(ls(), c("y","evec","x","z","gg","guild_mem",
                        "output_ccrem","bic.mat")))

setwd('.../crejm-code/spcov')
files.sources = list.files()
sapply(files.sources, source)
setwd('.../crejm-code/library')
source('rcodelib_stableplayer_postselection.R')
sourceCpp('codelib_postselection.cpp')
setwd('.../crejm-code')
#------------------------------------------------------------------------------------

model_postselection<-cleanup_selection(output_ccrem[[4]],y,x,z,gg)
Gamma1.0<-model_postselection$Gamma.1
Gamma2.0<-model_postselection$Gamma.2
Gamma3.0<-model_postselection$Gamma.3
Sigma1.0<-model_postselection$Sigma.1
Sigma2.0<-model_postselection$Sigma.2

x1<-model_postselection$x1
x2<-model_postselection$x2
x3<-model_postselection$x3
z1<-model_postselection$z1
z2<-model_postselection$z2
z3<-model_postselection$z3
gg1<-as.matrix(model_postselection$g1)
gg2<-as.matrix(model_postselection$g2)
gg3<-as.matrix(model_postselection$g3)

m<-16
n<-dim(x1)[1]/(m-1)
K<-dim(guild_mem)[2]
t.cut<-3
D<-500
itermax<-100


postselection.est<-dccjm_postselection_est(y,evec,x1,x2,x3,z1,z2,z3,gg1,gg2,gg3,
                    Gamma1.0,Gamma2.0,Gamma3.0,
                    Sigma1.0,Sigma2.0,guild_mem,
                    t.cut,n,D,itermax)

save(postselection.est,file='postselection.est.RData')

