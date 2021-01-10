
library(Rfast)
library(igraph)
library(ggb)
library(GLMMadaptive)
library(lme4)

setwd('.../crejm-code/spcov')
files.sources = list.files()
sapply(files.sources, source)
setwd('.../crejm-code/library')
source('rcodelib_stableplayer.R')
sourceCpp('codelib.cpp')
setwd('.../crejm-code')
#------------------ load real data ----------------------------
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
  idx<-which(temp2!=0,arr.ind = T)
  G[s:e,]<-matrix(gg[unique(idx[,1]),],ncol=pg,byrow = F)
}

#---- Prepare data for GLMM Adaptive ------------------------
xx<-cbind(x,G)
DF <- data.frame(xx)
DF$y<-matrix(0,n*(m-1),1)
DF$y[y>0]<-1
id <- rep(1:n, each = m-1)
DF$id <- factor(id)

m1<-mixed_model(fixed = y ~ X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26,
                random = ~  1 | id, data = DF,
                family = binomial())

DF <- data.frame(xx)
DF$y<-matrix(0,n*(m-1),1)
DF$y[y>0]<-log(y[y>0])
id <- rep(1:n, each = m-1)
DF$id <- factor(id)

m2<-lmer(y ~ X2+X3+X4++X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+(1 | id),
         data = DF)

coeff_m2<-summary(m2)


DF <- data.frame(xx)
DF$y<-matrix(0,n*(m-1),1)
DF$y[evec>0]<-1
id <- rep(1:n, each = m-1)
DF$id <- factor(id)

m3<-mixed_model(fixed = y ~ X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26,
                random = ~  1 | id, data = DF,
                family = binomial())


#----------- Initial Estimates --------------
t.cut<-3
D<-500
itermax<-50
Gamma1.0<-matrix(m1$coefficients,pfpc+pg,1)
Gamma2.0<-matrix(c(coeff_m2$coefficients[,1]/coeff_m2$sigma,1/coeff_m2$sigma),
                 pfpc+pg+1,1)
Gamma3.0<-matrix(m3$coefficients,pfpc+pg,1)
Sigma1.0<-posdef(3*pc)
Sigma2.0<-gen_bandedcov(m-1,t.cut,0.5)


init.est<-dccjm_est(y,evec,x,z,gg,Gamma1.0,Gamma2.0,Gamma3.0,
                    Sigma1.0,Sigma2.0,guild_mem,
                    t.cut,n,D,itermax)

save(init.est,file='init.est.RData')
