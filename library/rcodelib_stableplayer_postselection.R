
library(Rfast)
library(msos)
library(mvtnorm)
library(CVXR)


#------------------ Modelling scripts -----------------------------
cleanup_selection<-function(selection_ccrem,y,x,z,gg){
  
  Gamma1<-selection_ccrem$Gamma1.est
  Gamma2<-selection_ccrem$Gamma2.est
  Gamma3<-selection_ccrem$Gamma3.est
  Sigma1<-selection_ccrem$Sigma1.est
  Sigma2<-selection_ccrem$Sigma2.est
  
  idx<-which(Gamma1[1:21]!=0)
  x1<-x[,idx]
  idx<-which(Gamma2[1:21]!=0)
  x2<-x[,idx]
  idx<-which(Gamma3[1:21]!=0)
  x3<-x[,idx]
  
  idx<-which(Gamma1[22:26]!=0)
  g1<-gg[,idx]
  idx<-which(Gamma2[22:26]!=0)
  g2<-gg[,idx]
  idx<-which(Gamma3[22:26]!=0)
  g3<-gg[,idx]
  
  #load("init.est_1224.RData")
  
  idx<-which(Gamma1!=0)
  Gamma1.selected<-Gamma1[idx]#init.est$Gamma1.est[idx]#
  idx<-which(Gamma2!=0)
  Gamma2.selected<-Gamma2[idx]#init.est$Gamma2.est[idx]#
  idx<-which(Gamma3!=0)
  Gamma3.selected<-Gamma3[idx]#init.est$Gamma3.est[idx]#
  
  
  idx<-which(diag(Sigma1[1:18,1:18])>0.015)
  z1<-z[,idx]
  idx<-which(diag(Sigma1[19:36,19:36])>0.015)
  z2<-z[,idx]
  idx<-which(diag(Sigma1[37:54,37:54])>0.015)
  z3<-z[,idx]
  
  idx<-which(diag(Sigma1)>0.015)
  
  Sigma1.selected<-Sigma1[idx,]
  Sigma1.selected<-Sigma1.selected[,idx]
  
  
  return(list('g1'=g1,'z1'=z1,'x1'=x1,
              'g2'=g2,'z2'=z2,'x2'=x2,
              'g3'=g3,'z3'=z3,'x3'=x3,
              'y'=y,
              'Gamma.1'=Gamma1.selected,
              'Gamma.2'=Gamma2.selected,
              'Gamma.3'=Gamma3.selected,
              'Sigma.1'=Sigma1.selected,
              'Sigma.2'=Sigma2))
}

empbayes_rpql_0720<-function(B,yvec,evec,xmat,zmat,psi1,psi2,psi3,
                             Sigma2inv,Sigma3inv,sig1,s,e,q1,q2,q3){
  
  nD = dim(B)[1]
  
  Xb1 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi1)
  Xb2 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi2)
  Xb3 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi3)
  
  act = matrix(0,nD,(e-s+1))
  time = matrix(0,nD,(e-s+1))
  revgen = matrix(0,nD,(e-s+1))
  fb = matrix(0,nD,(e-s+1))
  
  alphavec = evec
  delvec = yvec
  delvec[delvec>0]<-1
  Yvec_1 = log(yvec)
  Yvec_1[Yvec_1==-Inf]<-0
  
  for (dd in 1:nD){
    
    b = B[dd,]
    
    temp = Xb1+as.matrix(zmat[(s-1):(e-1),])%*%b[1:18];
    f0 = -1*log(1+exp(temp));
    f1 = delvec[(s-1):(e-1)]*temp + f0;
    act[dd,] = cumsum(f1);
    
    mu = (Yvec_1[(s-1):(e-1)]-Xb2-as.matrix(zmat[(s-1):(e-1),])%*%b[19:36])^2;
    g2 = delvec[(s-1):(e-1)]*(-Yvec_1[(s-1):(e-1)]-0.5*(1/sig1)*mu)
    time[dd,] = cumsum(g2);
    
    temp = Xb3+as.matrix(zmat[(s-1):(e-1),])%*%b[37:54]
    f0 = -1*log(1+exp(temp));
    f1 = delvec[(s-1):(e-1)]*alphavec[(s-1):(e-1)]*temp + delvec[(s-1):(e-1)]*f0;
    revgen[dd,] = cumsum(f1);
    
    fb[dd,] = as.vector(-0.5*(t(b[36+q3])%*%Sigma3inv%*%b[36+q3])+
                          1-0.5*(t(b[18+q2])%*%Sigma2inv%*%b[18+q2]))*matrix(1,1,e-s+1);
    
    
  }
  out = act+time+revgen+fb;
  b_hat = matrix(0,e-s+1,dim(B)[2]);
  for (j in 1:(e-s+1)){
    idx = min(which(out[,j]==max(out[,j])))
    b_hat[j,] = B[idx,]
  }
  return(b_hat)
}

empbayes_ccrejm<-function(B,C,C1,C2,C3,yvec,evec,xmat,zmat,psi1,psi2,psi3,
                        Sigma1inv,Sigma2inv,sig1,s,e,q1,q2,q3){
  
  nD = dim(B)[1]
  K<-dim(C)[2]
  
  Xb1 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi1)
  Xb2 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi2)
  Xb3 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi3)
  
  act = matrix(0,nD,(e-s+1))
  time = matrix(0,nD,(e-s+1))
  revgen = matrix(0,nD,(e-s+1))
  fb = matrix(0,nD,(e-s+1))
  fc = matrix(0,nD,(e-s+1))
  
  alphavec = evec
  delvec = yvec
  delvec[delvec>0]<-1
  Yvec_1 = log(yvec)
  Yvec_1[Yvec_1==-Inf]<-0
  
  for (dd in 1:nD){
    
    b = B[dd,]
    
    temp = Xb1+as.matrix(zmat[(s-1):(e-1),])%*%(b[1:18])+C1[(s-1):(e-1),dd];
    f0 = -1*log(1+exp(temp));
    f1 = delvec[(s-1):(e-1)]*temp + f0;
    act[dd,] = cumsum(f1);
    
    mu = ((1/sqrt(sig1))*Yvec_1[(s-1):(e-1)]-Xb2-(1/sqrt(sig1))*as.matrix(zmat[(s-1):(e-1),])%*%b[19:36]-
            (1/sqrt(sig1))*C2[(s-1):(e-1),dd])^2;
    g2 = delvec[(s-1):(e-1)]*(-Yvec_1[(s-1):(e-1)]-0.5*mu)
    time[dd,] = cumsum(g2);
    
    temp = Xb3+as.matrix(zmat[(s-1):(e-1),])%*%b[37:54]+C3[(s-1):(e-1),dd]
    f0 = -1*log(1+exp(temp));
    f1 = delvec[(s-1):(e-1)]*alphavec[(s-1):(e-1)]*temp + delvec[(s-1):(e-1)]*f0;
    revgen[dd,] = cumsum(f1);
    
    bb<-B[dd,c(q1,18+q2,36+q3)]
    
    fb[dd,] = as.vector(-0.5*(t(bb)%*%Sigma1inv%*%bb))*matrix(1,1,e-s+1)
    
    fc[dd,] = (sum(sapply(1:K, function(i) -t(C[i,])%*%Sigma2inv%*%C[i,])))*matrix(1,1,e-s+1)
    
  }
  out = act+time+revgen+fb+fc;
  b_hat = matrix(0,e-s+1,dim(B)[2]);
  C1_hat = matrix(0,e-s+1,1);
  C2_hat = matrix(0,e-s+1,1);
  C3_hat = matrix(0,e-s+1,1);
  for (j in 1:(e-s+1)){
    idx = which(out[,j]==max(out[,j]))
    b_hat[j,] = B[idx,]
    C1_hat[j]<-(C1[(s-1):(e-1),])[j,idx]
    C2_hat[j]<-(C2[(s-1):(e-1),])[j,idx]
    C3_hat[j]<-(C3[(s-1):(e-1),])[j,idx]
  }
  return(list("B"=b_hat,"C1"=C1_hat,"C2"=C2_hat,"C3"=C3_hat))
}

empbayes_ccrejm_guildcorr<-function(B1,B2,C,C1,C2,C3,
                          yvec1,evec1,xmat1,zmat1,yvec2,evec2,xmat2,zmat2,
                          psi1,psi2,psi3,Sigma1inv,Sigma2inv,
                          sig1,q1,q2,q3,tt){
  nD = dim(B1)[1]
  K<-dim(C)[2]
  
  if(tt>2){

  Xb11 = as.matrix(xmat1)%*%as.vector(psi1)
  Xb21 = as.matrix(xmat1)%*%as.vector(psi2)
  Xb31 = as.matrix(xmat1)%*%as.vector(psi3)

  Xb12 = as.matrix(xmat2)%*%as.vector(psi1)
  Xb22 = as.matrix(xmat2)%*%as.vector(psi2)
  Xb32 = as.matrix(xmat2)%*%as.vector(psi3)
  }else{
    Xb11 = xmat1%*%as.vector(psi1)
    Xb21 = xmat1%*%as.vector(psi2)
    Xb31 = xmat1%*%as.vector(psi3)
    
    Xb12 = xmat2%*%as.vector(psi1)
    Xb22 = xmat2%*%as.vector(psi2)
    Xb32 = xmat2%*%as.vector(psi3)
  }

  act = matrix(0,nD,1)
  time = matrix(0,nD,1)
  revgen = matrix(0,nD,1)
  fb = matrix(0,nD,1)
  fc = matrix(0,nD,1)

  alphavec1 = evec1
  delvec1 = yvec1
  delvec1[delvec1>0]<-1
  Yvec_11 = log(yvec1)
  Yvec_11[Yvec_11==-Inf]<-0
  alphavec2 = evec2
  delvec2 = yvec2
  delvec2[delvec2>0]<-1
  Yvec_12 = log(yvec2)
  Yvec_12[Yvec_12==-Inf]<-0

  for (dd in 1:nD){

    b1 = B1[dd,]
    b2 = B2[dd,]

    temp = Xb11+matrix(zmat1,nrow=tt-1)%*%(b1[1:18])+C1[,dd];
    f0 = -1*log(1+exp(temp));
    f11 = delvec1*temp + f0;
    temp = Xb12+matrix(zmat2,nrow=tt-1)%*%(b2[1:18])+C1[,dd];
    f0 = -1*log(1+exp(temp));
    f12 = delvec2*temp + f0;
    act[dd] = sum(f11)+sum(f12);

    mu = ((1/sqrt(sig1))*Yvec_11-Xb21-(1/sqrt(sig1))*matrix(zmat1,nrow=tt-1)%*%b1[19:36]-
            (1/sqrt(sig1))*C2[,dd])^2;
    g21 = delvec1*(-Yvec_11-0.5*mu)
    mu = ((1/sqrt(sig1))*Yvec_12-Xb22-(1/sqrt(sig1))*matrix(zmat2,nrow=tt-1)%*%b2[19:36]-
            (1/sqrt(sig1))*C2[,dd])^2;
    g22 = delvec2*(-Yvec_12-0.5*mu)
    time[dd] = sum(g21)+sum(g22);

    temp = Xb31+matrix(zmat1,nrow=tt-1)%*%b1[37:54]+C3[,dd]
    f0 = -1*log(1+exp(temp));
    f11 = delvec1*alphavec1*temp + delvec1*f0;
    temp = Xb32+matrix(zmat2,nrow=tt-1)%*%b2[37:54]+C3[,dd]
    f0 = -1*log(1+exp(temp));
    f12 = delvec2*alphavec2*temp + delvec2*f0;
    revgen[dd] = sum(f11)+sum(f12);

    bb1<-B1[dd,c(q1,18+q2,36+q3)]
    bb2<-B2[dd,c(q1,18+q2,36+q3)]

    fb[dd] = -0.5*(t(bb1)%*%Sigma1inv%*%bb1)-0.5*(t(bb2)%*%Sigma1inv%*%bb2)

    fc[dd] = (sum(sapply(1:K, function(i) -t(C[i,])%*%Sigma2inv%*%C[i,])))

  }
  out = act+time+revgen+fb+fc;
  idxx = which(out==max(out))
  b1_hat = B1[idxx,]
  b2_hat = B2[idxx,]
  C1_hat<-C1[,idxx]
  C2_hat<-C2[,idxx]
  C3_hat<-C3[,idxx]
  return(list("B1"=b1_hat,"B2"=b2_hat,"C1"=C1_hat,"C2"=C2_hat,"C3"=C3_hat))
}

empbayes_ccrejm_noguildref<-function(B,yvec,evec,xmat,zmat,psi1,psi2,psi3,
                          Sigma1inv,sig1,s,e,q1,q2,q3){
  
  nD = dim(B)[1]
 
  Xb1 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi1)
  Xb2 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi2)
  Xb3 = as.matrix(xmat[(s-1):(e-1),])%*%as.vector(psi3)
  
  act = matrix(0,nD,(e-s+1))
  time = matrix(0,nD,(e-s+1))
  revgen = matrix(0,nD,(e-s+1))
  fb = matrix(0,nD,(e-s+1))
  fc = matrix(0,nD,(e-s+1))
  
  alphavec = evec
  delvec = yvec
  delvec[delvec>0]<-1
  Yvec_1 = log(yvec)
  Yvec_1[Yvec_1==-Inf]<-0
  
  for (dd in 1:nD){
    
    b = B[dd,]
   
    temp = Xb1+as.matrix(zmat[(s-1):(e-1),])%*%(b[1:18])
    f0 = -1*log(1+exp(temp));
    f1 = delvec[(s-1):(e-1)]*temp + f0;
    act[dd,] = cumsum(f1);
    
    mu = ((1/sqrt(sig1))*Yvec_1[(s-1):(e-1)]-Xb2-
            (1/sqrt(sig1))*as.matrix(zmat[(s-1):(e-1),])%*%b[19:36])^2;
    g2 = delvec[(s-1):(e-1)]*(-Yvec_1[(s-1):(e-1)]-0.5*mu)
    time[dd,] = cumsum(g2);
    
    temp = Xb3+as.matrix(zmat[(s-1):(e-1),])%*%b[37:54]
    f0 = -1*log(1+exp(temp));
    f1 = delvec[(s-1):(e-1)]*alphavec[(s-1):(e-1)]*temp + delvec[(s-1):(e-1)]*f0;
    revgen[dd,] = cumsum(f1);
    
    bb<-B[dd,c(q1,18+q2,36+q3)]
    
    fb[dd,] = as.vector(-0.5*(t(bb)%*%Sigma1inv%*%bb))*matrix(1,1,e-s+1)
    
  }
  out = act+time+revgen+fb
  b_hat = matrix(0,e-s+1,dim(B)[2]);
  for (j in 1:(e-s+1)){
    idx = which(out[,j]==max(out[,j]))
    b_hat[j,] = B[idx,]
  }
  return(list("B"=b_hat))
}



#1. Main script for parameter estimation using MCEM after selection
dccjm_postselection_est<-function(y,evec,x1,x2,x3,z1,z2,z3,gg1,gg2,gg3,
                    Gamma1.0,Gamma2.0,Gamma3.0,
                    Sigma1.0,Sigma2.0,guild_mem,
                    t.cut,n,D=2000,itermax){
  
  # x: n(m-1) by pf+pc design matrix 
  # z: n(m-1) by pc design matrix
  # g: K(m-1) by pg design matrix: first (m-1) in guild 1, 2nd (m-1) in guild 2 and so on
  # y: response n(m-1) vector that includes 0
  # evec: response n(m-1) vector for purchase indicator
  # guild.mem: 0, 1 matrix of dim n(m-1) by K
  
  
  #housekeeping
  lQ <- numeric()
  iter<-0
  err<-1
  tol<-0.001
  
  pc1<-dim(z1)[2]
  pc2<-dim(z2)[2]
  pc3<-dim(z3)[2]
  pg1<-dim(gg1)[2]
  pg2<-dim(gg2)[2]
  pg3<-dim(gg3)[2]
  m<-(dim(x1)[1]/n)+1
  K<- dim(guild_mem)[2]
  
  #4. Get the banded penalty matrix for Sigma2
  T<-m-1
  s<-matrix(0,T,T)
  for(i in 1:T){
    s[i,]<-1*(abs(i-c(1:T))<=t.cut)
  }
  e1<-cbind(s,s,s)
  e<-rbind(e1,e1,e1)
  P2 <- matrix(1, 3*T, 3*T)-e
  
  # Create G1,G2,G3 matrix first n(m-1) by pg
  G1<- matrix(0,n*(m-1),pg1)
  G2<- matrix(0,n*(m-1),pg2)
  G3<- matrix(0,n*(m-1),pg3)
  
 for(i in 1:n){
    
    s<-1+(i-1)*(m-1)
    e<-(m-1)*i
    temp1<-matrix(rep(c(guild_mem[s:e,]),pg1),ncol=pg1,byrow = F)
    temp2<-temp1*gg1
    idx<-which(temp2!=0,arr.ind = T)
    G1[s:e,]<-matrix(gg1[unique(idx[,1]),],ncol=pg1,byrow = F)
    temp1<-matrix(rep(c(guild_mem[s:e,]),pg2),ncol=pg2,byrow = F)
    temp2<-temp1*gg2
    idx<-which(temp2!=0,arr.ind = T)
    G2[s:e,]<-matrix(gg2[unique(idx[,1]),],ncol=pg2,byrow = F)
    
    temp1<-matrix(rep(c(guild_mem[s:e,]),pg3),ncol=pg3,byrow = F)
    temp2<-temp1*gg3
    idx<-which(temp2!=0,arr.ind = T)
    G3[s:e,]<-matrix(gg3[unique(idx[,1]),],ncol=pg3,byrow = F)
  }
  
  
  while(iter<=itermax){
    
    if(iter==0){
      out.lQ<-lQ_eval_est(y,evec,x1,x2,x3,z1,z2,z3,G1,G2,G3,
                          Gamma1.0,Gamma2.0,Gamma3.0,
                          Sigma1.0,Sigma2.0,
                          guild_mem,pc1,pc2,pc3,K,n,D,iter)
      lQ[iter+1] <- out.lQ$lQ
      cat("Initial Penalized Q fucntion value: ",lQ[iter+1],"\n",sep="")
      B1<-out.lQ$B1
      B2<-out.lQ$B2 
      B3<-out.lQ$B3
      B<-out.lQ$B #nD by pc1+pc2+pc3 matrix of player random effects for models 1 and 2
      C1<-out.lQ$C1 #n(m-1) by D matrix of Guild random effects for model 1
      C2<-out.lQ$C2 #n(m-1) by D matrix of Guild random effects for model 2
      C3<-out.lQ$C3 #n(m-1) by D matrix of Guild random effects for model 3
      C<-out.lQ$C #KD by 2(m-1) matrix of full Guild random effects
      w<- out.lQ$w
    }else{
      
      B1<-out.lQ$B1
      B2<-out.lQ$B2 
      B3<-out.lQ$B3
      B<-out.lQ$B  #nD by 2pc matrix of player random effects for models 1 and 2
      C1<-out.lQ$C1#n(m-1) by D matrix of Guild random effects for model 1
      C2<-out.lQ$C2#n(m-1) by D matrix of Guild random effects for model 2
      C3<-out.lQ$C3 #n(m-1) by D matrix of Guild random effects for model 3
      C<-out.lQ$C #KD by 2(m-1) matrix of full Guild random effects
      w<- out.lQ$w
    }
    gc()
    #Solve the optimization problems. Get new estimates here.
    
    Gamma1.new<- get_gdest_cpp(y,evec,x1,z1,G1,Gamma1.0,w,B1,C1,D,type=1,tol=0.01)
    
    Gamma2.new<- get_gdest_cpp(y,evec,x2,z2,G2,Gamma2.0,w,B2,C2,D,type=2,tol=0.01)
    
    Gamma3.new<- get_gdest_cpp(y,evec,x3,z3,G3,Gamma3.0,w,B3,C3,D,type=3,tol=0.01)
    
    Sigma1.new<-get_Sigma1_est(Sigma1.0,B,w,n,m,K,D,pc1+pc2+pc3,type=1)
    
    Sigma2.new<-get_Sigma2_est(Sigma2.0,P2,0.1,C,w,K,D)
    
    #Re-Calculate likelihood and error
    iter = iter+1
    out.lQ<-lQ_eval_est_new(y,evec,x1,x2,x3,z1,z2,z3,G1,G2,G3,
                        Gamma1.new,Gamma2.new,Gamma3.new,
                        Sigma1.new,
                        Sigma2.new,guild_mem,
                        w,B1,B2,B3,B,C1,C2,C3,C,
                        pc1,pc2,pc3,K,n,D,iter)
    lQ[iter+1] <- out.lQ$lQ
    cat("Iteration [",iter,"] Q fucntion: ",lQ[iter+1],"\n",sep="")
    err<-lQ[iter+1]-lQ[iter]
    
    Gamma1.0 = Gamma1.new
    Gamma2.0 = Gamma2.new
    Gamma3.0 = Gamma3.new
    Sigma1.0 = Sigma1.new
    Sigma2.0 = Sigma2.new
  
    
  }
  return(list('niter'=iter,'lQ'=lQ,'Gamma1.est'=Gamma1.0,
              'Gamma2.est'=Gamma2.0,'Gamma3.est'=Gamma3.0,
              'Sigma1.est'=Sigma1.0,
              'Sigma2.est'=Sigma2.0))
  
}

#2. Sigma.1 
get_Sigma1_est<-function(Sigma.init,randeff,wmat,n,m,K,D,pc,type){
  
  if(type==1){
    cust_grp<- matrix(rep(1:n,D),n*D,1)-1
    AA<-get_AA(n,pc,D,wmat,randeff,cust_grp)
  }else{
    guild_grp<- matrix(rep(1:K,D),K*D,1)-1
    AA<-get_CC(K,m-1,D,wmat,randeff,guild_grp)
  }
  AA<-polish_Sigma_est(AA)

  return(AA)
}

#3. call spcov with penalty for Sigma.2
get_Sigma2_est<-function(Sigma.init,P,lambda,randeff,wmat,K,D){
  
  #type 2: guild x time rand eff cov matrix
  mm<-(dim(randeff)[2])/3
  guild_grp<- matrix(rep(1:K,D),K*D,1)-1
  AA<-get_CC(K,mm,D,wmat,randeff,guild_grp)
  AA<-polish_Sigma_est(AA)
  
  idx = lambda*P;
  idx[idx>{10^6}]<-10^6
  Rho = idx;
  Sigma_out<-spcov(Sigma.init, AA, Rho, step.size=100)$Sigma
  Sigma_out<-0.5*(Sigma_out+t(Sigma_out))
  
  return(Sigma_out)
}

#4. Calculate the Q function (estimation) 
lQ_eval_est_new<-function(y,evec,x1,x2,x3,z1,z2,z3,G1,G2,G3,
                      Gamma1.est,Gamma2.est,Gamma3.est,
                      Sigma1.est,Sigma2.est,
                      guild.mem,w_old,B1_old,B2_old,B3_old,B_old,
                      C1_old,C2_old,C3_old,C_old,
                      pc1,pc2,pc3,K,n,D,iter){
  
  #Step1: sample n*D multivariate normal variables using the current estimates
  set.seed(iter) 
  B = sample_re(n*D,Sigma1.est);
  B1 = B[,1:pc1]
  B2 = B[,(pc1+1):(pc1+pc2)]
  B3 = B[,(pc1+pc2+1):(pc1+pc2+pc3)]
  
  #Step2: sample K*D multivariate normal variables using the current estimates
  set.seed(iter) 
  C = sample_re(K*D,Sigma2.est);
  C0 = process_guildre(C,guild.mem,n,K,D)
  ii<-1:(n*(m-1))
  C1 = C0[ii,]
  C2 = C0[ii+(n*(m-1)),]
  C3 = C0[ii+(2*n*(m-1)),]
  
  #Step3 : calculate wmat using the current estimates
  wmat<-get_wmat_cpp(y,evec,x1,x2,x3,z1,z2,z3,G1,G2,G3,
                     Gamma1.est,Gamma2.est,Gamma3.est,
                     B,C1,C2,C3,n,D)
  wmat[is.na(wmat)] = 0;
  rr = rowsums(wmat);
  w = wmat/rr;
  w[is.na(w)] = 0;
  
  #Step3: evaluate the Q function
  lq.value<- -nll_cpp(y,evec,x1,z1,G1,Gamma1.est,w_old,
              B1_old,C1_old,D,1)-nll_cpp(y,evec,x2,z2,G2,Gamma2.est,w_old,
              B2_old,C2_old,D,2)-
    nll_cpp(y,evec,x3,z3,G3,Gamma3.est,w_old,B3_old,C3_old,D,3)+
    gll_est(Sigma1.est,Sigma2.est,w_old,B_old,C_old,K)  
  
  return(list('lQ'=lq.value,'w'=w,'B1'=B1,'B2'=B2,'B3'=B3,
              'B'=B,'C1'=C1,'C2'=C2,'C3'=C3,'C'=C))
  
}

lQ_eval_est<-function(y,evec,x1,x2,x3,z1,z2,z3,G1,G2,G3,
                      Gamma1.est,Gamma2.est,Gamma3.est,
                          Sigma1.est,Sigma2.est,
                          guild.mem,pc1,pc2,pc3,K,n,D,iter){
  
  #Step1: sample n*D multivariate normal variables using the current estimates
  set.seed(iter) 
  B = sample_re(n*D,Sigma1.est);
  B1 = B[,1:pc1]
  B2 = B[,(pc1+1):(pc1+pc2)]
  B3 = B[,(pc1+pc2+1):(pc1+pc2+pc3)]
  
  #Step2: sample K*D multivariate normal variables using the current estimates
  set.seed(iter) 
  C = sample_re(K*D,Sigma2.est);
  C0 = process_guildre(C,guild.mem,n,K,D)
  ii<-1:(n*(m-1))
  C1 = C0[ii,]
  C2 = C0[ii+(n*(m-1)),]
  C3 = C0[ii+(2*n*(m-1)),]
  
  #Step3 : calculate wmat using the current estimates
  wmat<-get_wmat_cpp(y,evec,x1,x2,x3,z1,z2,z3,G1,G2,G3,
                     Gamma1.est,Gamma2.est,Gamma3.est,
                     B,C1,C2,C3,n,D)
  wmat[is.na(wmat)] = 0;
  rr = rowsums(wmat);
  w = wmat/rr;
  w[is.na(w)] = 0;
  
  #Step3: evaluate the Q function
  lq.value<- -nll_cpp(y,evec,x1,z1,G1,Gamma1.est,w,
                      B1,C1,D,1)-nll_cpp(y,evec,x2,z2,G2,Gamma2.est,w,
                                                B2,C2,D,2)-
    nll_cpp(y,evec,x3,z3,G3,Gamma3.est,w,B3,C3,D,3)+
    gll_est(Sigma1.est,Sigma2.est,w,B,C,K)  
  
  return(list('lQ'=lq.value,'w'=w,'B1'=B1,'B2'=B2,'B3'=B3,
              'B'=B,'C1'=C1,'C2'=C2,'C3'=C3,'C'=C))
  
}

#5. Gaussian Log Likelihood for the Q function
gll_est<-function(Sigma1,Sigma2,wmat,randeff1,randeff2,K){
  gg<-dmvnorm(randeff1,mean=rep(0,dim(Sigma1)[1]),sigma = Sigma1,log = TRUE)
  gll.1<-matrix(gg,ncol=dim(wmat)[2],byrow=F)*wmat
  gg<-dmvnorm(randeff2,mean=rep(0,dim(Sigma2)[1]),sigma = Sigma2,log = TRUE)
  w.d<-matrix(rep(colmeans(wmat),each=K),K*D,1)
  gll.2<-matrix(gg,ncol=1)*w.d
  return(sum(gll.1)+sum(gll.2))
}

#6. Polishing Sigma to ensure numerical stability
polish_Sigma_est<-function(Sigma_given){
  
  Sigma_clean<-Sigma_given
  min.eig <- min(eigen(Sigma_clean)$values)
  while(min.eig<10^{-2}){
    
    Sigma_clean = Sigma_clean+10^(-3)*diag(dim(Sigma_clean)[2]);
    min.eig <- min(eigen(Sigma_clean)$values)
  }
  return(Sigma_clean)
}



