
library(Rfast)
library(msos)
library(mvtnorm)
library(CVXR)


#---------------------- For simulation------------------------
get_simdata_est<-function(n,K,pf,pc,pg,t.cut,guild_mem,
                          Gamma.1,Gamma.2,Gamma.3,sigma,Sigma.1,
                          Sigma.2){
  
  
  pc<-(dim(Sigma.1)[2])/3
  m<-(dim(Sigma.2)[2])/3+1
  Sigma.11<-Sigma.1[1:pc,1:pc]
  Sigma.12<-Sigma.1[(pc+1):(2*pc),(pc+1):(2*pc)]
  Sigma.13<-Sigma.1[(2*pc+1):(3*pc),(2*pc+1):(3*pc)]
  #Random Effects
  set.seed(1)
  bi<-rmvnorm(n,rep(0,3*pc),Sigma.1)
  set.seed(2)
  ck<-rmvnorm(K,rep(0,3*(m-1)),Sigma.2)
  
  # Design matrices
  set.seed(1)
  x<-matrix(rnorm(n*(m-1)*(pc+pf)),n*(m-1),pf+pc)#matrix(runif(n*(m-1)*pc),n*(m-1),pc)
  set.seed(2)
  g<-matrix(rnorm(K*(m-1)*pg),K*(m-1),pg)#matrix(runif(n*pf),n,pf)
  z<-x[,(pf+1):(pf+pc)]
  
  # Create G matrix first n(m-1) by pg
  G<- matrix(0,n*(m-1),pg)
  for(i in 1:n){
    
    s<-1+(i-1)*(m-1)
    e<-(m-1)*i
    temp1<-matrix(rep(c(guild_mem[s:e,]),pg),ncol=pg,byrow = F)
    temp2<-temp1*g
    idx<-(temp2!=0)
    G[s:e,]<-matrix(g[idx],ncol=pg,byrow = F)
  }
  x_g<-cbind(x,G)
  
  #------------
  ci.1 = matrix(0,n*(m-1),1);
  ci.2 = matrix(0,n*(m-1),1);
  ci.3 = matrix(0,n*(m-1),1);
  
  for(i in 1:n){
    s<-1+(i-1)*(m-1)
    e<-(m-1)*i
    ci.1[s:e]<-rowsums(t(ck[,1:(m-1)])*guild_mem[s:e,])
    ci.2[s:e]<-rowsums(t(ck[,m:(2*(m-1))])*guild_mem[s:e,])
    ci.3[s:e]<-rowsums(t(ck[,(2*(m-1)+1):(3*(m-1))])*guild_mem[s:e,])
  }
  #------------------------------------
  
  # Response 1 (login indicator)
  y1.true<-matrix(0,n*(m-1))
  for(i in 1:n){
    
    s<- 1+(i-1)*(m-1)
    e<- i*(m-1)
    xi<- x_g[s:e,]
    zi<-z[s:e,]
    qi.1<-xi%*%Gamma.1+zi%*%bi[i,1:pc]+ci.1[s:e]
    set.seed(i)
    y1.true[s:e]<-rbinom(m-1,1,plogis(qi.1))
  }
  
  # Response 2 (playing time)
  y2.true<-matrix(0,n*(m-1))
  for(i in 1:n){
    s<- 1+(i-1)*(m-1)
    e<- i*(m-1)
    xi<- x_g[s:e,]
    zi<-z[s:e,]
    qi.2<-xi%*%(Gamma.2[1:(pf+pc+pg)]/Gamma.2[pf+pc+pg+1])+
      zi%*%bi[i,(pc+1):(2*pc)]+ci.2[s:e]
    set.seed(i)
    Zi.2<-rmvnorm(1,mean=qi.2,sigma=(sigma^2)*diag(m-1))
    y2.true[s:e]<- y1.true[s:e]*exp(Zi.2)
  }
  
  # Response 3 (purchase indicator)
  e.true<-matrix(0,n*(m-1))
  for(i in 1:n){
    
    s<- 1+(i-1)*(m-1)
    e<- i*(m-1)
    ei<-e.true[s:e]
    y1i<- y1.true[s:e]
    xi<- x_g[s:e,]
    zi<-z[s:e,]
    ci3<-ci.3[s:e]
    idx<-(y1i>0)
    qi.3<-xi[idx,]%*%Gamma.3+zi[idx,]%*%bi[i,(2*pc+1):(3*pc)]+ci3[idx]
    set.seed(i)
    ei[idx]<-rbinom(sum(1*idx),1,plogis(qi.3))
    e.true[s:e]<-ei
  }
  
  return(list('g'=g,'z'=z,'x'=x,'y'=y2.true,'e'=e.true,
              'Gamma1.true'=Gamma.1,
              'Gamma2.true'=Gamma.2,'Gamma3.true'=Gamma.3,
              'sigma'=sigma,
              'Sigma1.true'=Sigma.1,'Sigma2.true'=Sigma.2,
              't.cut'=t.cut))
}

gen_bandedcov<-function(t.val,t.cut,q){
  
  s<-matrix(0,t.val,t.val)
  for(i in 1:t.val){
    s[i,]<-1*(abs(i-c(1:t.val))<=t.cut)
  }
  e1<-cbind(s,s,s)
  e<-rbind(e1,e1,e1)
  
  g<-graph_from_adjacency_matrix(e-diag(diag(e)), mode = "undirected")
  b <- rep(1, 3*t.val)
  Sig <- q*generate_gb_covariance(g, b)
  return(Sig)
}

gen_guildmat<-function(n,t.val,K){
  
  guildmat<-matrix(0,n*t.val,K)
  
  for(i in 1:dim(guildmat)[1]){
    set.seed(i)
    j<-which(rmultinom(1, 1, rep(1/K,K))>0)#sample(K,1)
    guildmat[i,j]<-1
  }
  return(guildmat)
}

posdef <- function (n, ev = runif(n, 0, 1)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

#------------------ Modelling scripts -----------------------------

#1.1. Main script for parameter estimation using MCEM
dccjm_est<-function(y,evec,x,z,gg,Gamma1.0,Gamma2.0,Gamma3.0,
                    Sigma1.0,Sigma2.0,guild_mem,
                    t.cut,n,D=2000,itermax){
  
  # x: n(m-1) by pf+pc design matrix 
  # z: n(m-1) by pc design matrix
  # g: K(m-1) by pg design matrix: first (m-1) in guild 1, 2nd (m-1) in guild 2 and so on
  # y: response n(m-1) vector that includes 0
  # e: response n(m-1) vector for purchase indicator
  # guild.mem: 0, 1 matrix of dim n(m-1) by K
  
  
  #housekeeping
  lQ <- numeric()
  iter<-0
  err<-1
  tol<-0.001
  
  #pfpc<-dim(x)[2]
  pc<-dim(z)[2]
  pg<-dim(gg)[2]
  m<-(dim(x)[1]/n)+1
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
  
  while(iter<=itermax){
    
    if(iter==0){
      out.lQ<-lQ_eval_est(y,evec,x,z,G,Gamma1.0,Gamma2.0,Gamma3.0,
                          Sigma1.0,Sigma2.0,
                          guild_mem,K,n,D,iter)
      lQ[iter+1] <- out.lQ$lQ
      cat("Initial Penalized Q fucntion value: ",lQ[iter+1],"\n",sep="")
      B<-out.lQ$B #nD by 3pc matrix of player random effects for models 1 and 2
      C1<-out.lQ$C1 #n(m-1) by D matrix of Guild random effects for model 1
      C2<-out.lQ$C2 #n(m-1) by D matrix of Guild random effects for model 2
      C3<-out.lQ$C3 #n(m-1) by D matrix of Guild random effects for model 3
      C<-out.lQ$C #KD by 2(m-1) matrix of full Guild random effects
      w<- out.lQ$w
    }else{
      
      B<-out.lQ$B  #nD by 2pc matrix of player random effects for models 1 and 2
      C1<-out.lQ$C1#n(m-1) by D matrix of Guild random effects for model 1
      C2<-out.lQ$C2#n(m-1) by D matrix of Guild random effects for model 2
      C3<-out.lQ$C3 #n(m-1) by D matrix of Guild random effects for model 3
      C<-out.lQ$C #KD by 2(m-1) matrix of full Guild random effects
      w<- out.lQ$w
    }
    gc()
    #Solve the optimization problems. Get new estimates here.
    
    Gamma1.new<- get_gdest_cpp(y,evec,x,z,G,Gamma1.0,w,B,C1,D,type=1,tol=0.1)
    
    Gamma2.new<- get_gdest_cpp(y,evec,x,z,G,Gamma2.0,w,B,C2,D,type=2,tol=0.01)
    
    Gamma3.new<- get_gdest_cpp(y,evec,x,z,G,Gamma3.0,w,B,C3,D,type=3,tol=0.1)
    
    Sigma1.new<-get_Sigma_est(Sigma1.0,B,w,n,m,K,D,pc,type=1)
   
    Sigma2.new<-get_Sigma_est(Sigma2.0,C,w,n,m,K,D,pc,type=2)
    
    #Re-Calculate likelihood and error
    iter = iter+1
    out.lQ<-lQ_eval_est_new(y,evec,x,z,G,Gamma1.new,Gamma2.new,Gamma3.new,
                            Sigma1.new,
                            Sigma2.new,guild_mem,
                            w,B,C1,C2,C3,C,
                            K,n,D,iter)
    lQ[iter+1] <- out.lQ$lQ
    cat("Iteration [",iter,"] Q fucntion: ",lQ[iter+1],"\n",sep="")
    err<-lQ[iter+1]-lQ[iter]
    
    #Refresh estimates
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

#1.2. Main script for selection using MCEM that has a penalized Q function
dccjm_select<-function(y,evec,x,z,gg,Gamma1.0,Gamma2.0,Gamma3.0,
                       Sigma1.0,Sigma2.0,guild.mem,
                       A.1,a.1,A.2,a.2,A.3,a.3,
                       n,t.cut,lam1,lam2,lam3,D=2000,itermax){
  
  # x: n(m-1) by pf+pc design matrix 
  # z: n(m-1) by pc design matrix
  # gg: K(m-1) by pg design matrix: first (m-1) in guild 1, 2nd (m-1) in guild 2 and so on
  # y: response n(m-1) vector that includes 0
  # e: response n(m-1) vector for purchase indicator
  # guild.mem: 0, 1 matrix of dim n(m-1) by K
  # A.1 : q1 by (pf+pc+pg) constraint matrix
  # A.2 : q2 by (pf+pc+pg+1) constraint matrix
  # A.3 : q3 by (pf+pc+pg) constraint matrix
  # t.cut: time cut off for banded Sigma2. Must be <=(m-1)
  
  #housekeeping
  lQ <- numeric()
  iter<-0
  err<-1
  tol<-0.001
  
  pfpc<-dim(x)[2]
  pc<-dim(z)[2]
  pf<-pfpc-pc
  pg<-dim(gg)[2]
  m<-(dim(x)[1]/n)+1
  K<- dim(guild.mem)[2]
  
  # Create G matrix first n(m-1) by pg
  G<- matrix(0,n*(m-1),pg)
  for(i in 1:n){
    
    s<-1+(i-1)*(m-1)
    e<-(m-1)*i
    temp1<-matrix(rep(c(guild.mem[s:e,]),pg),ncol=pg,byrow = F)
    temp2<-temp1*gg
    idx<-which(temp2!=0,arr.ind = T)
    G[s:e,]<-matrix(gg[unique(idx[,1]),],ncol=pg,byrow = F)
  }
  
  while(iter<=itermax){
    
    #update the weight matrices
    mats<-get_mats(Gamma1.0,Gamma2.0,Gamma3.0,Sigma1.0,t.cut,m,pf,pc)
    W.1<- mats$W1
    W.2<-mats$W2
    W.3<-mats$W3
    P1.mat<- mats$P1
    P2.mat<- mats$P2
    
    if(iter==0){
      out.lQ<-lQ_eval_select_1(y,evec,x,z,G,Gamma1.0,Gamma2.0,Gamma3.0,
                               Sigma1.0,Sigma2.0,
                               guild.mem,W.1,W.2,W.3,P1.mat,P2.mat,lam1,lam2,lam3,
                               K,n,D,iter)
      lQ[iter+1] <- out.lQ$lQ
      cat("Initial Penalized Q fucntion value: ",lQ[iter+1],"\n",sep="")
      B<-out.lQ$B #nD by 3pc matrix of player random effects for models 1,2,3
      C<-out.lQ$C #KD by 3(m-1) matrix of full Guild random effects
      C1<-out.lQ$C1 #n(m-1) by D matrix of Guild random effects for model 1
      C2<-out.lQ$C2 #n(m-1) by D matrix of Guild random effects for model 2
      C3<-out.lQ$C3 #n(m-1) by D matrix of Guild random effects for model 3
      w<- out.lQ$w
    }else{
      
      B<-matrix(0,n*D,3*pc)
      b<-out.lQ$B
      B[,idx_1]<-b[,idx_1]  #nD by 3pc matrix of player random effects for models 1,2,3
      C<-out.lQ$C #KD by 3(m-1) matrix of full Guild random effects
      C1<-out.lQ$C1#n(m-1) by D matrix of Guild random effects for model 1
      C2<-out.lQ$C2#n(m-1) by D matrix of Guild random effects for model 2
      C3<-out.lQ$C3 #n(m-1) by D matrix of Guild random effects for model 3
      w<- out.lQ$w
    }
    
    #Solve the optimization problems. Get new estimates here.
    
    Gamma1.new<- get_Gamma(y,evec,x,z,G,Gamma1.0,w,W.1,A.1,a.1,B,C1,n,D,
                           lam1,type=1,tol=0.01)
    Gamma1.new[abs(Gamma1.new)<0.01]<-0
    
    Gamma2.new<- get_Gamma(y,evec,x,z,G,Gamma2.0,w,W.2,A.2,a.2,B,C2,n,D,
                           lam2,type=2,tol=0.01)
    idx_3<-(abs(Gamma2.new/Gamma2.new[pf+pc+pg+1])<0.01)
    Gamma2.new[idx_3]<-0
    
    Gamma3.new<- get_Gamma(y,evec,x,z,G,Gamma3.0,w,W.3,A.3,a.3,B,C3,n,D,
                           lam1,type=3,tol=0.01)
    Gamma3.new[abs(Gamma3.new)<0.01]<-0
    
    Sigma1.new<-get_Sigma(Sigma1.0,P1.mat,lam1,B,w,n,D,type=1)
    idx_1 = (abs(diag(Sigma1.new))>10^(-2))
    idx_2 = (abs(diag(Sigma1.new))<=10^(-2))
    Sigma1.new[idx_2,] = 10^(-3)
    Sigma1.new[,idx_2] = 10^(-3);
    Sigma1.new<-polish_Sigma_est(Sigma1.new)
    
    Sigma2.new<-get_Sigma(Sigma2.0,P2.mat,lam3,C,w,K,D,type=2)
    Sigma2.new<-polish_Sigma_est(Sigma2.new)
    
    #Re-Calculate likelihood and error
    iter = iter+1
    out.lQ<-lQ_eval_select_2(y,evec,x,z,G,Gamma1.new,Gamma2.new,Gamma3.new,
                             Sigma1.new,
                             Sigma2.new,guild.mem,w,B,C1,C2,C3,C,
                             W.1,W.2,W.3,P1.mat,P2.mat,lam1,lam2,lam3,
                             K,n,D,iter)
    
    lQ[iter+1] <- out.lQ$lQ
    cat("Iteration [",iter,"] Penalized Q function: ",lQ[iter+1],"\n",sep="")
    
    Gamma1.0 = Gamma1.new
    Gamma2.0 = Gamma2.new
    Gamma3.0 = Gamma3.new
    Sigma1.0 = Sigma1.new
    Sigma2.0 = Sigma2.new
    err<-lQ[iter+1]-lQ[iter]
    
  }
  return(list('niter'=iter,'lQ'=lQ,'Gamma1.est'=Gamma1.0,
              'Gamma2.est'=Gamma2.0,'Gamma3.est'=Gamma3.0,
              'Sigma1.est'=Sigma1.0,
              'Sigma2.est'=Sigma2.0,
              'idx_1'=idx_1,'idx2'=idx_2,'lQ.val'=lQ[iter+1]))
  
}

#2. Get Gamma(s) estimates from proximal gradient descent (for selection)
get_Gamma<-function(y,evec,x,z,G,Gamma.k,wmat,W,A,a,B,CC,n,D,lam,type,tol){
  
  Gamma.new<-Gamma.k
  
  #Do proximal gradient descent with line search
  pgd_err<-1
  pgd_iter<-0
  pgd_itermax<-10
  while(pgd_err>tol){
    
    fb = nll_cpp(y,evec,x,z,G,Gamma.k,wmat,B,CC,D,type)
    df.Gamma<-df_Gamma_cpp(y,evec,x,z,G,Gamma.k,wmat,B,CC,D,type)
    df.Gamma<-df.Gamma/(1+norm(df.Gamma,type="2"))
    
    # proximal gradient descent
    lr.new<-1
    U<- Gamma.k-lr.new*df.Gamma
    Gamma.new<-proxgrad(U,W,A,a,lam,lr.new,length(U))
    fb.new<-nll_cpp(y,evec,x,z,G,Gamma.new,wmat,B,CC,D,type)
    temp<- (0.5/lr.new)*sum((Gamma.new-Gamma.k)^2) 
    ee = fb+sum(df.Gamma*(Gamma.new-Gamma.k))+temp
    
    # do line search & update
    while(fb.new-ee>0.01){
      lr.new = 0.2*lr.new;
      U<- Gamma.k-lr.new*df.Gamma
      Gamma.new = proxgrad(U,W,A,a,lam,lr.new,length(U))
      
      #update line search condition
      fb.new<-nll_cpp(y,evec,x,z,G,Gamma.new,wmat,B,CC,D,type)
      temp<- (0.5/lr.new)*sum((Gamma.new-Gamma.k)^2) 
      ee = fb+sum(df.Gamma*(Gamma.new-Gamma.k))+temp
    }
    
    # error update
    pgd_err<-sqrt(sum(Gamma.new-Gamma.k)^2)
    Gamma.k<- Gamma.new
    pgd_iter<-pgd_iter+1
  }
  Gamma.hat<-Gamma.new
  return(Gamma.hat)
}

#3. Get proximal gradient estimates from CVX
proxgrad<-function(U,W,A,a,lam,t,p){
  
  gam <- Variable(p)
  obj<- (0.5/t)*sum_squares(gam-U)+lam*sum(abs(W%*%gam))
  constr <- list(A%*%gam <= a)
  prob <- Problem(Minimize(obj),constr)
  result <- solve(prob)
  gam.out<-result$getValue(gam)
  return(gam.out)
}

#4.1. Sigma.1 and Sigma.2 estimation
get_Sigma_est<-function(Sigma.init,randeff,wmat,n,m,K,D,pc,type){
  
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

#4.2. call spcov with penalty (selection)
get_Sigma<-function(Sigma.init,P,lambda,randeff,wmat,n,D,pc,type){
  
  #type 1: player rand eff cov matrix
  #type 2: guild x time rand eff cov matrix
  if(type==1){
    pc<-(dim(randeff)[2])/3
    cust_grp<- matrix(rep(1:n,D),n*D,1)-1
    AA<-get_AA(n,pc,D,wmat,randeff,cust_grp)
  }else{
    K<-n
    mm<-(dim(randeff)[2])/3
    guild_grp<- matrix(rep(1:K,D),K*D,1)-1
    AA<-get_CC(K,mm,D,wmat,randeff,guild_grp)
  }
  AA<-polish_Sigma_est(AA)
  
  idx = lambda*P;
  idx[idx>{10^6}]<-10^6
  Rho = idx;
  Sigma_out<-spcov(Sigma.init, AA, Rho, step.size=100)$Sigma
  Sigma_out<-0.5*(Sigma_out+t(Sigma_out))
  
  return(Sigma_out)
}

#5. Calculate the Q function (estimation) 
lQ_eval_est_new<-function(y,evec,x,z,G,Gamma1.est,Gamma2.est,Gamma3.est,
                          Sigma1.est,Sigma2.est,
                          guild.mem,w_old,B_old,C1_old,C2_old,C3_old,C_old,
                          K,n,D,iter){
  
  #Step1: sample n*D multivariate normal variables using the current estimates
  set.seed(iter) 
  B = sample_re(n*D,Sigma1.est);
  
  #Step2: sample K*D multivariate normal variables using the current estimates
  set.seed(iter) 
  C = sample_re(K*D,Sigma2.est);
  C0 = process_guildre(C,guild.mem,n,K,D)
  ii<-1:(n*(m-1))
  C1 = C0[ii,]
  C2 = C0[ii+(n*(m-1)),]
  C3 = C0[ii+(2*n*(m-1)),]
  
  #Step3 : calculate wmat using the current estimates
  wmat<-get_wmat_cpp(y,evec,x,z,G,Gamma1.est,Gamma2.est,Gamma3.est,B,C1,C2,C3,n,D)
  wmat[is.na(wmat)] = 0;
  rr = rowsums(wmat);
  w = wmat/rr;
  w[is.na(w)] = 0;
  
  #Step3: evaluate the Q function
  lq.value<- -nll_cpp(y,evec,x,z,G,Gamma1.est,w_old,
                      B_old,C1_old,D,1)-nll_cpp(y,evec,x,z,G,Gamma2.est,w_old,
                                                B_old,C2_old,D,2)-
    nll_cpp(y,evec,x,z,G,Gamma3.est,w_old,B_old,C3_old,D,3)+
    gll_est(Sigma1.est,Sigma2.est,w_old,B_old,C_old,K)  
  
  return(list('lQ'=lq.value,'w'=w,'B'=B,'C1'=C1,'C2'=C2,'C3'=C3,'C'=C))
  
}

lQ_eval_est<-function(y,evec,x,z,G,Gamma1.est,Gamma2.est,Gamma3.est,
                      Sigma1.est,Sigma2.est,
                      guild.mem,K,n,D,iter){
  
  #Step1: sample n*D multivariate normal variables using the current estimates
  set.seed(iter) 
  B = sample_re(n*D,Sigma1.est);
  
  #Step2: sample K*D multivariate normal variables using the current estimates
  set.seed(iter) 
  C = sample_re(K*D,Sigma2.est);
  C0 = process_guildre(C,guild.mem,n,K,D)
  ii<-1:(n*(m-1))
  C1 = C0[ii,]
  C2 = C0[ii+(n*(m-1)),]#C0[-ii,]
  C3 = C0[ii+(2*n*(m-1)),]
  
  #Step3 : calculate wmat using the current estimates
  wmat<-get_wmat_cpp(y,evec,x,z,G,Gamma1.est,Gamma2.est,Gamma3.est,B,C1,C2,C3,n,D)
  wmat[is.na(wmat)] = 0;
  rr = rowsums(wmat);
  w = wmat/rr;
  w[is.na(w)] = 0;
  
  #Step3: evaluate the Q function
  lq.value<- -nll_cpp(y,evec,x,z,G,Gamma1.est,w,
                      B,C1,D,1)-nll_cpp(y,evec,x,z,G,Gamma2.est,w,
                                        B,C2,D,2)-
    nll_cpp(y,evec,x,z,G,Gamma3.est,w,B,C3,D,3)+
    gll_est(Sigma1.est,Sigma2.est,w,B,C,K)  
  
  return(list('lQ'=lq.value,'w'=w,'B'=B,'C1'=C1,'C2'=C2,'C3'=C3,'C'=C))
  
}

#6. Calculate the Q function (selection) 
lQ_eval_select_2<-function(y,evec,x,z,G,Gamma1.est,Gamma2.est,Gamma3.est,
                           Sigma1.est,Sigma2.est,
                           guild.mem,w_old,B_old,C1_old,C2_old,C3_old,C_old,
                           W1,W2,W3,P1,P2,lam1,lam2,lam3,
                           K,n,D,iter){
  
  #Step1: sample n*D multivariate normal variables using the current estimates
  set.seed(iter) 
  B = sample_re(n*D,Sigma1.est);
  
  #Step2: sample K*D multivariate normal variables using the current estimates
  set.seed(iter) 
  C = sample_re(K*D,Sigma2.est);
  C0 = process_guildre(C,guild.mem,n,K,D)
  ii<-1:(n*(m-1))
  C1 = C0[ii,]
  C2 = C0[ii+(n*(m-1)),]#C0[-ii,]
  C3 = C0[ii+(2*n*(m-1)),]
  
  #Step3 : calculate wmat using the current estimates
  wmat<-get_wmat_cpp(y,evec,x,z,G,Gamma1.est,Gamma2.est,Gamma3.est,B,C1,C2,C3,n,D)
  wmat[is.na(wmat)] = 0;
  rr = rowsums(wmat);
  w = wmat/rr;
  w[is.na(w)] = 0;
  
  #Step4: get the penalty value
  # penalty.value<-penalty(Gamma1.est,Gamma2.est,Gamma3.est,Sigma1.est,Sigma2.est,
  #                        W1,W2,W3,P1,P2,lam1,lam2,lam3)
  
  #Step5: evaluate the Penalized Q function
  lq.value<- -nll_cpp(y,evec,x,z,G,Gamma1.est,w_old,
                      B_old,C1_old,D,1)-nll_cpp(y,evec,x,z,G,Gamma2.est,w_old,
                                                B_old,C2_old,D,2)-
    nll_cpp(y,evec,x,z,G,Gamma3.est,w_old,B_old,C3_old,D,3)+
    gll_est(Sigma1.est,Sigma2.est,w_old,B_old,C_old,K)#-penalty.value  
  
  return(list('lQ'=lq.value,'w'=w,'B'=B,'C1'=C1,'C2'=C2,'C3'=C3,'C'=C))
  
}

lQ_eval_select_1<-function(y,evec,x,z,G,Gamma1.est,Gamma2.est,Gamma3.est,
                           Sigma1.est,Sigma2.est,
                           guild.mem,
                           W1,W2,W3,P1,P2,lam1,lam2,lam3,
                           K,n,D,iter){
  
  #Step1: sample n*D multivariate normal variables using the current estimates
  set.seed(iter) 
  B = sample_re(n*D,Sigma1.est);
  
  #Step2: sample K*D multivariate normal variables using the current estimates
  set.seed(iter) 
  C = sample_re(K*D,Sigma2.est);
  C0 = process_guildre(C,guild.mem,n,K,D)
  ii<-1:(n*(m-1))
  C1 = C0[ii,]
  C2 = C0[ii+(n*(m-1)),]#C0[-ii,]
  C3 = C0[ii+(2*n*(m-1)),]
  
  #Step3 : calculate wmat using the current estimates
  wmat<-get_wmat_cpp(y,evec,x,z,G,Gamma1.est,Gamma2.est,Gamma3.est,B,C1,C2,C3,n,D)
  wmat[is.na(wmat)] = 0;
  rr = rowsums(wmat);
  w = wmat/rr;
  w[is.na(w)] = 0;
  
  #Step4: get the penalty value
  # penalty.value<-penalty(Gamma1.est,Gamma2.est,Gamma3.est,Sigma1.est,Sigma2.est,
  #                        W1,W2,W3,P1,P2,lam1,lam2,lam3)
  
  #Step5: evaluate the Q function
  lq.value<- -nll_cpp(y,evec,x,z,G,Gamma1.est,w,
                      B,C1,D,1)-nll_cpp(y,evec,x,z,G,Gamma2.est,w,
                                        B,C2,D,2)-
    nll_cpp(y,evec,x,z,G,Gamma3.est,w,B,C3,D,3)+
    gll_est(Sigma1.est,Sigma2.est,w,B,C,K)#-penalty.value   
  
  return(list('lQ'=lq.value,'w'=w,'B'=B,'C1'=C1,'C2'=C2,'C3'=C3,'C'=C))
  
}


#6. Gaussian Log Likelihood for the Q function
gll_est<-function(Sigma1,Sigma2,wmat,randeff1,randeff2,K){
  gg<-dmvnorm(randeff1,mean=rep(0,dim(Sigma1)[1]),sigma = Sigma1,log = TRUE)
  gll.1<-matrix(gg,ncol=dim(wmat)[2],byrow=F)*wmat
  gg<-dmvnorm(randeff2,mean=rep(0,dim(Sigma2)[1]),sigma = Sigma2,log = TRUE)
  w.d<-matrix(rep(colmeans(wmat),each=K),K*D,1)
  gll.2<-matrix(gg,ncol=1)*w.d
  return(sum(gll.1)+sum(gll.2))
}

#7. Polishing Sigma to ensure numerical stability
polish_Sigma_est<-function(Sigma_given){
  
  Sigma_clean<-Sigma_given
  min.eig <- min(eigen(Sigma_clean)$values)
  while(min.eig<10^{-2}){
    
    Sigma_clean = Sigma_clean+10^(-3)*diag(dim(Sigma_clean)[2]);
    min.eig <- min(eigen(Sigma_clean)$values)
  }
  return(Sigma_clean)
}

#8. Penalty Matrices
get_mats<-function(Gamma1.k,Gamma2.k,Gamma3.k,Sigma1.k,t.val,m,pf,pc){
  
  #1. Get A1.mat for model 1
  w1<-c((Gamma1.k)^{-2})
  W1.mat<-diag(w1)
  W1.mat[W1.mat>10^3]<-10^3
  
  #2. Get A2.mat for model 2
  nparams<-dim(Gamma2.k)[1]
  w2<-c((Gamma2.k)^{-2})
  w2[nparams]<-0
  W2.mat<-diag(w2)
  W2.mat[W2.mat>10^3]<-10^3
  
  #3. Get A3.mat for model 3
  w3<-c((Gamma3.k)^{-2})
  W3.mat<-diag(w3)
  W3.mat[W3.mat>10^3]<-10^3
  
  #3. Get the penalty matrix for Sigma1
  compeffects<- c(1:14,16:18,21)
  w1<-diag(W1.mat)[compeffects]
  w2<-diag(W2.mat)[compeffects]
  w3<-diag(W3.mat)[compeffects]
  i1<-(c(w1,w2,w3)>=10^3)
  weights.c<-c(w1,w2,w3)*(1/(diag(Sigma1.k))^2)
  weights.c[i1]<-10^6
  P1<-diag(weights.c)
  
  #4. Get the banded penalty matrix for Sigma2
  T<-m-1
  s<-matrix(0,T,T)
  for(i in 1:T){
    s[i,]<-1*(abs(i-c(1:T))<=t.val)
  }
  e1<-cbind(s,s,s)
  e<-rbind(e1,e1,e1)
  P2 <- matrix(1, 3*T, 3*T)-e
  
  return(list('P1'=P1,'P2'=P2,'W1'=W1.mat,'W2'=W2.mat,'W3'=W3.mat))
}

#9. Penalty (to be done in c++)
penalty<-function(Gamma1.est,Gamma2.est,Gamma3.est,Sigma1.est,Sigma2.est,
                  W1,W2,W3,P1,P2,lam1,lam2,lam3){
  
  a1<-W1%*%Gamma1.est
  a2<-W2%*%Gamma2.est
  a3<-W3%*%Gamma3.est
  p1<-lam1*(sum(abs(a1))+(lam2/lam1)*sum(abs(a2))+
              sum(abs(a3)))
  p2<- lam1*sum(abs(diag(P1)*Sigma1.est))+lam3*sum(abs(diag(P2)*Sigma2.est))
  return(p1+p2)
}


