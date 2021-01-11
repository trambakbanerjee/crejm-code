// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
using namespace std;
using namespace Rcpp;
using namespace arma;

// SAMPLING RANDOM EFFECTS FROM ZERO MEAN GAUSSIANS
// [[Rcpp::export]]
mat sample_re(int n, const mat & Sigma) {
  int ncols = Sigma.n_cols;
  mat Y = randn(n, ncols);
  return Y * chol(Sigma);
}

// PROCESS GUILD RANDOM EFFECTS
// [[Rcpp::export]]
mat process_guildre(const mat & C,const mat & guild_mem,
                    int n,int K,int D) {
  int ncols = C.n_cols;
  int a = ncols/3;
  mat C_1 = zeros<mat>(n*a,D);
  mat C_2 = zeros<mat>(n*a,D);
  mat C_3 = zeros<mat>(n*a,D);
  
  IntegerVector r_1 = rep(seq(0,a-1),n);
  IntegerVector r_2 = seq(0,K-1);
  uvec a_1 = as<uvec>(r_1);
  uvec a_2 = as<uvec>(r_2);
  
  for(int d=0;d<D;d++){
    int s = d*K;
    int e = K-1+d*K;
    mat tc_1 = trans(C(span(s,e),span(0,a-1)));
    mat tc_2 = trans(C(span(s,e),span(a,2*a-1)));
    mat tc_3 = trans(C(span(s,e),span(2*a,ncols-1)));
    mat t_1 = tc_1(a_1,a_2);
    mat t_2 = tc_2(a_1,a_2);
    mat t_3 = tc_3(a_1,a_2);
    C_1(span::all,d) = sum(t_1%guild_mem,1);
    C_2(span::all,d) = sum(t_2%guild_mem,1);
    C_3(span::all,d) = sum(t_3%guild_mem,1);
  }
  return join_cols(C_1,C_2,C_3);
}

// get AA matrix
// [[Rcpp::export]]
mat get_AA(int n, int pc, int D, const mat & wmat, 
           const mat & randeff,const mat & cust_grp)
{
  
  mat AA = zeros<mat>(3*pc,3*pc);
  mat M = zeros<mat>(D,3*pc);
  mat b1 = zeros<mat>(1,3*pc);
  IntegerVector r_2 = seq(0,3*pc-1);
  uvec a_2 = as<uvec>(r_2);
  for(int i=0; i<n; i++){
    uvec a_1 = find(cust_grp==i);
    M = randeff(a_1, a_2);
    for(int d=0;d<D;d++){
      b1 = M(d,span::all);
      AA = (b1.t()*b1)*wmat(i,d)+AA;
    }
  }
  return AA/n;
}

// get CC matrix
// [[Rcpp::export]]
mat get_CC(int K, int pc, int D, const mat & wmat, 
           const mat & randeff,const mat & guild_grp)
{
  
  mat CC = zeros<mat>(3*pc,3*pc);
  mat M = zeros<mat>(D,3*pc);
  mat b1 = zeros<mat>(1,3*pc);
  IntegerVector r_2 = seq(0,3*pc-1);
  uvec a_2 = as<uvec>(r_2);
  //int n = wmat.n_rows;
  for(int i=0; i<K; i++){
    uvec a_1 = find(guild_grp==i);
    M = randeff(a_1, a_2);
    for(int d=0;d<D;d++){
      b1 = M(d,span::all);
      double ww = mean(wmat(span::all,d));
      CC = (b1.t()*b1)*ww+CC;
    }
  }
  return CC/K;
}

// get derivatives
// [[Rcpp::export]]
mat df_Gamma_cpp(const NumericVector & y, const mat & evec,
                 const mat & x, const mat & z,  
                 const mat & G, const vec & Gamma_k,
                 const mat & w, const mat & B,const mat & C,
                 int D, int type)
{
  // Here C is n(m-1) by D. These need to be created
  int n = w.n_rows;
  int mn = y.length();
  int m = (mn/n)+1;
  int p = Gamma_k.size();
  int pc = z.n_cols;
  
  mat dnll = zeros<mat>(p,1);
  NumericVector alpha(y.size());
  LogicalVector idx = (y>0);
  alpha[idx]= 1;
  vec alphavec = as<arma::vec>(alpha);
  IntegerVector r_1 = rep_each(seq(0,n-1),m-1);
  IntegerVector r_2 = seq(0,pc-1);
  uvec a_1 = as<uvec>(r_1);
  uvec a_2 = as<uvec>(r_2);
  if(type==1){
    mat X_G = join_rows(x,G);
    mat Xb_1 = zeros<mat>(n,m-1);
    for(int i=0;i<n;i++){
      int s = i*(m-1);
      int e = m-2+i*(m-1);
      mat temp_1 = X_G(span(s,e),span::all)*Gamma_k;
      Xb_1(i,span::all) = trans(temp_1);
    }
    
    mat part_1=zeros<mat>(p,1);
    mat t_3 = zeros<mat>(p,D);
    mat t_0 = zeros<mat>(n*(m-1),p);
    mat t_2 = zeros<mat>(n*(m-1),p);
    mat w_1 = zeros<mat>(n*(m-1),p);
    
    for(int d=0;d<D;d++){
      int s = d*n;
      int e = n-1+d*n;
      mat bb = B(span(s,e),span(0,pc-1));
      mat b = bb(a_1,a_2);
      vec w_d = w(span::all,d);
      vec w_0 = w_d(a_1);
      t_0.each_col() = alphavec%w_0;
      vec t_1 = -trans(sum(X_G%t_0,0));
      w_1.each_col() = w_0;
      vec t_4 = reshape(trans(Xb_1),n*(m-1),1)+sum(z%b,1)+C(span::all,d);
      t_2.each_col() = 1/(1+exp(-t_4));
      t_3(span::all,d) = trans(sum(w_1%X_G%t_2,0))+t_1;
    }
    part_1 = sum(t_3,1);
    dnll = part_1;
  
  }else if(type==2){
    NumericVector yy(y.size());
    for(int i=0;i<(n*(m-1));i++){
      if(y[i]>0){
        yy[i] = log(y[i]);
      }
    }
    vec yvec = as<arma::vec>(yy);
    mat X_G = join_rows(x,G);
    
    vec beta_2 = Gamma_k(span(0,p-2));
    double tau = Gamma_k(p-1);
    mat Xb_2 = zeros<mat>(n,m-1);
    for(int i=0;i<n;i++){
      int s = i*(m-1);
      int e = m-2+i*(m-1);
      mat temp_1 = X_G(span(s,e),span::all)*beta_2;
      Xb_2(i,span::all) = trans(temp_1);
    }
    //part 1 for [beta gamma]
    mat part_1 = zeros<mat>(p-1,1);
    mat t_4 = zeros<mat>(p-1,D);
    mat t_0 = zeros<mat>(n*(m-1),p-1);
    mat t_3 = zeros<mat>(n*(m-1),p-1);
    for(int d=0;d<D;d++){
      int s = d*n;
      int e = n-1+d*n;
      mat bb = B(span(s,e),span(pc,2*pc-1));
      mat b = bb(a_1,a_2);
      vec w_d = w(span::all,d);
      vec w_0 = w_d(a_1);
      t_0.each_col() = alphavec%w_0;
      mat t_1 = X_G%t_0;
      vec t_2 = tau*yvec-reshape(trans(Xb_2),n*(m-1),1)-tau*sum(z%b,1)-tau*C(span::all,d);
      t_3.each_col() = t_2;
      t_4(span::all,d) = -trans(sum(t_3%t_1,0));
    }
    part_1 = sum(t_4,1);
    
    //part 2 for tau
    mat t_4_1 = zeros<mat>(D,1);
    mat part_2 = zeros<mat>(1,1);
    for(int d=0;d<D;d++){
      int s = d*n;
      int e = n-1+d*n;
      mat bb = B(span(s,e),span(pc,2*pc-1));
      mat b = bb(a_1,a_2);
      vec w_d = w(span::all,d);
      vec w_0 = w_d(a_1);
      vec t_0 = alphavec%w_0;
      double t_1 = as_scalar(-(1/tau)*sum(t_0));
      vec t_2 = yvec-sum(z%b,1)-C(span::all,d);
      vec t_3 = (tau*t_2-reshape(trans(Xb_2),n*(m-1),1))%t_2%t_0;
      t_4_1(d,0) = sum(t_3)+t_1;
    }
    part_2(0,0) = as_scalar(sum(t_4_1));
    dnll = join_cols(part_1,part_2);
  
  }else{
    mat X_G = join_rows(x,G);
    mat Xb_1 = zeros<mat>(n,m-1);
    for(int i=0;i<n;i++){
      int s = i*(m-1);
      int e = m-2+i*(m-1);
      mat temp_1 = X_G(span(s,e),span::all)*Gamma_k;
      Xb_1(i,span::all) = trans(temp_1);
    }

    mat part_1=zeros<mat>(p,1);
    mat t_3 = zeros<mat>(p,D);
    mat t_0 = zeros<mat>(n*(m-1),p);
    mat t_2 = zeros<mat>(n*(m-1),p);
    mat w_1 = zeros<mat>(n*(m-1),p);

    for(int d=0;d<D;d++){
      int s = d*n;
      int e = n-1+d*n;
      mat bb = B(span(s,e),span(2*pc,3*pc-1));
      mat b = bb(a_1,a_2);
      vec w_d = w(span::all,d);
      vec w_0 = w_d(a_1);
      t_0.each_col() = alphavec%evec%w_0;
      vec t_1 = -trans(sum(X_G%t_0,0));
      w_1.each_col() = w_0;
      vec t_4 = reshape(trans(Xb_1),n*(m-1),1)+sum(z%b,1)+C(span::all,d);
      t_2.each_col() = alphavec%(1/(1+exp(-t_4)));
      t_3(span::all,d) = trans(sum(w_1%X_G%t_2,0))+t_1;
    }
    part_1 = sum(t_3,1);
    dnll = part_1;
  }
  return dnll;
}

// get nll
// [[Rcpp::export]]
double nll_cpp(const NumericVector & y, const mat & evec,
               const mat & x, const mat & z,  
               const mat & G, const vec & Gamma_k,
               const mat & w, const mat & B,const mat & C,
               int D, int type)
{
  double nll = 0.0;
  // Here C is n(m-1) by D. These need to be created
  int n = w.n_rows;
  int mn = y.length();
  int m = (mn/n)+1;
  int p = Gamma_k.size();
  int pc = z.n_cols;
  static const double pi = 3.14159265;
  
  NumericVector alpha(y.size());
  LogicalVector idx = (y>0);
  alpha[idx]= 1;
  vec alphavec = as<arma::vec>(alpha);
  mat alphamat = trans(reshape(alphavec,m-1,n));
  
  IntegerVector r_1 = rep_each(seq(0,n-1),m-1);
  IntegerVector r_2 = seq(0,pc-1);
  uvec a_1 = as<uvec>(r_1);
  uvec a_2 = as<uvec>(r_2);
  
  mat X_G = join_rows(x,G);
  
  if(type==1){
    
    mat ff = zeros<mat>(n,D);
    mat Xb_1 = zeros<mat>(n,m-1);
    for(int i=0;i<n;i++){
      int s = i*(m-1);
      int e = m-2+i*(m-1);
      mat temp_1 = X_G(span(s,e),span::all)*Gamma_k;
      Xb_1(i,span::all) = trans(temp_1);
    }
    for(int d=0;d<D;d++){
      int s = d*n;
      int e = n-1+d*n;
      mat bb = B(span(s,e),span(0,pc-1));
      mat b = bb(a_1,a_2);
      
      mat temp_1 = sum(z%b,1)+C(span::all,d);
      mat temp_2 = Xb_1 + trans(reshape(temp_1,m-1,n));
      mat f_0 = log(1+exp(temp_2));
      uvec ii = find_nonfinite(f_0);
      f_0(ii) = temp_2(ii);
      vec f_1= sum(f_0,1);
      vec f_2= -sum(alphamat%temp_2,1);
      ff(span::all,d) = f_1+f_2;
    }
    nll = accu(ff%w);
    
  }else if(type==2){
    mat g = zeros<mat>(n,D);
    NumericVector yy(y.size());
    for(int i=0;i<(n*(m-1));i++){
      if(y[i]>0){
        yy[i] = log(y[i]);
      }
    }
    vec yvec = as<arma::vec>(yy);
    mat ymat = trans(reshape(yvec, m-1,n));
    vec beta_2 = Gamma_k(span(0,p-2));
    double tau = Gamma_k(p-1);
    mat Xb_2 = zeros<mat>(n,m-1);
    for(int i=0;i<n;i++){
      int s = i*(m-1);
      int e = m-2+i*(m-1);
      mat temp_1 = X_G(span(s,e),span::all)*beta_2;
      Xb_2(i,span::all) = trans(temp_1);
    }
    
    vec g_1 = -sum(alphamat%(ymat+0.5*log(2*pi)),1);
    vec g_2 = sum(log(tau)*alphamat,1);
    for(int d=0;d<D;d++){
      int s = d*n;
      int e = n-1+d*n;
      mat bb = B(span(s,e),span(pc,2*pc-1));
      mat b = bb(a_1,a_2);
      vec temp_1 = tau*(sum(z%b,1)+C(span::all,d));
      mat temp_2 = Xb_2+trans(reshape(temp_1,m-1,n));
      mat mu_1 = (tau*ymat-temp_2)%(tau*ymat-temp_2);
      vec g_3 = -0.5*sum(alphamat%mu_1,1);
      g(span::all,d) = (-g_1-g_2-g_3);
    }
    nll = accu(g%w);//sum(g*w)
  }else{
    mat epsmat = trans(reshape(evec,m-1,n));
    mat ff = zeros<mat>(n,D);
    mat Xb_1 = zeros<mat>(n,m-1);
    for(int i=0;i<n;i++){
      int s = i*(m-1);
      int e = m-2+i*(m-1);
      mat temp_1 = X_G(span(s,e),span::all)*Gamma_k;
      Xb_1(i,span::all) = trans(temp_1);
    }
    for(int d=0;d<D;d++){
      int s = d*n;
      int e = n-1+d*n;
      mat bb = B(span(s,e),span(2*pc,3*pc-1));
      mat b = bb(a_1,a_2);
      
      mat temp_1 = sum(z%b,1)+C(span::all,d);
      mat temp_2 = Xb_1 + trans(reshape(temp_1,m-1,n));
      mat f_0 = log(1+exp(temp_2));
      uvec ii = find_nonfinite(f_0);
      f_0(ii) = temp_2(ii);
      vec f_1= sum(alphamat%f_0,1);
      vec f_2= -sum(alphamat%epsmat%temp_2,1);
      ff(span::all,d) = f_1+f_2;
    }
    nll = accu(ff%w);
  }
  return nll;
}

// Backtracking line search
// [[Rcpp::export]]
vec linesearch_cpp(const NumericVector & y, const mat & evec,
                   const mat & x, const mat & z,  
                   const mat & G,
                   const vec & Gamma_old, const vec & df,
                   const mat & w, const mat & B, const mat & C,
                   int D, int type)
{
  double lr_new = 1.0;
  double t_1 = 1+norm(df);
  vec df_Gamma = df;
  int p = Gamma_old.size();
  
  double fb = nll_cpp(y,evec,x,z,G,Gamma_old,w,B,C,D,type);
  vec Gamma_new = Gamma_old-lr_new*(df_Gamma/t_1);
  //---------------------------
  if(type==2 && Gamma_new(p-1)<0.25){
    
    Gamma_new(p-1) = 0.25;
    df_Gamma(p-1) = (Gamma_old(p-1)-Gamma_new(p-1))/lr_new;
    t_1 = 1+norm(df_Gamma);
  }
  //---------------------------------------
  double fb_new = nll_cpp(y,evec,x,z,G,Gamma_new,w,B,C,D,type);
  double temp = 1.0;//accu(df_Gamma%df_Gamma);
  double ee = fb-lr_new*0.5*temp;
  while(fb_new>ee){
    lr_new = 0.2*lr_new;
    Gamma_new = Gamma_old-lr_new*(df_Gamma/t_1);
    if(type==2 && Gamma_new(p-1)<0.25){
      Gamma_new(p-1) = 0.25;
      df_Gamma(p-1) = (Gamma_old(p-1)-Gamma_new(p-1))/lr_new;
      t_1 = 1+norm(df_Gamma);
      temp = 1;//accu(df_Gamma%df_Gamma);
    }
    fb_new = nll_cpp(y,evec,x,z,G,Gamma_new,w,B,C,D,type);
    ee = fb-lr_new*0.5*temp;
  }
  // double lr = lr.new;
  return Gamma_new;
}

// get gradient descent estimates
// [[Rcpp::export]]
vec get_gdest_cpp(const NumericVector & y, const mat & evec, 
                  const mat & x, const mat & z,  
                  const mat & G, const vec & Gamma_k,
                  const mat & w, const mat & B,const mat & C,
                  int D, int type,double tol)
{
  
  int p = Gamma_k.size();
  
  vec gdest(p);
  if(type<=3){
    vec df_Gamma(p);
    vec Gamma_old = Gamma_k;
    double err = 1.0;
    
    while(err > tol){
      
      df_Gamma = df_Gamma_cpp(y,evec,x,z,G,Gamma_old,w,B,C,D,type);
      
      //do line search & update
      vec Gamma_new = linesearch_cpp(y,evec,x,z,G,Gamma_old,df_Gamma,w,B,C,D,type);
      
      //error update
      err = norm(Gamma_new-Gamma_old,2);
      Gamma_old = Gamma_new;
      gdest = Gamma_new;
    }
  }else{
    int n = w.n_rows;
    int mn = y.length();
    int m = (mn/n)+1;
    int pc = z.n_cols;
    mat Amat = zeros<mat>(p,p);
    mat cmat = zeros<mat>(p,D);
    mat t_0 = zeros<mat>(n*(m-1),p);
    mat t_1 = zeros<mat>(m-1,p);
    
    NumericVector alpha(y.size());
    LogicalVector idx = (y>0);
    alpha[idx]= 1;
    vec alphavec = as<arma::vec>(alpha);
    mat X_G = join_rows(x,G);
    vec Xb_1 = X_G*Gamma_k;
    
    IntegerVector r_1 = rep_each(seq(0,n-1),m-1);
    IntegerVector r_2 = seq(0,pc-1);
    uvec a_1 = as<uvec>(r_1);
    uvec a_2 = as<uvec>(r_2);
    
    for(int d=0;d<D;d++){
      
      int s = d*n;
      int e = n-1+d*n;
      mat bb = B(span(s,e),span(2*pc,3*pc-1));
      mat b = bb(a_1,a_2);
      
      vec temp_1 = Xb_1+sum(z%b,1)+C(span::all,d);
      vec f_0 = 1/(1+exp(-temp_1));
      vec w_d = w(span::all,d);
      vec w_0 = w_d(a_1);
      t_0.each_col() = alphavec%(evec-f_0)%w_0;
      mat f_1= t_0%X_G;//sum((t_0%X_G),0);
      cmat(span::all,d) = trans(sum(f_1));
      for(int i=0;i<n;i++){
        int s_1 = i*(m-1);
        int e_1 = m-2+i*(m-1);
        t_1.each_col() = alphavec(span(s_1,e_1));
        mat temp_2 = X_G(span(s_1,e_1),span::all)%t_1;
        Amat = (trans(temp_2)*temp_2)*w(i,d)+Amat;
      }
    }
    gdest = Gamma_k+inv(Amat)*sum(cmat,1);
  }
  return gdest;
}

// get w matrix
// [[Rcpp::export]]
mat get_wmat_cpp(const NumericVector & y, const mat & evec, 
                 const mat & x, const mat & z,  
                 const mat & G, const vec & Gamma_1,const vec & Gamma_2,
                 const vec & Gamma_3,const mat & B,
                 const mat & C_1,const mat & C_2,const mat & C_3,
                 int n, int D)
{
  
  mat w(n,D);
  int mn = y.length();
  int m = (mn/n)+1;
  int p_2 = Gamma_2.size();
  int pc = z.n_cols;
  static const double pi = 3.14159265;
  
  NumericVector alpha(y.size());
  LogicalVector idx = (y>0);
  alpha[idx]= 1;
  vec alphavec = as<arma::vec>(alpha);
  mat alphamat = trans(reshape(alphavec,m-1,n));
  mat emat = trans(reshape(evec,m-1,n));
  NumericVector yy(y.size());
  for(int i=0;i<(n*(m-1));i++){
    if(y[i]>0){
      yy[i] = log(y[i]);
    }
  }
  vec yvec = as<arma::vec>(yy);
  mat ymat = trans(reshape(yvec, m-1,n));
  mat X_G = join_rows(x,G);
  
  IntegerVector r_1 = rep_each(seq(0,n-1),m-1);
  IntegerVector r_2 = seq(0,3*pc-1);
  uvec a_1 = as<uvec>(r_1);
  uvec a_2 = as<uvec>(r_2);
  
  vec beta_2 = Gamma_2(span(0,p_2-2));
  double tau = Gamma_2(p_2-1);
  mat Xb_1 = zeros<mat>(n,m-1);
  mat Xb_2 = zeros<mat>(n,m-1);
  mat Xb_3 = zeros<mat>(n,m-1);
  for(int i=0;i<n;i++){
    int s = i*(m-1);
    int e = m-2+i*(m-1);
    mat temp = X_G(span(s,e),span::all)*Gamma_1;
    Xb_1(i,span::all) = trans(temp);
    temp = X_G(span(s,e),span::all)*beta_2;
    Xb_2(i,span::all) = trans(temp);
    temp = X_G(span(s,e),span::all)*Gamma_3;
    Xb_3(i,span::all) = trans(temp);
  }
  
  vec g_1 = -sum(alphamat%(ymat+0.5*log(2*pi)),1);
  vec g_2 = sum(log(tau)*alphamat,1);
  
  for(int d=0;d<D;d++){
    
    int s = d*n;
    int e = n-1+d*n;
    mat bb = B(span(s,e),span(0,3*pc-1));
    mat b = bb(a_1,a_2);
    
    // login Indicator
    vec temp_1 = sum(z%b(span::all,span(0,pc-1)),1)+C_1(span::all,d);
    mat temp_2 = Xb_1 + trans(reshape(temp_1,m-1,n));
    mat f_0 = -log(1+exp(temp_2));
    uvec ii = find_nonfinite(f_0);
    f_0(ii) = -temp_2(ii);
    vec f_1= sum(f_0,1);
    vec f_2= sum(alphamat%temp_2,1);
    vec ff = f_1+f_2;
    
    //Play Time
    temp_1 = tau*(sum(z%b(span::all,span(pc,2*pc-1)),1)+C_2(span::all,d));
    temp_2 = Xb_2+trans(reshape(temp_1,m-1,n));
    mat mu_1 = (tau*ymat-temp_2)%(tau*ymat-temp_2);
    vec g_3 = -0.5*sum(alphamat%mu_1,1);
    vec gg = g_1+g_2+g_3;
    
    // Purchase Indicator
    temp_1 = sum(z%b(span::all,span(2*pc,3*pc-1)),1)+C_3(span::all,d);
    temp_2 = Xb_3 + trans(reshape(temp_1,m-1,n));
    mat e_0 = -log(1+exp(temp_2));
    uvec jj = find_nonfinite(e_0);
    e_0(jj) = -temp_2(jj);
    vec e_1= sum(alphamat%e_0,1);
    vec e_2= sum(alphamat%emat%temp_2,1);
    vec ee = e_1+e_2;
    
    w(span::all,d) = exp(ff+gg+ee);
  }
  // vec rr = sum(w,1);
  // mat rs(n,D);
  // rs.each_col() = 1/rr;
  // mat wmat = w%rs;
  return w;//wmat;
  
}
