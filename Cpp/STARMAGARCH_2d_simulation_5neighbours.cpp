#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
arma::cube simulate_2d(Rcpp::NumericVector par, Rcpp::IntegerVector dim, int burnin_time, int burnin_space){
  int n = dim(2)+burnin_time;
  int m1 = dim(0)+2*burnin_space;
  int m2 = dim(1)+2*burnin_space;
  //arma::cube sigma = arma::ones(m1,m2,n);
  arma::cube x(m1,m2,n);
  arma::cube y(m1,m2,n);
  for(int t=0;t<n;t++){
    for(int k = 0; k < m1; k++){
      for(int l = 0; l < m2; l++){
      x(k,l,t)= 0.0 ;//par(7)/(1-3*par(8)-3*par());    
      y(k,l,t)= 0.0;
      }
    }
  }
  
  arma::cube sigma =x;
  arma::cube returx(dim(0),dim(1),dim(2));
  
  
  for(int t = 1;t<n;t++) {
    for(int k=1;k<m1-1;k++){
      for(int l=1;l<m2-1;l++){
      sigma(k,l,t) = par(7)+ par(8)*x(k,l,t-1)*x(k,l,t-1)+par(9)*(x(k-1,l,t-1)*x(k-1,l,t-1)+x(k+1,l,t-1)*x(k+1,l,t-1))+par(10)*(x(k,l-1,t-1)*x(k,l-1,t-1)+x(k,l+1,t-1)*x(k,l+1,t-1))+
         par(11)*sigma(k,l,t-1)+par(12)*(sigma(k-1,l,t-1)+sigma(k+1,l,t-1))+par(13)*(sigma(k,l-1,t-1)+sigma(k,l+1,t-1));
        x(k,l,t) = sqrt(sigma(k,l,t))*rnorm(1)[0];
      }
    }
  }
  for(int t = 1;t<n;t++) {
    for(int k=1;k<m1-1;k++){
      for(int l=1;l<m2-1;l++){
      y(k,l,t) =  par(1)*y(k,l,t-1)+par(2)*(y(k-1,l,t-1)+y(k+1,l,t-1))+par(3)*(y(k,l-1,t-1)+y(k,l+1,t-1))+
        par(4)*x(k,l,t-1)+par(5)*(x(k-1,l,t-1)+x(k+1,l,t-1))+par(6)*(x(k,l-1,t-1)+x(k,l+1,t-1))+x(k,l,t);
      }
    }
  }
  
  for(int t = 0;t<dim(2);t++) {
    for(int k=0;k<dim(0);k++){
      for(int l=0;l<dim(1);l++){
           returx(k,l,t) = y(burnin_space+k,burnin_space+l, burnin_time+t)+par(0);
      }
    }
  }
  return returx;
}