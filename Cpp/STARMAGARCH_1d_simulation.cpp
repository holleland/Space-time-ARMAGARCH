#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix simulate_1d(Rcpp::NumericVector par, Rcpp::IntegerVector dim, int burnin_time, int burnin_space){
  int n = dim(1)+burnin_time;
  int m = dim(0)+2*burnin_space;
  //arma::cube sigma = arma::ones(m1,m2,n);
  NumericMatrix x(m,n);
  NumericMatrix y(m,n);  
  NumericMatrix sigma(m,n);
  NumericMatrix returx(dim(0),dim(1));
  
  for(int t=0;t<n;t++){
    for(int k = 0; k < m; k++){
      x(k,t)= 0.0;//sqrt(par(5)/(1-3*par(6)-3*par(7)));    
      y(k,t)=0.0;
      sigma(k,t)=0.0;//par(5)/(1-3*par(6)-3*par(7));
    }
  }

  for(int t = 1;t<n;t++) {
    for(int k=1;k<m-1;k++){
      sigma(k,t) = par(5)+ par(6)*(x(k-1,t-1)*x(k-1,t-1)+x(k,t-1)*x(k,t-1)+x(k+1,t-1)*x(k+1,t-1))+par(7)*(sigma(k-1,t-1)+sigma(k,t-1)+sigma(k+1,t-1));
        x(k,t) = sqrt(sigma(k,t))*rnorm(1)[0];
    }
  }
  for(int t = 1;t<n;t++) {
    for(int k=1;k<m-1;k++){
      y(k,t) =  par(1)*y(k,t-1)+par(2)*(y(k-1,t-1)+y(k+1,t-1))+ par(3)*x(k,t-1)+par(4)*(x(k-1,t-1)+x(k+1,t-1))+x(k,t);
    }
  }
  
  for(int t = 0;t<dim(1);t++) {
    for(int k=0;k<dim(0);k++){
           returx(k,t) = y(burnin_space+k, burnin_time+t)+par(0);
      }
  }
  return returx;
}
