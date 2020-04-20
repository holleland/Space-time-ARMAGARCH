
// DETTE MÅ FIKSES PÅ: NOEN FEIL MED DIM
#include <TMB.hpp>                                // Links in the TMB libraries
#include <math.h>





template<class Type>
  Type objective_function<Type>::operator() ()
{
  // DATA 
  DATA_ARRAY(x);
  DATA_SCALAR(mut);
  
  
  // PARAMETERS
  PARAMETER(mu);
  PARAMETER(phi1);
  PARAMETER(phi2);
  PARAMETER(theta1);
  PARAMETER(theta2);
  PARAMETER(sigma);
  
  vector<int> x_dim = x.dim; // dimension of x
  array<Type> u(x_dim(0),x_dim(1));
  array<Type> ell(x_dim(0),x_dim(1)-1);
  
  // subtracting the mean parameter
  for(int k=0;k<x_dim(0);k++){
    for(int t=0; t<x_dim(1); t++){
      x(k,t) -= mu;
    }
    u(k,0)=0.0;
  }
  
  
  Type f;
 
  f=0;
  // AR PART
  for(int t=0;t<x_dim(1)-1;t++){
    for(int k=0;k<x_dim(0);k++){
      // 0<k<m1-1
      if( (k > 0) && (k < x_dim(0)-1)){
        u(k,t+1) =  x(k,t+1)  - phi1*x(k,t)-phi2*(x(k-1,t)+x(k+1,t))-theta1*u(k,t)-theta2*(u(k-1,t)+u(k+1,t));
        // k=0 
      } else if(k ==0) {
        u(k,t+1) = x(k,t+1) -phi1*x(k,t)-phi2*(x(x_dim(0)-1,t)+x(k+1,t))- theta1*u(k,t)-theta2*(u(x_dim(0)-1,t)+u(k+1,t));
        // k=m1-1 
      } else if(k == x_dim(0)-1){
        u(k,t+1) = x(k,t+1)  -phi1*x(k,t)-phi2*(x(k-1,t)+x(0,t))- theta1*u(k,t)-theta2*(u(k-1,t)+u(0,t));
      }
      //if(t>10){
        ell(k,t) = -dnorm(u(k,t+1),mut,sigma,true);
        f -= dnorm(u(k,t+1),mut,sigma,true);
        
        //}
    }
  }
  ADREPORT(ell);
  return f;
  }
