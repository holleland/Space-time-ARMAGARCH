
// DETTE MÅ FIKSES PÅ: NOEN FEIL MED DIM
#include <TMB.hpp>                                // Links in the TMB libraries
#include <math.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA ARRAY/MATRIX
  //DATA_MATRIX(x);
  DATA_ARRAY(x);
  DATA_SCALAR(mut);
  DATA_VECTOR(ustart);
  //DATA_SCALAR(theta1);
  //DATA_SCALAR(theta2);
 // DATA_MATRIX(w); //neighbour matrix size m x m
  //DATA_SCALAR(h0);

  // PARAMETERS
  PARAMETER(mu);
  PARAMETER(phi1);
  PARAMETER(phi2);
  PARAMETER(theta1);
  PARAMETER(theta2);
  PARAMETER(omega);
  PARAMETER(alpha);
  PARAMETER(beta);
  //PARAMETER(gam1);

  //PARAMETER(phi4);
  

  //making sure parameters are positive
  //omega = exp(omega);
  //alpha = exp(alpha);
  //beta  = exp(beta);
  
  //Type h0 = omega;
  //Type x0 = omega;
  //
  vector<int> x_dim = x.dim; // dimension of x
  //REPORT(x_dim);
  array<Type> u(x_dim(0),x_dim(1));
  array<Type> ell(x_dim(0),x_dim(1)-1);
  array<Type> sigma(x_dim(0),x_dim(1));
  array<Type> x2(x_dim(0),x_dim(1));
  
  
  // subtracting the mean parameter
  for(int k=0;k<x_dim(0);k++){
    for(int t=0; t<x_dim(1); t++){
        x(k,t) -= mu;
      }
      u(k,0)=ustart(k);
    }
  
  
  //sigma.setOnes();
  //sigma = h0*sigma;

  // AR PART
  for(int t=0;t<x_dim(1)-1;t++){
    for(int k=0;k<x_dim(0);k++){
      // 0<k<m1-1
        if( (k > 0) && (k < x_dim(0)-1)){
          u(k,t+1) =  x(k,t+1) - phi1*x(k,t)-phi2*(x(k-1,t)+x(k+1,t))-theta1*u(k,t)-theta2*(u(k-1,t)+u(k+1,t));
      // k=0 
        } else if(k ==0) {
          u(k,t+1) = x(k,t+1) -phi1*x(k,t)-phi2*(x(x_dim(0)-1,t)+x(k+1,t))- theta1*u(k,t)-theta2*(u(x_dim(0)-1,t)+u(k+1,t));
      // k=m1-1 
        } else if(k == x_dim(0)-1){
          u(k,t+1) = x(k,t+1) -phi1*x(k,t)-phi2*(x(k-1,t)+x(0,t))- theta1*u(k,t)-theta2*(u(k-1,t)+u(0,t));
      }
    }
  }
    
  // GARCH part

  for(int k=0;k<x_dim(0);k++){
      for(int t = 0;t<x_dim(1);t++) { //x_dim(2)-1
         x2(k,t) = u(k,t)*u(k,t);
        }
        sigma(k,0) =x2(k,0);// omega;///(1-3*(alpha+beta)); //x2(k,0);
      }

  
  // Declare the "objective function" (neg. log. likelihood)
  Type f;
  f=0;
  for (int t = 0; t < x_dim(1)-1; t++) { 
         for(int k=0;k<x_dim(0);k++){
           if( (k > 0) && (k < x_dim(0)-1) ){
           sigma(k,t+1) = omega + alpha*(x2(k-1,t)+x2(k,t)+x2(k+1,t))+beta*(sigma(k-1,t)+sigma(k,t)+sigma(k+1,t));
           } else if(k ==0){
             sigma(k,t+1) = omega + alpha*(x2(x_dim(0)-1,t)+x2(k,t)+x2(k+1,t))+beta*(sigma(x_dim(0)-1,t)+sigma(k,t)+sigma(k+1,t));
           } else if(k == x_dim(0)-1){
             sigma(k,t+1) = omega + alpha*(x2(k-1,t)+x2(k,t)+x2(0,t))+ beta*(sigma(k-1,t)+sigma(k,t)+sigma(0,t)+sigma(k-1,t));
             }
           //if(t>10) {
            ell(k,t)=-dnorm(u(k,t+1),mut,sqrt(sigma(k,t+1)),true);
            f -= dnorm(u(k,t+1),mut,sqrt(sigma(k,t+1)),true);
            //}
          }
      
  }
    ADREPORT(ell);
    return f;
  }
