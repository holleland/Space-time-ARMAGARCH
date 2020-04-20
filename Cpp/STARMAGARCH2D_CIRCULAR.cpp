
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
  DATA_ARRAY(ustart);
  //DATA_SCALAR(theta1);
  //DATA_SCALAR(theta2);
 // DATA_MATRIX(w); //neighbour matrix size m x m
  //DATA_SCALAR(h0);

  // PARAMETERS
  PARAMETER(mu);
  PARAMETER(phi1);
  PARAMETER(phi2);
  //PARAMETER(phi3);
  PARAMETER(theta1);
  PARAMETER(theta2);
  //PARAMETER(theta3);
  PARAMETER(omega);
  PARAMETER(alpha1);
  PARAMETER(alpha2);
  //PARAMETER(alpha3);
  PARAMETER(beta1);
  PARAMETER(beta2);
  //PARAMETER(beta3);
  //PARAMETER(gam1);
  
  Type phi3 = phi2;
  Type theta3 =theta2;
  Type alpha3 =alpha2;
  Type beta3 =beta2;
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
  array<Type> u(x_dim(0),x_dim(1),x_dim(2));
  array<Type> sigma(x_dim(0),x_dim(1),x_dim(2));
  array<Type> u2(x_dim(0),x_dim(1),x_dim(2));
  
  
  // subtracting the mean parameter and setting initial values
  for(int k=0;k<x_dim(0);k++){
    for(int l=0;l<x_dim(1);l++){
    for(int t=0; t<x_dim(2); t++){
        x(k,l,t) -= mu;
      }
      u(k,l,0)=ustart(k,l);
      u2(k,l,0)=u(k,l,0)*u(k,l,0);
      sigma(k,l,0) =u2(k,l,0);
    }
    }
  
  
  //sigma.setOnes();
  //sigma = h0*sigma;

  // AR PART
  for(int t=0;t < x_dim(2)-1;t++){
    for(int k=0;k < x_dim(0);k++){
      for(int l=0;l < x_dim(1);l++){
      // 0<k<m1-1 0<l<m2-1
        if( (k > 0) && (k < x_dim(0)-1)&& (l > 0) && (l < x_dim(1)-1)){
          u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(k-1,l,t)+x(k+1,l,t))-phi3*(x(k,l-1,t)+x(k,l+1,t))-theta1*u(k,l,t)-theta2*(u(k-1,l,t)+u(k+1,l,t))-theta3*(u(k,l-1,t)+u(k,l+1,t));
        }else if( (k == 0) && (l > 0) && (l < x_dim(1)-1)){
          u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(x_dim(0)-1,l,t)+x(k+1,l,t))-phi3*(x(k,l-1,t)+x(k,l+1,t))-theta1*u(k,l,t)-theta2*(u(x_dim(0)-1,l,t)+u(k+1,l,t))-theta3*(u(k,l-1,t)+u(k,l+1,t));
        }else if( (k == x_dim(0)-1)&& (l > 0) && (l < x_dim(1)-1)){
            u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(k-1,l,t)+x(0,l,t))-phi3*(x(k,l-1,t)+x(k,l+1,t))-theta1*u(k,l,t)-theta2*(u(k-1,l,t)+u(0,l,t))-theta3*(u(k,l-1,t)+u(k,l+1,t));
        }else if( (k > 0) && (k < x_dim(0)-1)&& (l ==0 )){
          u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(k-1,l,t)+x(k+1,l,t))-phi3*(x(k,x_dim(1)-1,t)+x(k,l+1,t))-theta1*u(k,l,t)-theta2*(u(k-1,l,t)+u(k+1,l,t))-theta3*(u(k,x_dim(1)-1,t)+u(k,l+1,t));
        }else if( (k > 0) && (k < x_dim(0)-1)&& (l == x_dim(1)-1)){
          u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(k-1,l,t)+x(k+1,l,t))-phi3*(x(k,l-1,t)+x(k,0,t))-theta1*u(k,l,t)-theta2*(u(k-1,l,t)+u(k+1,l,t))-theta3*(u(k,l-1,t)+u(k,0,t));
        }else if( (k == 0) && (l == 0)){
          u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(x_dim(0)-1,l,t)+x(k+1,l,t))-phi3*(x(k,x_dim(1)-1,t)+x(k,l+1,t))-theta1*u(k,l,t)-theta2*(u(x_dim(0)-1,l,t)+u(k+1,l,t))-theta3*(u(k,x_dim(1)-1,t)+u(k,l+1,t));
        }else if( (k == 0) && (l == x_dim(1)-1)){
          u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(x_dim(0)-1,l,t)+x(k+1,l,t))-phi3*(x(k,l-1,t)+x(k,0,t))-theta1*u(k,l,t)-theta2*(u(x_dim(0)-1,l,t)+u(k+1,l,t))-theta3*(u(k,l-1,t)+u(k,0,t));
        }else if( (k == x_dim(0)-1)&& (l == 0) ){
          u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(k-1,l,t)+x(0,l,t))-phi3*(x(k,x_dim(1)-1,t)+x(k,l+1,t))-theta1*u(k,l,t)-theta2*(u(k-1,l,t)+u(0,l,t))-theta3*(u(k,x_dim(1)-1,t)+u(k,l+1,t));
        }else if( (k == x_dim(0)-1)&& (l == x_dim(1)-1)){
          u(k,l,t+1) =  x(k,l,t+1) - phi1*x(k,l,t)-phi2*(x(k-1,l,t)+x(0,l,t))-phi3*(x(k,l-1,t)+x(k,0,t))-theta1*u(k,l,t)-theta2*(u(k-1,l,t)+u(0,l,t))-theta3*(u(k,l-1,t)+u(k,0,t));
        }
        u2(k,l,t+1) = u(k,l,t+1)*u(k,l,t+1);
        
      }
    }
  }
    
  // GARCH part

  // Declare the "objective function" (neg. log. likelihood)
  Type f;
  f=0;
  for (int t = 0; t < x_dim(2)-1; t++) { 
         for(int k=0;k<x_dim(0);k++){
           for(int l=0;l<x_dim(1);l++){
             if( (k > 0) && (k < x_dim(0)-1)&& (l > 0) && (l < x_dim(1)-1)){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(k-1,l,t)+u2(k+1,l,t))+
                 alpha3*(u2(k,l-1,t)+u2(k,l+1,t))+
                 beta1*sigma(k,l,t)+beta2*(sigma(k-1,l,t)+sigma(k+1,l,t))+
                 beta3*(sigma(k,l-1,t)+sigma(k,l+1,t));
             }else if( (k == 0) && (l > 0) && (l < x_dim(1)-1)){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(x_dim(0)-1,l,t)+u2(k+1,l,t))+alpha3*(u2(k,l-1,t)+u2(k,l+1,t))+beta1*sigma(k,l,t)+beta2*(sigma(x_dim(0)-1,l,t)+sigma(k+1,l,t))+beta3*(sigma(k,l-1,t)+sigma(k,l+1,t));
             }else if( (k == x_dim(0)-1)&& (l > 0) && (l < x_dim(1)-1)){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(k-1,l,t)+u2(0,l,t))+alpha3*(u2(k,l-1,t)+u2(k,l+1,t))+beta1*sigma(k,l,t)+beta2*(sigma(k-1,l,t)+sigma(0,l,t))+beta3*(sigma(k,l-1,t)+sigma(k,l+1,t));
             }else if( (k > 0) && (k < x_dim(0)-1)&& (l ==0 )){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(k-1,l,t)+u2(k+1,l,t))+alpha3*(u2(k,x_dim(1)-1,t)+u2(k,l+1,t))+beta1*sigma(k,l,t)+beta2*(sigma(k-1,l,t)+sigma(k+1,l,t))+beta3*(sigma(k,x_dim(1)-1,t)+sigma(k,l+1,t));
             }else if( (k > 0) && (k < x_dim(0)-1)&& (l == x_dim(1)-1)){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(k-1,l,t)+u2(k+1,l,t))+alpha3*(u2(k,l-1,t)+u2(k,0,t))+beta1*sigma(k,l,t)+beta2*(sigma(k-1,l,t)+sigma(k+1,l,t))+beta3*(sigma(k,l-1,t)+sigma(k,0,t));
             }else if( (k == 0) && (l == 0)){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(x_dim(0)-1,l,t)+u2(k+1,l,t))+alpha3*(u2(k,x_dim(1)-1,t)+u2(k,l+1,t))+beta1*sigma(k,l,t)+beta2*(sigma(x_dim(0)-1,l,t)+sigma(k+1,l,t))+beta3*(sigma(k,x_dim(1)-1,t)+sigma(k,l+1,t));
             }else if( (k == 0) && (l == x_dim(1)-1)){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(x_dim(0)-1,l,t)+u2(k+1,l,t))+alpha3*(u2(k,l-1,t)+u2(k,0,t))+beta1*sigma(k,l,t)+beta2*(sigma(x_dim(0)-1,l,t)+sigma(k+1,l,t))+beta3*(sigma(k,l-1,t)+sigma(k,0,t));
             }else if( (k == x_dim(0)-1)&& (l == 0) ){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(k-1,l,t)+u2(0,l,t))+alpha3*(u2(k,x_dim(1)-1,t)+u2(k,l+1,t))+beta1*sigma(k,l,t)+beta2*(sigma(k-1,l,t)+sigma(0,l,t))+beta3*(sigma(k,x_dim(1)-1,t)+sigma(k,l+1,t));
             }else if( (k == x_dim(0)-1)&& (l == x_dim(1)-1)){
               sigma(k,l,t+1) = omega + alpha1*u2(k,l,t)+alpha2*(u2(k-1,l,t)+u2(0,l,t))+alpha3*(u2(k,l-1,t)+u2(k,0,t))+beta1*sigma(k,l,t)+beta2*(sigma(k-1,l,t)+sigma(0,l,t))+beta3*(sigma(k,l-1,t)+sigma(k,0,t));
             }
               
// if(t>10) {
            f -= dnorm(u(k,l,t+1),mut,sqrt(sigma(k,l,t+1)),true);
             //}
          }
         }
  }
    return f;
  }
