simulate.arma.garch<-function(mu=5,phi1=.35,phi2=-.09,th1=.12,th2=-.05,
                              omg=.5,a=.13,b=.1,burnin = 1000, pdim = c(30,200),ncspace = 50){
  x <- u<- h<- array(0, dim=pdim+c(2*ncspace,burnin))
  z <- array(rnorm(prod(dim(x))), dim=dim(x))
  for(t in 1:(dim(x)[2]-1)){
    for(k in 1:dim(x)[1]){
      h[k,t+1]<- sqrt(omg + a * (u[(k-2)%%pdim[1]+1,t]^2+u[k,t]^2+u[k%%pdim[1]+1,t]^2)+
                        b * (h[(k-2)%%pdim[1]+1,t]^2+h[k,t]^2+h[k%%pdim[1]+1,t]^2))
      u[k,t+1]<- h[k,t+1]*z[k,t+1]
    }
  }
  
  for(t in 1:(dim(x)[2]-1)){
    for(k in 1:dim(x)[1]){
      x[k,t+1]<- phi1*x[k,t]+phi2*(x[(k-2)%%pdim[1]+1,t]+x[k%%pdim[1]+1,t])+
        th1*u[k,t]+th2*(u[(k-2)%%pdim[1]+1,t]+u[k%%pdim[1]+1,t])+u[k,t+1]
    }
  }
  
  #removing burn-in
  x<-x[ncspace+1:pdim[1],-(1:burnin)];#h<-h[ncspace+1:pdim[1],-(1:burnin)];u<-u[ncspace+1:pdim[1],-(1:burnin)];
  #z<-z[ncspace+1:pdim[1],-(1:burnin)]
  
  return(x+mu)
}

arma.garch.bootstrap <- function(lambda, burnin, pdim,ncspace, B=200){
  boot.res <- matrix(NA_real_, ncol = length(lambda), nrow=B)
  for(i in 1:B){
    #x <- simulate.arma.garch(mu=lambda[1],phi1=lambda[2],phi2=lambda[3],th1=lambda[4],th2=lambda[5],
                             #omg=lambda[6],a=lambda[7],b=lambda[8],
                             #burnin=burnin,pdim=pdim,ncspace=ncspace) 
    x <- simulate_1d(par = lambda,dim = pdim, burnin_time = burnin, burnin_space = ncspace)
    f = MakeADFun(data=list(x=x, mut=0.0, y = x*0, ustart = rep(0,dim(x)[1]),
                            gam1=0),
                  parameters=list(mu = mean(as.vector(x)), phi1 = .25, phi2=-0.1, theta1=0.1, theta2=0.1,
                                  omega=.33, alpha=.02, beta=.02))
    fit <- nlminb(f$par,f$fn,f$gr, f$he, lower=c(-200,
                                                 -5,-5,
                                                 -5,-5,
                                                 0.0001,0.0001,0.0001),
                  upper=c(500,
                          5,5,
                          5,5,
                          1000,5,5))
    boot.res[i,] <- fit$par
  }
  return(colMeans(boot.res))
}

run.arma.garch.bootstrap <- function(par = c(5.032,.354,-.153,.122,-.051,.513,.142,.178),
                                     burnin = 1000, pdim = c(30,200),ncspace = 50,
                                     B=200){
  
  par.names <- c("mu","phi1","phi2","th1","th2","omg","a","b")
  mu <- par[1]; phi1 <- par[2];phi2 <- par[3]; th1 <- par[4]; th2 <- par[5]
  omg <- par[6];a <- par[7]; b <- par[8]
  
  # Estimation
  
  x <- simulate_1d(par=c(mu,phi1,phi2,th1,th2,omg,a,b), dim=pdim, burnin_time = burnin,
                   burnin_space = ncspace)
    #simulate.arma.garch(mu,phi1,phi2,th1,th2,omg,a,b,burnin,pdim,ncspace) 
  
  #compile("Cpp/STARMAGARCH1D_CIRCULAR.cpp")         # Compile the C++ file
  f = MakeADFun(data=list(x=x, mut=0.0, ustart = rep(0,pdim[1])),
                parameters=list(mu = mean(as.vector(x)), phi1 = .25, phi2=-0.1, theta1=0.1, theta2=0.1,
                                omega=.33, alpha=.02, beta=.02))
  fit <- nlminb(f$par,f$fn,f$gr, f$he, 
                lower=c(-200,-5,-5,
                        -5,-5,
                        0.0001,0.0001,0.0001),
                upper=c(500,
                        5,5,
                        5,5,
                        1000,5,5))
  #fit
  Infomat<-f$he(fit$par)
  pbb<-arma.garch.bootstrap(fit$par, burnin=burnin, 
                            pdim=pdim,ncspace=ncspace,B=B)
  #pbb2 <- arma.garch.bootstrap(2*fit$par-pbb, burnin = burnin, pdim = pdim, ncspace=ncspace, B=B)
  #pbb3 <- arma.garch.bootstrap(2*pbb-pbb2, burnin = burnin, pdim = pdim, ncspace=ncspace, B=B)
  res.gg<-cbind(fit$par, sqrt(diag(solve(Infomat))),pbb)#,pbb2,pbb3)
  colnames(res.gg)<-c("gg", "SDgg", "pbb")#, "pbb2")
  rownames(res.gg)<- par.names[1:8]
  c(res.gg)
}
