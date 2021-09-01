# ------------------------
# --  CSTA   vs   CSTAG --
# ------ experiment ------
# ------------------------
pcks<- c("TMB", "doParallel")
lapply(pcks, require, character.only = TRUE)
compile("Cpp/STARMAGARCH1D_CIRCULAR.cpp")
compile("Cpp/STARMA1D_CIRCULAR.cpp")       

run.arma.garch.comparison <- function(true.par,burnin = 1000, pdim = c(30,1500),
                                      seed = NULL){
  # -- Parameters:  --
  mu   = true.par[1]
  phi1 = true.par[2]
  phi2 = true.par[3]
  th1  = true.par[4]
  th2  = true.par[5]
  omg  = true.par[6]
  a    = true.par[7]
  b    = true.par[8]
  
  par.names <- c("mu","phi1","phi2","th1","th2","omg","a","b","sigma")
  
  if(!is.null(seed))
  set.seed(seed)
  
  # -- Simulation --
  # x  - ARMA-GARCH
  # x2 - ARMA
  
  x <-x2<- u<- h<- array(0, dim=pdim)
  z <- array(rnorm(prod(pdim)), dim=pdim)
  z2 <- sqrt(omg/(1-sum(a,b)*3))*z
  
  for(t in 1:(pdim[2]-1)){
    for(k in 1:pdim[1]){
      h[k,t+1]<- sqrt(omg + a * (u[(k-2)%%pdim[1]+1,t]^2+u[k,t]^2+u[k%%pdim[1]+1,t]^2)+
                        b * (h[(k-2)%%pdim[1]+1,t]^2+h[k,t]^2+h[k%%pdim[1]+1,t]^2))
      u[k,t+1]<- h[k,t+1]*z[k,t+1]
    }
  }
  
  for(t in 1:(pdim[2]-1)){
    for(k in 1:pdim[1]){
      x[k,t+1]<- phi1*x[k,t]+phi2*(x[(k-2)%%pdim[1]+1,t]+x[k%%pdim[1]+1,t])+
        th1*u[k,t]+th2*(u[(k-2)%%pdim[1]+1,t]+u[k%%pdim[1]+1,t])+u[k,t+1]
      x2[k,t+1]<- phi1*x2[k,t]+phi2*(x2[(k-2)%%pdim[1]+1,t]+x2[k%%pdim[1]+1,t])+
        th1*z2[k,t]+th2*(z2[(k-2)%%pdim[1]+1,t]+z2[k%%pdim[1]+1,t])+z2[k,t+1]
    }
  }
  
  # -- Removing burn-in --
  x<-x[,-(1:burnin)];x2<-x2[,-(1:burnin)];h<-h[,-(1:burnin)];u<-u[,-(1:burnin)];
  z<-z[,-(1:burnin)];z2<-z2[,-(1:burnin)];
  
  x  <- x  + mu
  x2 <- x2 + mu

  # -- Estimation   --
  # -- CSTARMAGARCH --
  
  # .. AG-AG ..
  dyn.load(dynlib("Cpp/STARMAGARCH1D_CIRCULAR")) 
  f <- MakeADFun(data = list(x=x, mut=0.0, ustart = rep(0,pdim[1])),
                parameters = list(mu = mean(as.vector(x)), phi1 = .25, phi2=-0.1, theta1=0.1, theta2=0.1,
                                  omega=.33, alpha=.02, beta=.02), silent =TRUE, DLL = "STARMAGARCH1D_CIRCULAR")
  fit <- nlminb(f$par,f$fn,f$gr, f$he, lower=c(-200,
                                               -5,-5,
                                               -5,-5,
                                               0.0001,0.0001,0.0001),
                upper=c(500,
                        5,5,
                        5,5,
                        1000,5,5))
  Infomat <- f$he(fit$par)
  res.gg  <- cbind(fit$par, sqrt(diag(solve(Infomat))))
  colnames(res.gg) <- c("gg", "SDgg")
  rownames(res.gg) <- par.names[1:8]
  
  # .. A-AG ..
  
  f <- MakeADFun(data=list(x=x2, mut=0.0, ustart = rep(0,pdim[1])),
                          parameters=list(mu = mean(as.vector(x2)), phi1 = .25, phi2=-0.1, 
                                          theta1=0.1, theta2=0.1,
                                          omega=.33, alpha=.02, beta=.02), silent =TRUE,
                DLL = "STARMAGARCH1D_CIRCULAR")
  fit <- nlminb(f$par,f$fn,f$gr, f$he, lower=c(-200,
                                               -5,-5,
                                               -5,-5,
                                               0.0001,0.0001,0.0001),
                upper = c(500,
                          5,5,
                          5,5,
                          1000,5,5))
  Infomat2 <- f$he(fit$par)
  res.ag   <- cbind(fit$par, sqrt(diag(solve(f$he(fit$par)))))
  colnames(res.ag) <- c("ag", "SDag")
  rownames(res.ag) <- par.names[1:8]
  
  
  dyn.unload(dynlib("Cpp/STARMAGARCH1D_CIRCULAR")) 
  
  # -- CSTARMA --
  
  # .. AG-A ..
  dyn.load(dynlib("Cpp/STARMA1D_CIRCULAR")) 
  f <- MakeADFun(data=list(x=x, mut=0.0, y =x*0, gam1=0.0),
                 parameters=list(mu = mean(as.vector(x)),  
                                 phi1 = .25, phi2=-.1, theta1=0.1, theta2=0.1, #gam1=0.1, 
                                 sigma = 5), silent =TRUE,
                 DLL = "STARMA1D_CIRCULAR")
  fit.arma <- nlminb(f$par,f$fn,f$gr, f$he, lower=c(-200,-5,-5,-5,-5,0.1),
                     upper=c(500,5,5,5,5,
                             5000))
  
  theta.arma <- fit.arma$par
  # ---------------------------------
  # -- SANDWICH FOR AG-A SITUATION --
  # ---------------------------------
  
  Infomat3 <- f$he(fit.arma$par)
  if(FALSE){
  f = MakeADFun(ADreport = TRUE, data=list(x=x, mut=0.0, y =x*0, gam1=0.0),
                parameters=list(mu = mean(as.vector(x)),  
                                phi1 = .25, phi2=-.1, theta1=0.1, theta2=0.1, #gam1=0.1, 
                                sigma = 5), silent =TRUE,
                DLL = "STARMA1D_CIRCULAR")
  
  ell.gr.A <- matrix(rowSums(apply(f$gr(theta.arma), 1, function(x) c(x%*%t(x)))), ncol = length(theta.arma))
  }
  res.ga   <- cbind(fit.arma$par, sqrt(diag(solve(Infomat3))))#sqrt(diag(solve(Infomat3)%*% ell.gr.A %*% solve(Infomat3))))
  colnames(res.ga) <- c("ga", "SDga")
  rownames(res.ga) <- par.names[c(1:5,9)]
  
  # .. A-A ..
  f <- MakeADFun(data=list(x=x2, mut=0.0, y =x*0, gam1=0.0),
                 parameters=list(mu = mean(as.vector(x)),  
                                 phi1 = .25, phi2=-.1, theta1=0.1, theta2=0.1, #gam1=0.1, 
                                 sigma = 5), silent =TRUE,
                 DLL = "STARMA1D_CIRCULAR")
  fit.arma <- nlminb(f$par,f$fn,f$gr, f$he, lower=c(-200,-5,-5,-5,-5,0.1),
                     upper=c(500,5,5,5,5,5000))
  
  theta.arma <- fit.arma$par
  Infomat4   <- f$he(fit.arma$par)
  res.aa     <- cbind(theta.arma, sqrt(diag(solve(Infomat4))))
  
  colnames(res.aa) <- c("aa", "SDaa")
  rownames(res.aa) <- par.names[c(1:5,9)]
  dyn.unload(dynlib("Cpp/STARMA1D_CIRCULAR")) 
  
  # ------------------------
  # -- Combining results: --
  # ------------------------
  resG <- c(res.gg,res.ag)
  names(resG) <- paste(c(rep("GG-QMLE",8),rep("GG-SD",8), rep("AG-QMLE",8),rep("AG-SD",8)),
                       rep(par.names[1:8],4), sep = "-")
  resA <- c(res.ga, res.aa)
  names(resA) <- paste(c(rep("GA-QMLE",6),rep("GA-SD",6), rep("AA-QMLE",6),rep("AA-SD",6)),
                       rep(par.names[c(1:5,9)],4), sep = "-")
  c(mean(x),mean(x2), resG, resA)
}

# --------------------
# -- RUN EXPERIMENT --
# --------------------
cl <- makeCluster(detectCores()-1)
cl <- makeCluster(20)
clusterEvalQ(cl,expr = {
    library(TMB) 
    })
true.par <- c(5.032,.354,-.153,.122,-.051,.513,.142,.171)
t1<-Sys.time()
seeds <- 1:2000
#run.arma.garch.comparison(true.par, pdim = c(30,1200)) #test run
sim.res <- parSapply(cl, seeds, run.arma.garch.comparison, true.par=true.par, burnin = 1000,
                     pdim = c(30, 1200))
cat(t2<- Sys.time()-t1,"\n")

stopCluster(cl)

# -----------------
# -- Make table: --
# -----------------

means<-rowMeans(sim.res[-(1:2),])

A <- matrix(means[1:32],ncol=8,nrow=4, byrow=TRUE) # GG and AG 
B <- matrix(means[33:(length(means))], ncol = 6, nrow=4,byrow=TRUE) # GA and AA
C <- matrix(NA_real_, ncol = 9, nrow= 8)
C[1:4,1:8]      <- A
C[5:8,c(1:5,9)] <- B
C               <- rbind(c(true.par,sqrt(true.par[6]/(1-sum(true.par[7:8])*3))),C)
C[seq(3,9,2),]  <- sqrt(200*50)*C[seq(3,9,2),]

colnames(C) <-  c("$\\mu$","$\\phi_1$", "$\\phi_2$", "$\\eta_1$", "$\\eta_2$",
                  "$\\omega$", "$\\alpha$","$\\beta$",
                  "$\\sigma$")
rownames(C) <- c("$\\theta_0$","QMLE-GG", "$\\sqrt{N}$SD-GG", "QMLE-AG", 
                 "$\\sqrt{N}$SD-AG","QMLE-GA", "$\\sqrt{N}$SD-GA", "QMLE-AA", "$\\sqrt{N}$SD-AA")
A       <- C[4:5,]
C[4:5,] <- C[6:7,]
C[6:7,] <- A
rownames(C)[4:7] <- rownames(C)[c(6:7,4:5)]
D       <- c(true.par[1],mean(sim.res[1,]),sqrt(200*30)*sd(sim.res[1,]),NA,NA,NA,NA, 
             mean(sim.res[2,]),sqrt(200*30)*sd(sim.res[2,]))
C[c(2,6),9] <- sqrt(C[c(2,6),6]/(1- 3*(C[c(2,6),7]+C[c(2,6),8])))

# -- Making table: --
library(xtable)
cap <- "Maximum likelihood estimates with corresponding standard deviations of 
CSTA with either iid WN$(0, \\sigma^2)$ (\\iidwn) or GARCH (\\garchwn) errors, 
where $\\sigma$ is set by \\autoref{eq:unconditonal_variance}. The values are 
means of 2000 repeated simulations of sample size $30\\times200$ ($N=6000$). 
The blue $\\sigma$ values are derived from the estimated GARCH parameters 
through \\autoref{eq:unconditonal_variance}, i.e. 
$\\{\\hat\\omega/(1-3\\hat\\alpha-3\\hat\\beta)\\}^{1/2}$."
print(xtable(C, digits = 3, caption = cap, 
             label = "tab:meansim_GARCHvsARMA"),
      hline.after = c(-1,0,1,3,5,7,9),sanitize.text.function = function(x){x}, 
      file = "Tables/resG_mean_of_simulations_6k_4juni.tex")

print(C)

str(sim.res)
