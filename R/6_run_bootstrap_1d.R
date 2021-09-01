# --Adjust the spatial sample size
# spatial.size %in% c(5,20,50)
spatial.size <- 5

# ----------------------------
pcks<- c("TMB", 
         "Rcpp", "doParallel")
lapply(pcks, require, character.only = TRUE)
compile("Cpp/STARMAGARCH1D_CIRCULAR.cpp")
Rcpp::sourceCpp('Cpp/STARMAGARCH_1d_simulation.cpp')
#compile("Cpp/STARMA1D_CIRCULAR.cpp")       

#true.par <- c(5,.35,-.09,.12,-.05,.5,.13,.1)
true.par <- c(5.032,.354,-.153,.122,-.051,.513,.142,.171)
par.names <- c("mu","phi1","phi2","th1","th2","omg","a","b")
names(true.par)<-par.names
mu <- 5
cores <- detectCores()-1

source("R/bootstrap_1d.R")
cl <- makeCluster(cores)
clusterEvalQ(cl,expr = {library(TMB)
  dyn.load(dynlib("Cpp/STARMAGARCH1D_CIRCULAR"))
  Rcpp::sourceCpp('Cpp/STARMAGARCH_1d_simulation.cpp')
  source("R/bootstrap_1d.R")})

(t1<-Sys.time())
t <- rep(mu, 1000)
sim.res <- parSapply(cl, t,run.arma.garch.bootstrap, B=200, pdim = c(spatial.size, 200), par = true.par)
cat(t2<- Sys.time()-t1,"\n")
save(sim.res, file=paste("Data/bootstrap_1d_1000_B200_pdim", spatial.size, "200.RData", sep = ""))
stopCluster(cl)

boot <- t(apply(sim.res, 2, function(x) 2*x[1:8]-x[17:24]))
sdres <- apply(rbind(sim.res[c(1:8),],t(boot)),1,sd)

coverage <- matrix(rowMeans(apply(rbind(sim.res,t(boot)), 2, function(x){
  x[c(1:8,25:32)]-qnorm(.975)*sdres < rep(true.par,times=2) & 
    x[c(1:8,25:32)]+qnorm(.975)*sdres>rep(true.par,times=2)
  })), ncol=8,byrow=TRUE)


matrix(c(apply(sim.res[c(1:8,17:24),],1,sd),rowMeans(sim.res[9:16,])),byrow=TRUE,ncol=8)
org <- t(sim.res[1:8,])
res <- cbind(org, boot)
mres<-matrix(c(true.par,rowMeans(sim.res[1:8,]),sdres[1:8],rowMeans(sim.res[1:8,])-true.par,
               (rowMeans(sim.res[1:8,])-true.par)/true.par,coverage[1,],
               colMeans(boot),apply(boot,2,sd),colMeans(boot)-true.par, 
               (colMeans(boot)-true.par)/true.par,
               coverage[2,]), ncol = 8, byrow=TRUE)
colnames(mres)<-colnames(org)<-
  colnames(boot)<-c("$\\mu$","$\\phi_1$", "$\\phi_2$", "$\\eta_1$", "$\\eta_2$",
                                                                 "$\\omega$", "$\\alpha$","$\\beta$")
rownames(mres)<-c("$\\theta_0$","$\\widehat\\theta_n$","SD","Bias","RBias","Cover",
                  "$\\widehat\\theta_{\\text{BC}}$","SD","Bias","RBias", "Cover")
mres2 <- data.frame(s = " ", names = as.character(rownames(mres)), mres)
colnames(mres2)<-c("s", "r", colnames(mres))
rownames(mres2)<-1:nrow(mres2)
library(xtable)
print(xtable(mres2, digits=3, align = c("r", rep("c",10))),
      #caption = "Simulated sample: $130\\times 1200$. Truncated sample: $30\\times200$. "),
      include.rownames=FALSE, sanitize.text.function = function(x){x}, 
      file = paste("Tex/2_table2_pbbc_200boot_", spatial.size, "x200.tex", sep = ""))
