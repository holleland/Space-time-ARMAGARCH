# -------------------------------------------------------
# ---------------- Running experiements -----------------
# -------------------------------------------------------

# -------------------------
# -- CIRCULAR EXPERIMENT --
# -------------------------
library(spdep)
library(TMB)
library(Rcpp)
library(starmagarch)
library(doParallel)
batch.size <- 1#250
#compile("Cpp/STARMAGARCH.cpp")
#dyn.load(dynlib("Cpp/STARMAGARCH"))
# --------------------
true.par  <- c(5, 
               .25,.0625,
               -.15,-.3/4,
               .5, .3,.12/4,
               .15,.10/4)
W <- create.neighbourhood.array(m = c(5,5), sp = 2, type = "rook", torus = TRUE, sum.to.one = FALSE)
parameters = list(
  mu = true.par[1],
  phi = matrix(true.par[2:3], ncol = 1),
  theta = matrix(true.par[4:5], ncol = 1),
  omega = true.par[6],
  alpha= matrix(true.par[7:8], ncol = 1),
  beta = matrix(true.par[9:10], ncol = 1)
)
run_simest <- function(i=1, parameters, W){
  x <- simSTARMAGARCH(parameters=parameters,
                    n = 1000,
                    m = c(5,5),
                    burnin=1000, 
                    W=W)
  f <- CreateLikelihood(data = x, W=W, parameters = list(
        mu = mean(x), 
        phi = matrix(c(.2,.05), ncol = 1), 
        theta = matrix(c(-.1,-.05), ncol = 1),
        omega = 0.4, 
        alpha = matrix(c(.2,.01), ncol = 1),
        beta = matrix(c(.1,.01), ncol = 1)
  ), silent = TRUE)
  return(fitSTARMAGARCH(f, data = x, simple=TRUE))
}
# test 
run_simest(1, parameters = parameters, W=W)

cores <- 1#detectCores()-1
cl <- makeCluster(cores)
clusterEvalQ(cl,expr = {
  library(TMB);library(Rcpp); library(starmagarch)
  #source("STARMAGARCH/estimation.R")
  #source("STARMAGARCH/simulation.R")
  #dyn.load(dynlib("Cpp/STARMAGARCH"))
})

# -- Batch 1: --
(t1 <- Sys.time())
cres1 <- parSapply(cl, 1:batch.size, run_simest, parameters = parameters, W=W)
cat("Circular batch 1 of 4 took:\t ", Sys.time()-t1,"\n")
save(cres1, file = "Data/cres_5x5x1000_with_sd_1.RData")

# -- Batch 2: --
(t2 <- Sys.time())
cres2 <- parSapply(cl, 1:batch.size, run_simest, parameters = parameters, W=W)
cat("Circular batch 2 of 4 took:\t ", Sys.time()-t2,"\n")
save(cres2, file = "Data/cres_5x5x1000_with_sd_2.RData")

# -- Batch 3: --
(t3 <- Sys.time())
cres3 <- parSapply(cl, 1:batch.size, run_simest, parameters = parameters, W=W)
cat("Circular batch 3 of 4 took:\t ", Sys.time()-t3,"\n")
save(cres3, file = "Data/cres_5x5x1000_with_sd_3.RData")

# -- Batch 4: --
(t4 <- Sys.time())
cres4 <- parSapply(cl, 1:batch.size, run_simest, parameters = parameters, W=W)
cat("Circular batch 4 of 4 took:\t ", Sys.time()-t4,"\n")
save(cres4, file = "Data/cres_5x5x1000_with_sd_4.RData")

stopCluster(cl)

cres <- cbind(cres1, cres2, cres3, cres4)

# ------------------------------------------
# -- Parametric Bootstrap Bias Correction --
# ------------- Experiment -----------------
# ------------------------------------------

W <- create.neighbourhood.array(m = c(5,5), sp = 2, type = "rook", torus = TRUE, sum.to.one = FALSE)
Rcpp::sourceCpp('Cpp/STARMAGARCH_2d_simulation_5neighbours.cpp')
cl <- makeCluster(cores)
clusterEvalQ(cl,expr = {
  library(TMB);library(Rcpp); library(starmagarch)
  #On Loke: 
  Rcpp::sourceCpp('Cpp/STARMAGARCH_2d_simulation_5neighbours.cpp')
  #source("STARMAGARCH/estimation.R")
  #source("STARMAGARCH/simulation.R")
  #dyn.load(dynlib("Cpp/STARMAGARCH"))
  boot<-function(j=1, par, dim=c(5,5,1000), W){
   x <- simulate_2d(par = par[c(1:3,3,4:5,5,6:8,8,9:10,10)], dim = dim, burnin_time = 1000, burnin_space = 50)
   f <- CreateLikelihood(data = apply(x, 3, c), W=W, parameters = list(
      mu = mean(x), 
      phi = matrix(c(.2,.05), ncol = 1), 
      theta = matrix(c(0.01,.01), ncol = 1),
      omega = 0.3, 
      alpha = matrix(c(.2,.01), ncol = 1),
      beta = matrix(c(.1,.01), ncol = 1)
    ), silent = TRUE)
  fitSTARMAGARCH(f, data = x2, simple=TRUE)[1:10]
}
})
# bootstrap 2.0
run_simest <- function(i=1, par, W, B = 200){
  thetao <- boot(par = par, W = W)
  thetab <- 2*thetao - rowMeans(sapply(1:B, boot, par = thetao, W=W))
  return(c(thetao,thetab))
}

# test
run_simest(1, true.par, W, B =2)

# -- Batch 1: --
(t1 <- Sys.time())
bres1 <- parSapply(cl, 1:batch.size, run_simest, par = true.par, W=W, B = 200)
cat("PBBC batch 1 of 4 took:\t ", Sys.time()-t1,"\n")
save(bres1, file = "Data/bres_5x5x1000_B200_1.RData")

# -- Batch 2: --
(t2 <- Sys.time())
bres2 <- parSapply(cl, 1:batch.size, run_simest, par = true.par, W=W, B = 200)
cat("PBBC batch 2 of 4 took:\t ", Sys.time()-t2,"\n")
save(bres2, file = "Data/bres_5x5x1000_B200_2.RData")

# -- Batch 3: --
(t3 <- Sys.time())
bres3 <- parSapply(cl, 1:batch.size, run_simest, par = true.par, W=W, B = 200)
cat("PBBC batch 3 of 4 took:\t ", Sys.time()-t3,"\n")
save(bres3, file = "Data/bres_5x5x1000_B200_3.RData")

# -- Batch 4: --
(t4 <- Sys.time())
bres4 <- parSapply(cl, 1:batch.size, run_simest, par = true.par, W=W, B = 200)
cat("PBBC batch 4 of 4 took:\t ", Sys.time()-t4,"\n")
save(bres4, file = "Data/bres_5x5x1000_B200_4.RData")

# ------------

stopCluster(cl)

bres <- cbind(bres1, bres2, bres3, bres4)

# ---------------------------------------------------
# --------------------- Analysis --------------------
# ---------------------------------------------------

# -------------------
# -- MAKING TABLES --
# -------------------

# -- Circular-circular --
cmeans<-rowMeans(cres[1:10,])
cbias <- cmeans-true.par
csd <- apply(cres[1:10,],1,sd)

means <- rowMeans(bres)
bias.each <- apply(bres,2, function(x)x-rep(true.par,2))
bias <- means-rep(true.par,2)
rowMeans(apply(bres,2,function(x)x-rep(true.par,2)))-bias
sdres<-apply(bres,1, sd)
bias.mat <- rbind(bias[1:10], bias[11:20], cbias)
apply(bias.mat, 2, function(x)which.min(abs(x)))

coverage <- matrix(rowMeans(apply(rbind(cres[1:10,],bres), 2, function(x){
  x-qnorm(.975)*c(csd, sdres) < rep(true.par,times=3) & 
    x+qnorm(.975)*c(csd, sdres)>rep(true.par,times=3)
})), ncol=10,byrow=TRUE)
crmse <- sqrt(cbias^2+csd^2)
rmse <- sqrt(bias^2+sdres^2)




res<-matrix(c(true.par,
              cmeans,csd,cbias,crmse, cbias/true.par, coverage[1,],
              means[1:10],sdres[1:10],bias[1:10],rmse[1:10],bias[1:10]/true.par, coverage[2,],
              means[11:20],sdres[11:20],bias[11:20],rmse[11:20],bias[11:20]/true.par, coverage[3,]),
            byrow=TRUE, ncol = 10)



colnames(res)<-paste("$\\",c("mu", "phi_{0}","phi_{1}","eta_{0}", "eta_{1}", "omega", "alpha_{0}", "alpha_{1}", "beta_{0}", "beta_{1}"),"$", sep="")
rownames(res)<-c("$\\lambda_0$","Est.","SD CQMLE","Bias CQMLE","CRMSE","RBias CQMLE","Cover",
                 "1Est.", "1SD","1Bias","1RMSE","1RBias","Cover1",
                 "2Est.","2SD", "2Bias","2RMSE", "2RBias","Cover2")
round(res,3)

library(xtable)
print(xtable(res, digits =3,
             caption="Mean circular estimation results of 1000 repeated simulations of a two dimensional STARMA$(1,1)$\\,-\\,GARCH$(1,1)$ process based on $5\\times 5\\times 1000$ observations."),
      sanitize.text.function = function(x)x,
      file="Tex/3_boot_2D_5x5x1000_RMSE.tex")

# -------------------
# -- MAKING FIGURE --
# -------------------
library("reshape2")
library("ggplot2")
QMLE <- bres[c(1:10),]
PBBC <- bres[c(11:20),]
CQMLE <- cres[c(1:10),]
par.names <- c("mu","phi[0]","phi[1]","eta[0]","eta[1]",
               "omega","alpha[0]","alpha[1]","beta[0]","beta[1]")
df.truepar <- data.frame(par = par.names, value = true.par)
rownames(QMLE)<- par.names
rownames(PBBC)<- par.names
rownames(CQMLE)<- par.names
PBBC.melt <- melt(t(PBBC))
names(PBBC.melt)<-c("sim", "par", "est")
PBBC.melt$type <- "Non-circular data with PBBC"
QMLE.melt <- melt(t(QMLE))
names(QMLE.melt)<-c("sim", "par", "est")
QMLE.melt$type <- "Non-circular data"
CQMLE.melt <- melt(t(CQMLE))
names(CQMLE.melt)<-c("sim", "par", "est")
CQMLE.melt$type <- "Circular data"

df.res <- rbind(QMLE.melt, PBBC.melt, CQMLE.melt)
colorvalues <- c(rgb(142, 169, 105, maxColorValue = 255),
                 rgb(221,188,75,maxColorValue= 255),
                 rgb(78,160,183,maxColorValue = 255),
                 rgb(229, 79, 71,maxColorValue = 255),
                 rgb(138,103,149,maxColorValue= 255))

ggplot(data = df.res, aes(x = est))+stat_density(aes(color = type, fill = type))+
  scale_fill_manual(values=colorvalues)+
  scale_color_manual(values=colorvalues)+
  facet_wrap(~par, ncol = 5, scales="free", labeller = label_parsed)+
  theme_minimal()+xlab(NULL)+ylab("Density")+theme(strip.text = element_text(size = 14),
                                                   legend.title=element_blank(),
                                                   legend.position = "bottom")+
  geom_vline(data=df.truepar, aes(xintercept = value), lty=2, col = "red")+
  theme(#plot.background = element_rect(fill = rgb(242,237,234, maxColorValue = 255),
    #                               color = rgb(242,237,234, maxColorValue = 255)),
    panel.grid = element_blank(),
    strip.text = element_text(color = rgb(41,52,81,max=255), size = 16),
    axis.text = element_text(color = rgb(41,52,81,max=255)),
    axis.title = element_text(color = rgb(41,52,81,max=255), face="bold", size =14),
    legend.text = element_text(color = rgb(41,52,81,max=255), size = 14))

ggsave("Figures/7_pbbc_comparison.pdf", width = 12, height = 5)

