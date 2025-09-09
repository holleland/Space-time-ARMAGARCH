# -- This code is not optimized and may take a very long time to run! --
setwd("C:/Users/s15052/Documents/Space-time-ARMAGARCH")
# -- loading data --
load("Data/SSTA_spatially_differenced_without_land.RData")

# -- Estimation --
m <- c(20, 14)
library(starmagarch)
library(spdep)
library(TMB)
W <- create.neighbourhood.array(m = m, sp = 2, type = "queen", torus =TRUE)




phi <- matrix(0, ncol = 1, nrow = 2)
theta <- matrix(0, ncol =1, nrow = 2)
alpha <- matrix(0, ncol = 1, nrow = 2)
beta <- matrix(0, ncol = 1, nrow =2)
omega <- 0.1

initial.parameters <- list(
  mu    =  mean(dat.m),
  phi   = phi,
  theta = theta,
  omega = omega,
  alpha = alpha,
  beta  = beta
)
map <- list(
  beta <- matrix(factor(c(1,NA)), ncol = 1, nrow =2)
)

#y <- apply(ww, 3, function(x) {cbind(1:length(c(x)), reshape2::melt(x)$value)})
#str(y)
f <- CreateLikelihood(dat.m, W=W,
                      init = apply(dat.m,1,var), parameters=initial.parameters, map = map)

fit <- fitSTARMAGARCH(f, data = dat.m, print = FALSE)
plot(fit)
plot_garch(fit)

#x <- f$report()$x
#h <- f$report()$sigma
# -------------------------------------------------
# --- If the sandwich estimator should be used: ---
# -------------------------------------------------
if(FALSE){
f2 = MakeADFun(data=list(y=dat.m, W=W, W2=W2,
                        init = apply(dat.m, 1, sd)^2),
              parameters=initial.parameters, 
              map = map, DLL = "STARMAGARCH_paper",
              ADreport = TRUE)
ell.gr <- f2$gr(fit$par)
str(ell.gr)
(sandwich <- sqrt(diag(solve(f$he(fit$par))%*%matrix(
  rowSums(apply(ell.gr, 1, function(x) x%*%t(x))), 
  ncol = length(fit$par))%*%solve(f$he(fit$par)))))

res$sand <- sandwich
res$Tscore.sand <- res$V1/res$sand
res$Pvalue.sand <- 2*pnorm(abs(res$Tscore.sand), lower.tail=FALSE )
}
rownames(res) <- c("$\\phi_{0,1}$", "$\\phi_{1,2}$", "$\\phi_{5,0}$", "$\\phi_{11,0}$","$\\phi_{12,1}$",
                   "$\\eta_{2,1}$", "$\\eta_{3,1}$","$\\eta_{11,0}$", "$\\eta_{12,1}$",
                   "$\\omega$", "$\\alpha_0$", "$\\alpha_1$", "$\\alpha_2$", "$\\alpha_3$", 
                   "$\\beta_0$", "$\\beta_1$")
#colnames(res)[c(1,5:7)] <- c("Estimates", "Standard Error", "T-score", "P-value")
colnames(res)<- c("Estimates", "Standard Error", "T-score", "P-value")
res[,1] <- paste("$", round(res[, 1]*c(1e1,1e1,1e2,1e1,1e1,1e1,1e1,1e1,1e1,1e4,1e2,1e2,1e2,1e2,1e1,1e1),2),
                 "\\cdot 10^{-",log10(c(1e1,1e1,1e2,1e1,1e1,1e1,1e1,1e1,1e1,1e4,1e2,1e2,1e2,1e2,1e1,1e1)), "}$", sep="")
res[,2] <- paste("$", round(res[,2]*c(1e3,1e2,1e3,1e3,1e2,1e3,1e3,1e3,1e2,1e5,1e3,1e3,1e3,1e3,1e2,1e2),2),
                 "\\cdot 10^{-",c(3,2,3,3,2,3,3,3,2,5,3,3,3,3,2,2), "}$", sep="")
#res[,5] <- paste("$", round(res[,5]*c(1e3,1e2,1e3,1e3,1e2,1e3,1e3,1e3,1e2,1e4,1e3,1e3,1e3,1e3,1e1,1e1),2),
#                 "\\cdot 10^{-",c(3,2,3,3,2,3,3,3,2,4,3,3,3,3,1,1), "}$", sep="")

#tab.for.print <- cbind(res[1:9,c(1,5)], c(row.names(res)[10:16], NA, NA), rbind(as.matrix(res[10:16, c(1,5)]), matrix(NA, ncol = 2, nrow= 2)))
tab.for.print <- cbind(res[1:9,c(1,2)], c(row.names(res)[10:16], NA, NA), rbind(as.matrix(res[10:16, c(1,2)]), matrix(NA, ncol = 2, nrow= 2)))
names(tab.for.print)[3]<- ""
library(xtable)
print(xtable(tab.for.print,
             caption = "Estimation results for SST data.", label = "tab:SST_example"),
      sanitize.text.function = function(x)x)
