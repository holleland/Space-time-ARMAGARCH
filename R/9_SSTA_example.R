# -- This code is not optimized and may take a very long time to run! --

# -- loading data --
load("Data/SSTA_spatially_differenced_without_land.RData")

# -- Estimation --
m <- c(20, 14)
library(starmagarch)
W <- create.neighbourhood.array(m = m, sp = 3, type = "queen", torus =TRUE)

W2 <- array(0, dim = c(prod(m), prod(m), 4))
W2[,,1] <- diag(prod(m))
A1 <- nb2mat(nblag(cell2nb(m[2], 1, torus = TRUE, type ="rook"),2)[[1]])*4-diag(m[2])
for(j in 1:m[1])
  W2[(j-1)*m[2]+1:m[2],(j-1)*m[2]+1:m[2],2] <- A1/sum(A1[1,])
B2 <- create.neighbourhood.array(m = m, sp = 3, type = "rook", torus =TRUE)[,,2]*2
W2[,,3] <- B2-W2[,,2]
B3 <- create.neighbourhood.array(m = m, sp = 3, type = "queen", torus =TRUE)[,,2]*8
W2[,,4] <- (B3 - W2[,,3]*2 - W2[,,2]*2)/4

phi <- matrix(0, ncol = 12, nrow = 3)
theta <- matrix(0, ncol = 12, nrow = 2)
alpha <- matrix(0, ncol = 1, nrow = 4)
beta <- matrix(0, ncol = 1, nrow =2)

phi[1,1] <- 0.361
phi[3,1] <- 0.225
phi[1,5] <- 0.026
phi[1,11]<- 0.360
phi[2,12]<--0.248
theta[2,2] <- 0.163
theta[2,3] <- 0.143
theta[1,11]<--0.378
theta[2,12]<- 0.205

omega <- 5.06 * 1e-4
alpha[1,1] <- 6.25 * 1e-2
alpha[2,1] <- 6.57 * 1e-2
alpha[3,1] <- 3.68 * 1e-2
alpha[4,1] <- 4.05 * 1e-2

beta[1,1] <- 2.66 * 1e-1
beta[2,1] <- 5.09 * 1e-1
A1 <- W2[,,2]
W2[,,2] <- W2[,,3]
W2[,,3] <- A1

initial.parameters <- list(
  mu    = 0,
  phi   = phi,
  theta = theta,
  omega = omega,
  alpha = alpha,
  beta  = beta
)
#theta[2,12]<- 0
phi.map <- matrix(1:prod(dim(phi)), ncol = ncol(phi))
theta.map <- matrix(prod(dim(phi))+1:prod(dim(theta)), ncol = ncol(theta))
phi.map[phi==0] <- NA
theta.map[theta==0] <- NA

map <- list(
  mu = factor(NA),
  phi = factor(phi.map),
  theta = factor(theta.map)
  
)
compile("Cpp/STARMAGARCH_paper.cpp")
dyn.load(dynlib("Cpp/STARMAGARCH_paper"))

f = MakeADFun(data=list(y=dat.m, W=W, W2=W2,
                        init = apply(dat.m, 1, sd)^2),
              parameters=initial.parameters, 
              map = map, DLL = "STARMAGARCH_paper")
fit = nlminb(f$par,f$fn,f$gr, f$he, 
             lower=c(rep(-5,length(which(names(f$par) %in% c("phi","theta")))),
                     rep(1e-20, length(which(names(f$par)%in% c("omega","alpha","beta"))))),
             upper=c(rep(5,length(f$par))))

sd <- sqrt(diag(solve(f$he(fit$par))))
res <- as.data.frame(cbind(fit$par, sd))
res$Tscore <- res$V1/res$sd
res$Pvalue <- 2*pnorm(abs(res$V1/res$sd), lower.tail=FALSE )

#x <- f$report()$x
#h <- f$report()$sigma
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

rownames(res) <- c("$\\phi_{0,1}$", "$\\phi_{1,2}$", "$\\phi_{5,0}$", "$\\phi_{11,0}$","$\\phi_{12,1}$",
                   "$\\eta_{2,1}$", "$\\eta_{3,1}$","$\\eta_{11,0}$", "$\\eta_{12,1}$",
                   "$\\omega$", "$\\alpha_0$", "$\\alpha_1$", "$\\alpha_2$", "$\\alpha_3$", 
                   "$\\beta_0$", "$\\beta_1$")
colnames(res)[c(1,5:7)] <- c("Estimates", "Standard Error", "T-score", "P-value")
res[,1] <- paste("$", round(res[, 1]*c(1e1,1e1,1e2,1e1,1e1,1e1,1e1,1e1,1e1,1e4,1e2,1e2,1e2,1e2,1e1,1e1),2),
                 "\\cdot 10^{-",log10(c(1e1,1e1,1e2,1e1,1e1,1e1,1e1,1e1,1e1,1e4,1e2,1e2,1e2,1e2,1e1,1e1)), "}$", sep="")
res[,5] <- paste("$", round(res[,5]*c(1e3,1e2,1e3,1e3,1e2,1e3,1e3,1e3,1e2,1e4,1e3,1e3,1e3,1e3,1e1,1e1),2),
                 "\\cdot 10^{-",c(3,2,3,3,2,3,3,3,2,4,3,3,3,3,1,1), "}$", sep="")
tab.for.print <- cbind(res[1:9,c(1,5)], c(row.names(res)[10:16], NA, NA), rbind(as.matrix(res[10:16, c(1,5)]), matrix(NA, ncol = 2, nrow= 2)))
names(tab.for.print)[3]<- ""
library(xtable)
print(xtable(tab.for.print,
             caption = "Estimation results for SST data.", label = "tab:SST_example"),
      sanitize.text.function = function(x)x)
