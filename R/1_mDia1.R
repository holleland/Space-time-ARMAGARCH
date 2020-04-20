library(reshape2)
library(ggplot2)
library(spdep)
library(starma)
source("R/Functions.R")


df <- read_cell("Data/Chan2_1L_imActmap.txt") 

any(is.na(df))
sm <- rowMeans(df)
tm <- colMeans(df)

# ------------------- #
# --- Detrending ---- #
# ------------------- #
df.melt <- melt(df)
names(df.melt)[1:2] <- c("u","t")
p <- 4
lmfit<-lm(value ~poly(u,p),
          data = df.melt)
summary(lmfit)
#plot(sm)
u <- 1:nrow(df)
sfit <- predict(lmfit, newdata=poly(u,p))

# ------------------- #
# ----- plotting ---- #
# ------------------- #
means.df <- data.frame(means = c(sm,tm), x = c(1:nrow(df), (1:ncol(df))*10), 
                       type = c(rep("Space", nrow(df)),rep("Time (sec)", ncol(df))))
detrend.df <- data.frame(fit = sfit, x = 1:nrow(df),type = rep("Space", length(sfit)))
intercept.df <- data.frame(yint = lmfit$coefficients[1], type = "Time (sec)")
ggplot(means.df, aes(x, means))+geom_point(color = "blue")+
  facet_wrap(~type, ncol = 2, scales="free_x", 
             strip.position = "bottom")+
  theme_minimal()+
  theme(#plot.title = element_text(size = 20, face = "bold", hjust = .5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(colour = "black"))+
  xlab(NULL)+ylab("Means")+
  geom_line(data=detrend.df, aes(x, fit), color = "red", lwd = .6)+
  geom_hline(data=intercept.df, aes(yintercept = yint), color = "red", lwd = .6)
ggsave("Figures/9_mdia1_space-time-means.pdf", 
       width = 9, height = 4)

df.melt$fitted <- as.numeric(fitted(lmfit))


df.fit <- acast(df.melt, u~t, value.var = "fitted")

df.melt$residuals <- resid(lmfit)
df2<-acast(df.melt, u~t, value.var = "residuals")
plotting_df(df2, main = NULL,midpoint = 0,
            guides=FALSE, path = "Figures/9_mdia1_detrended.pdf")


# ----------------------
# -- Fitting CSTARMA: --
# ----------------------
#dyn.unload(dynlib("Cpp/STARMA"))
compile("Cpp/STARMA.cpp")
dyn.load(dynlib("Cpp/STARMA"))
parameters <- list(mu = 0.0, 
                   phi   = matrix(0, ncol = 3, nrow=2),
                   theta = matrix(0, ncol =2, nrow = 2),
                   sigma = sd(df))
                  # alpha = matrix(c(0.06), ncol = 1),
                  # beta  = matrix(c(0.82, 0.13), ncol = 1))
W <- create.neighbourhood.array(m = c(nrow(df), 1), sp = 4, torus = TRUE, type = "rook")
phimap <- matrix(c(1:6), ncol = 3, nrow = 2)
thetamap <- matrix(6+c(1:4), ncol = 2, nrow = 2)
phimap[4]<-NA_real_
thetamap[4]<-NA_real_
map <- list(mu = factor(NA),
            phi=factor(phimap),
           theta = factor(thetamap))
f <- MakeADFun(data = list(y=df2, W=W,init = apply(df2,1,var)),
               parameters = parameters, silent = TRUE, map = map, DLL = "STARMA")
lower <- c(rep(-5,length(which(names(f$par)%in% c("phi","theta")))),1e-8)
upper <- c(rep(5,length(which(names(f$par)%in% c("phi","theta")))),1e5)

fit <- stats::nlminb(f$par,f$fn,f$gr, f$he,
                     lower = lower,
                     upper = upper)
# ---------------
# -- SANDWICH: --
# ---------------
J <- f$he(fit$par)
#
# With ADreport = TRUE, TMB calculate the gradient of ell(t,u) (instead of L)
#
f <- MakeADFun(data = list(y=df2, W=W,init = apply(df2,1,var)),
               parameters = parameters, silent = TRUE, map = map,
               ADreport = TRUE, DLL = "STARMA")
S <- f$gr(fit$par)
Lhat <-matrix(rowSums(matrix(apply(S, 1, function(k) c(-k %*% -t(k))), 
                             nrow = length(fit$par)^2)), ncol = length(fit$par))
SIG <- solve(J) %*% Lhat %*% solve(J)

# ------------------------------
# -- Making table of results: --
# ------------------------------
armapar <- fit$par
matcoef <- data.frame(Estimates = fit$par,
                      SD = sqrt(diag(SIG)))
matcoef$Zscore <- matcoef$Estimates/matcoef$SD
matcoef$Pvalue <- 2*stats::pnorm(abs(matcoef$Zscore), lower.tail=FALSE)
rownames(matcoef)<-correct.names(names(f$par))
rownames(matcoef)<-c("$\\phi_{10}$","$\\phi_{11}$", paste("$\\phi_{", c(2,3,3),c(0,0:1),"}$", sep = ""),
                    "$\\eta_{10}$","$\\eta_{11}$",paste("$\\eta_{", c(2),c(0),"}$", sep=""), "$\\sigma$")
library(xtable)
print(xtable(matcoef, digits = 3,
             caption = "Parameter estimates with subscripts $jv$ where $j$ is temporal- and $v$ is spatial lag,  for a CSTARMA$(2,2)$ model based on the detrended mDia1 levels in the first layer of Cell no. 2."), 
      sanitize.text.function =function(x)x,
      file = "Tex/mdia1_ARMA_estimates.tex")

dyn.unload(dynlib("Cpp/STARMA"))
cat("CSTARMA complete.\n")

# ---------------------------
# -- Fitting CSTARMAGARCH: --
# ---------------------------
compile("Cpp/STARMAGARCH.cpp")
dyn.load(dynlib("Cpp/STARMAGARCH"))
parameters <- list(mu = 0.0,#mean(df2),
                   phi   = matrix(c(armapar[c(1:3)],0,armapar[4:5]), ncol = 3, nrow=2),
                   theta = matrix(c(armapar[c(6:8)],0), ncol = 2, nrow = 2),
                   omega = 400,
                   alpha = matrix(0, ncol = 1, nrow=2),
                   beta  = matrix(0, #c(0, 0,0,0,0,0,.5,0,0), 
                                  ncol=1,nrow=3))#c(,0,.04,rep(0,3),.05,0,0), ncol = 3, nrow=3))
W <- create.neighbourhood.array(m = c(nrow(df), 1), sp = 4, torus = TRUE, type = "rook")
alphamap <-matrix(1:2, ncol = 1, nrow =2)
alphamap[1]<-NA_real_
betamap <- matrix(12+1:3, ncol = 1, nrow =3)
betamap[c(2,3)]<-NA_real_
map = list(mu = factor(NA_real_),phi = factor(phimap),theta = factor(thetamap),
           #alpha = factor(alphamap),
           beta = factor(betamap))#,
f <- MakeADFun(data = list(y=df2, W=W,init = rep(armapar[length(armapar)]^2, nrow(df2))),#apply(u,1,var)),#rep(armapar[length(armapar)]^2, nrow(df2))),#apply(df2,1,var)),
               parameters = parameters, silent = TRUE, map = map, DLL = "STARMAGARCH")
lower <- c(rep(-5,length(which(names(f$par)%in% c("phi","theta")))),1e-8,
           rep(1e-8, length(which(names(f$par)%in% c("alpha","beta")))))
upper <- c(rep(5,length(which(names(f$par)%in% c("phi","theta")))),1e5,
           rep(1, length(which(names(f$par)%in% c("alpha","beta")))))

fit <- stats::nlminb(f$par,f$fn,f$gr, f$he,
                     lower = lower,
                     upper = upper)
(garchfit <- fit$par)

# ---------------
# -- SANDWICH: --
# ---------------
J <- f$he(fit$par)
f <- MakeADFun(data = list(y=df2, W=W,init = apply(df2,1,var)),
               parameters = parameters, silent = TRUE, map = map,
               ADreport = TRUE, DLL = "STARMAGARCH")

S <- f$gr(fit$par)
Lhat <-matrix(rowSums(matrix(apply(S, 1, function(k) c(-k %*% -t(k))), 
                             nrow = length(fit$par)^2)), ncol = length(fit$par))
SIG <- solve(J) %*% Lhat %*% solve(J)

# -----------------------
# -- Table of results: --
# -----------------------
matcoef2 <- data.frame(Estimates = fit$par,
                      SD = sqrt(diag(SIG)))
matcoef2$Zscore <- matcoef2$Estimates/matcoef2$SD
matcoef2$Pvalue <- 2*stats::pnorm(abs(matcoef2$Zscore), lower.tail=FALSE)
rownames(matcoef2)<-correct.names(names(f$par))
rownames(matcoef2)[1:(nrow(matcoef)-1)]<-rownames(matcoef)[1:(nrow(matcoef)-1)]
rownames(matcoef2)[nrow(matcoef):nrow(matcoef2)]<-c("$\\omega$", 
                                        paste("$\\alpha_{", c(1,1),c(0,1),"}$", sep=""),
                                        paste("$\\beta_{", c(1), c(0),"}$", sep=""))
print(xtable(matcoef2, digits = 3,
                                    caption = "Parameter estimates for a CSTARMA$(2,2)$-GARCH$(1,3)$ model based on the data from Cell 2's mDia1 levels in the first layer."), 
                     sanitize.text.function =function(x)x, hline.after = c(-1,0,nrow(matcoef)-1,nrow(matcoef2)),
      file = "Tex/mdia1_ARMAGARCH_estimates.tex")
# ---------------------
# -- Combined table: --
# ---------------------
matcoef3 <- rbind(matcoef2[1:(nrow(matcoef)-1),], rep(NA,4), matcoef2[nrow(matcoef):nrow(matcoef2),])
rownames(matcoef3)[nrow(matcoef)] <- rownames(matcoef)[nrow(matcoef)]
matcoef3 <- cbind(matcoef3,matcoef3)
matcoef3[1:nrow(matcoef),1:4]<-matcoef
matcoef3[(nrow(matcoef)+1):nrow(matcoef3),1:4]<- NA
tab <- matcoef3[,-c(3:4,7:8)]
tab <- round(tab,3)
tab[,c(1,3)]<-apply(tab[,c(1,3)], 2, function(x) ifelse(x>0, paste("\\phantom{-}",round(x,3),sep=""),round(x,3)))
tab[1:5,]<-apply(tab[1:5,], 2, function(x) paste("\\phic ",x, sep=""))
tab[6:8,]<-apply(tab[6:8,], 2, function(x) paste("\\etac ",x, sep=""))
tab[11:12,3:4]<-apply(tab[11:12,3:4], 2, function(x) paste("\\alphac ",x, sep=""))
tab[13,3:4]<-apply(tab[13,3:4], 2, function(x) paste("\\betac ",x, sep=""))
rownames(tab)<- paste(c(rep("\\phic",5),rep("\\etac",3),"","",rep("\\alphac",2),"\\betac"), rownames(tab))
print(xtable(tab, digits = 3,
             caption = "Parameter estimates for a CSTARMA$(2,2)$-GARCH$(1,3)$ model based on the data from Cell 2's mDia1 levels in the first layer."), 
      sanitize.text.function =function(x)x, hline.after = c(-1,0,nrow(matcoef),nrow(matcoef3)),
      file = "Tex/mdia1_combined_estimates.tex")
cat("CSTARMAGARCH complete.\n")

# ---------------------
# -- Making figures: --
# ---------------------
u <- f$report()$x
mean(u)
sd(u)
h <- f$report()$sigma
yhat <- df2-u
u2 <- u; yhat2 <- yhat; h2 <- h;
plotting_df(df = u2, main = NULL, path = "Figures/10_mdia1_CSTARMAGARCH_residuals.pdf")
plotting_df(df = yhat2[,-1], main = NULL, path = "Figures/10_mdia1_CSTARMAGARCH_fitted_values.pdf")
plotting_df(df = sqrt(h2), main = NULL, path = "Figures/10_mdia1_CSTARMAGARCH_volatility.pdf",
            midpoint = sd(df2), limits = range(sqrt(c(h2))))
dyn.unload(dynlib("Cpp/STARMAGARCH"))



# ------------
# -- Means: --
# ------------
ytilde <- apply(yhat, 2, function(x)x+sfit)
means.df <- data.frame(means = c(sm,rowMeans(ytilde),tm, colMeans(ytilde)), 
                       x = c(rep(1:nrow(df),2),rep(1:ncol(df),2)), 
                       data = c(rep("C", nrow(df)),rep("widehat(C)", nrow(df)),
                                rep("C", ncol(df)),rep("widehat(C)", ncol(df))),
                       type = c(rep("Spatial location", 2*nrow(df)),rep("Time point", 2*ncol(df))))
ggplot(means.df, aes(x, means))+geom_point(aes(color = data), pch = 19)+#geom_line(aes(color=data))+
  facet_wrap(~type, ncol = 2, scales="free_x", 
             strip.position = "bottom")+
  theme_minimal()+labs(colour = NULL)+
  scale_colour_manual(name = NULL, values = c("blue", "red"), 
                     labels = expression(C, widehat(C)))+
  theme(legend.position = "right",
    axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(colour = "black"))+
  xlab(NULL)+ylab("Means")
ggsave("Figures/x_space-time-means_C_vs_Chat.pdf", 
       width = 9, height = 4)

# -------------------
# -- Arc of means: --
# -------------------
tu <- cos(seq(0,pi,length.out = 142))
tv <- sin(seq(0,pi,length.out = 142))
ddddf<-data.frame(tu=tu,tv=tv, sm=sm)

ggplot(ddddf, aes(tu,tv))+geom_line(aes(col=sm), lwd = 5)+theme_minimal()+
  scale_color_gradient2(low=rgb(78,160,183, max = 255), 
                        mid=rgb(242,237,234,max = 255), high=rgb(229,79,71, max = 255), 
                        midpoint = mean(sm), name = "Spatial Mean",
                        guide = guide_colorbar(barheight = 12))+
geom_label(aes(x=0,y=.2),label = "mDia1", size = 24, fill = rgb(78,160,183, max = 255),
           color = rgb(242,237,234,max = 255))+
  #geom_label(aes(x=0,y=.5),label = "mDia1", size = 24, fill = rgb(78,160,183, max = 255),
  #           color = rgb(242,237,234,max = 255))+
  theme(#plot.background = element_rect(fill=rgb(242,237,234,max = 255), color = rgb(242,237,234,max = 255)),
    panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 14), #color = rgb(41,52,81, max=255), 
                                  # face = "bold"),# rgb(78,160,183, max = 255)),
        legend.title = element_text(size = 14, #color = rgb(41,52,81, max = 255), 
                                    face = "bold"))
ggsave("Figures/8_mdia1_cell_illustration.pdf",
       width = 9, height =4)
cat("mDia1 complete.\n")