library(xtable)
library(TMB)
library(reshape2)
library(ggplot2)
library(spdep)
library(starma)
source("R/Functions.R")
# ----------------
# -- Load data: --
# ----------------
df <- read_cell("Data/Chan0_1L_imActmap.txt")
# -- Removing NAs: 
df <- df[,-1]

# ------------------- #
# --- Detrending ---- #
# ------------------- #
df.melt <- melt(df)
names(df.melt)[1:2] <- c("u","t")

# ------------------------------------------- #
# --- Plotting spatial- and temporal means -- #
# ------------------------------------------- #
sm <- rowMeans(df) # spatial
tm <- colMeans(df) # temporal 
means.df <- data.frame(means = c(sm,tm), 
                       x = c(1:nrow(df), (1:ncol(df))*10), 
                       type = c(rep("Space", nrow(df)),rep("Time (sec)", ncol(df))))
ggplot(means.df, aes(x, means))+geom_point(color = "blue")+
  facet_wrap(~type, ncol = 2, scales="free_x", 
             strip.position = "bottom")+
  theme_minimal()+
  theme(axis.title = element_text(size = 14),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(colour = "black"))+
  xlab(NULL)+ylab("Means")+
  geom_hline(aes(yintercept = 0), color = "red", lwd = .6)
ggsave("Figures/9_edge_velocity_space-time-means.pdf", 
       width = 9, height = 4)

# No trends: df2 = df
df2 <- df
plotting_df(df2, main = NULL,midpoint = 0,#limits = c(-1,1)*max(abs(c(df2,df3))),#"(d) Detrended Cell 2 mDia1 levels at Layer 1.", 
            path = "Figures/9_edge_velocity_rawdata.pdf")

# ----------------------
# -- Fitting CSTARMA: --
# ----------------------
compile("Cpp/STARMA.cpp")
dyn.load(dynlib("Cpp/STARMA"))
parameters <- list(mu = 0.0, 
                   phi   = matrix(0, ncol = 3, nrow=2),
                   theta = matrix(0, ncol = 2, nrow = 2),
                   sigma = sd(df))
W <- create.neighbourhood.array(m = c(nrow(df), 1), sp = 5, torus = TRUE, type = "rook")
phimap <- matrix(c(1:6), ncol = 3, nrow = 2)
phimap[c(1,3,4)]<-NA_real_
thetamap <- matrix(6+c(1:4), ncol = 2, nrow = 2)
thetamap[4]<-NA_real_
map <- list(#mu = factor(NA),
            phi=factor(phimap),
            theta = factor(thetamap))
f <- MakeADFun(data = list(y=df2, W=W,init = apply(df2,1,var)),
               parameters = parameters, silent = TRUE, map = map, DLL = "STARMA")
lower <- c(-1e4,rep(-5,length(which(names(f$par)%in% c("phi","theta")))),1e-8)
upper <- c(1e4,rep(5,length(which(names(f$par)%in% c("phi","theta")))),1e5)

fit <- stats::nlminb(f$par,f$fn,f$gr, f$he,
                     lower = lower,
                     upper = upper)
# -------------------------
# -- Sandwich estimator: --
# -------------------------
J <- f$he(fit$par)
#
# With ADreport = TRUE, TMB calculate the gradient of ell(t,u) (instead of L)
#
f <- MakeADFun(data = list(y=df2, W=W,init = apply(df2,1,var)),
               parameters = parameters, silent = TRUE, map = map,
               ADreport = TRUE, DLL = "STARMA")

S <- f$gr(fit$par) # Calculate the gradients
Lhat <-matrix(rowSums(matrix(apply(S, 1, function(k) c(-k %*% -t(k))), 
                              nrow = length(fit$par)^2)), ncol = length(fit$par))
SIG <- solve(J) %*% Lhat %*% solve(J) # Sandwich
# -----------------------------------------
# -- Making table of estimation results: --
# -----------------------------------------
armapar <- fit$par
matcoef <- data.frame(Estimates = fit$par,
                      SD = sqrt(diag(SIG)))
matcoef$Zscore <- matcoef$Estimates/matcoef$SD
matcoef$Pvalue <- 2*stats::pnorm(abs(matcoef$Zscore), lower.tail=FALSE)
rownames(matcoef)<-correct.names(names(f$par))
matcoef

rownames(matcoef)<-c("$\\mu$",paste("$\\phi_{", c(1,3,3),c(1,0,1),"}$", sep = ""),
                    paste("$\\eta_{", c(1,1,2),c(0,1,0),"}$", sep=""), "$\\sigma$")

print(xtable(matcoef, digits = 3,
             caption = "Parameter estimates with subscripts $jv$ where $j$ is temporal- and $v$ is spatial lag,  for a CSTARMA$(2,2)$ model based on the detrended mDia1 levels in the first layer of Cell no. 2."), 
      sanitize.text.function =function(x)x,
      file = "Tex/Cell2_Chan0_1L_ARMA_estimates.tex")
cat("CSTARMA complete.\n")
# ---------------------------
# -- Fitting CSTARMAGARCH: --
# ---------------------------
dyn.unload(dynlib("Cpp/STARMA"))
compile("Cpp/STARMAGARCH.cpp")
dyn.load(dynlib("Cpp/STARMAGARCH"))
parameters <- list(mu = armapar[1],#mean(df2),
                   phi   = matrix(c(0,armapar[2],0,0,armapar[3:4]), ncol = 3, nrow=2),
                   theta = matrix(c(armapar[5:7],0), ncol = 2, nrow = 2),
                   omega = armapar[8]^2,
                   alpha = matrix(0, ncol = 2, nrow=2),
                   beta  = matrix(0, 
                                  ncol=2,nrow=2))#c(,0,.04,rep(0,3),.05,0,0), ncol = 3, nrow=3))
W <- create.neighbourhood.array(m = c(nrow(df), 1), sp = 4, torus = TRUE, type = "rook")

alphamap <-matrix(18+1:4, ncol = 2, nrow =2)
alphamap[4]<-NA_real_
betamap <- matrix(47+1:4, ncol = 2, nrow =2)
betamap[c(2:3)]<-NA_real_
map$alpha <- factor(alphamap)
map$beta <- NULL
map$beta <- factor(betamap)

f <- MakeADFun(data = list(y=df2, W=W,init = rep(sd(df2)^2, nrow(df2))),#rep(armapar[length(armapar)]^2, nrow(df2))),#apply(df2,1,var)),
               parameters = parameters, silent = TRUE, map = map, DLL = "STARMAGARCH")
lower <- c(-1e4,rep(-5,length(which(names(f$par)%in% c("phi","theta")))),1e-8,
           rep(1e-8, length(which(names(f$par)%in% c("alpha","beta")))))
upper <- c(1e4,rep(5,length(which(names(f$par)%in% c("phi","theta")))),1e5,
           rep(1, length(which(names(f$par)%in% c("alpha","beta")))))

fit <- stats::nlminb(f$par,f$fn,f$gr, f$he,
                     lower = lower,
                     upper = upper)
(garchfit <- fit$par)
# -------------------------
# -- Sandwich estimator: --
# -------------------------
J <- f$he(fit$par)
#
# With ADreport = TRUE, TMB calculate the gradient of ell(t,u) (instead of L)
#
f <- MakeADFun(data = list(y=df2, W=W,init = rep(sd(df2)^2, nrow(df2))),#rep(armapar[length(armapar)]^2, nrow(df2))),#apply(df2,1,var)),
               parameters = parameters, silent = TRUE, map = map,
               ADreport = TRUE, DLL = "STARMAGARCH")

S <- f$gr(fit$par)
Lhat <-matrix(rowSums(matrix(apply(S, 1, function(k) c(-k %*% -t(k))), 
                             nrow = length(fit$par)^2)), ncol = length(fit$par))
SIG <- solve(J) %*% Lhat %*% solve(J)

# ------------------------------
# -- Making resulting tables: --
# ------------------------------
matcoef2 <- data.frame(Estimates = fit$par,
                      SD = sqrt(diag(SIG)))
matcoef2$Zscore <- matcoef2$Estimates/matcoef2$SD
matcoef2$Pvalue <- 2*stats::pnorm(abs(matcoef2$Zscore), lower.tail=FALSE)
rownames(matcoef2) <- correct.names(names(f$par))
rownames(matcoef2)[1:(nrow(matcoef)-1)]<-rownames(matcoef)[1:(nrow(matcoef)-1)]
rownames(matcoef2)[nrow(matcoef):nrow(matcoef2)]<-c("$\\omega$", 
                                        paste("$\\alpha_{", c(1,1,2),c(0,1,0),"}$", sep=""),
                                        paste("$\\beta_{", c(1,2), c(0,1),"}$", sep=""))
print(xtable(matcoef2, digits = 3,
                                    caption = "Parameter estimates for a CSTARMA$(2,2)$-GARCH$(1,3)$ model based on the data from Cell 2's mDia1 levels in the first layer."), 
                     sanitize.text.function =function(x)x, hline.after = c(-1,0,nrow(matcoef)-1,nrow(matcoef2)),
      file = "Tex/edge_velocity_ARMAGARCH_estimates.tex")

# ----------------------------
# -- Making combined table: --
# ----------------------------
matcoef3 <- rbind(matcoef2[1:(nrow(matcoef)-1),], rep(NA,4), matcoef2[nrow(matcoef):nrow(matcoef2),])
rownames(matcoef3)[nrow(matcoef)] <- rownames(matcoef)[nrow(matcoef)]
matcoef3 <- cbind(matcoef3,matcoef3)
matcoef3[1:nrow(matcoef),1:4]<-matcoef
matcoef3[(nrow(matcoef)+1):nrow(matcoef3),1:4]<- NA
tab <- matcoef3[,-c(3:4,7:8)]
tab <- round(tab,3)
tab[,c(1,3)]<-apply(tab[,c(1,3)], 2, function(x) ifelse(x>0, paste("\\phantom{-}",round(x,3),sep=""),round(x,3)))
tab[2:4,]<-apply(tab[2:4,], 2, function(x) paste("\\phic ",x, sep=""))
tab[5:7,]<-apply(tab[5:7,], 2, function(x) paste("\\etac ",x, sep=""))
tab[10:12,3:4]<-apply(tab[10:12,3:4], 2, function(x) paste("\\alphac ",x, sep=""))
tab[13:14,3:4]<-apply(tab[13:14,3:4], 2, function(x) paste("\\betac ",x, sep=""))
rownames(tab)[-c(1,8)]<- paste(rep(c("\\phic","\\etac","\\alphac","\\betac"), each = 3)[-12], rownames(tab)[-c(1,8)])
print(xtable(tab, digits = 3,
             caption = "Parameter estimates for a CSTARMA$(2,2)$-GARCH$(1,3)$ model based on the data from Cell 2's mDia1 levels in the first layer."), 
      sanitize.text.function =function(x) x, hline.after = c(-1,0,nrow(matcoef),nrow(matcoef3)),
      file = "Tex/edge_velocity_ARMAGARCH_combined_estimates.tex")
cat("CSTARMAGARCH complete.\n")
# -------------
# -- FIGURES --
# -------------
u <- f$report()$x
h <- f$report()$sigma
yhat <- df2-u
plotting_df(df = u, main = NULL,
            path = "Figures/10_edge_velocity_CSTARMAGARCH_residuals.pdf")
plotting_df(df = yhat[,-1], main = NULL,
            path = "Figures/10_edge_velocity_CSTARMAGARCH_fitted_values.pdf")
plotting_df(df = sqrt(h), main = NULL,
            path = "Figures/10_edge_velocity_CSTARMAGARCH_volatility.pdf",
            midpoint = sd(df2),
            limits =c(range(sqrt(c(h)))))

dyn.unload(dynlib("Cpp/STARMAGARCH"))

# ------------
# -- Means: --
# ------------
ytilde <- yhat
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
ggsave("Figures/x_edge_velocity_space-time-means_C_vs_Chat.pdf", 
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
geom_label(aes(x=0,y=.2),label = "Edge velocity", size = 24, fill = rgb(78,160,183, max = 255),
           color = rgb(242,237,234,max = 255))+
  theme(#plot.background = element_rect(fill=rgb(242,237,234,max = 255), color = rgb(242,237,234,max = 255)),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 14), #color = rgb(41,52,81, max=255), 
                                   #face = "bold"),# rgb(78,160,183, max = 255)),
        legend.title = element_text(size = 14, #color = rgb(41,52,81, max = 255), 
                                    face = "bold"))

ggsave("Figures/8_edge_velocity_cell_illustration.pdf",
       width = 9, height =4)
cat("Edge velocity complete.\n")