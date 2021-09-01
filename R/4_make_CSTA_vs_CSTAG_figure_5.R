seed <- 578656187
set.seed(seed)

# -- Parameters: --
true.par <- c(5.032,.354,-.153,.122,-.051,.513,.142,.171)

mu   <- true.par[1]
phi1 <- true.par[2]
phi2 <- true.par[3]
th1  <- true.par[4]
th2  <- true.par[5]
omg  <- true.par[6]
a    <- true.par[7]
b    <- true.par[8]
burnin <- 980
pdim   <- c(30, 1180)

x  <- x2 <- u <- h <- array(0, dim=pdim)
z  <- array(rnorm(prod(pdim)), dim=pdim)
z2 <- sqrt(omg / (1 - sum(a, b) * 3)) * z

for(t in 1:(pdim[2]-1)){
  for(k in 1:pdim[1]){
    h[k,t+1]<- sqrt(omg + a * (u[(k-2)%%pdim[1] + 1, t]^2 + u[k,t]^2 + u[k%%pdim[1] + 1, t]^2)+
                      b * (h[(k-2)%%pdim[1] + 1, t]^2 + h[k,t]^2 + h[k%%pdim[1] + 1, t]^2))
    u[k,t+1]<- h[k,t+1] * z[k,t+1]
  }
}

for(t in 1:(pdim[2]-1)){
  for(k in 1:pdim[1]){
    x[k,t+1]  <- phi1 * x[k,t] + phi2 * (x[(k-2)%%pdim[1] + 1, t] + x[k%%pdim[1] + 1, t])+
      th1 * u[k,t] + th2 * (u[(k-2)%%pdim[1] + 1, t] + u[k%%pdim[1] + 1, t]) + u[k, t+1]
    x2[k,t+1] <- phi1 * x2[k,t] + phi2 * (x2[(k-2)%%pdim[1] + 1, t] + x2[k%%pdim[1] + 1, t])+
      th1 * z2[k,t] + th2 * (z2[(k-2)%%pdim[1] + 1, t] + z2[k%%pdim[1] + 1, t]) + z2[k,t+1]
  }
}

# -- Removing burn-in --
x <- x[,-(1:burnin)]; x2 <- x2[,-(1:burnin)]; h <- h[,-(1:burnin)]; u <- u[,-(1:burnin)];
z <- z[,-(1:burnin)]; z2 <- z2[,-(1:burnin)];

x  <-  x + mu
x2 <- x2 + mu

# ---------------------
# -- Making figure:  --
# ---------------------
library(ggplot2)
library(reshape2)
xdata <- melt(x - mu); x2data <- melt(x2 - mu); hdata <- melt(h); udata <- melt(u); zdata<-melt(z)
names(xdata) <- names(x2data) <- names(hdata) <- names(udata) <- names(zdata) <- c("Space", "Time", "value")
xdata$type   <- "CSTAG"; hdata$type <- "Volatility";
udata$type   <- "GARCH"; zdata$type <- "Innovations";
x2data$type  <- "CSTA"; 

data<-rbind(xdata,x2data)
data$Space <- data$Space - 1
data$Time  <- data$Time  - 1
ggplot(data, aes(Time,Space))+geom_raster(aes(fill=value), interpolate = FALSE) +
  facet_wrap(~ type, ncol = 1, dir = "v") +
  scale_fill_gradient2(low="darkblue", mid = "grey99", high="darkred",
                        na.value = "darkgreen", limits = c(-1,1) * max(abs(data$value))) +
                        scale_x_discrete(expand = c(0,0), limits = c(0,seq(50,200,50))) +
  scale_y_discrete(expand = c(0,0), limits = c(0,seq(5,25,5)))+
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), 
        plot.background = element_blank(),
        legend.key.height = unit(20,"mm"), 
        legend.text = element_text(size = 11),
        legend.title = element_blank(), 
        axis.text = element_text(size = 11),
        strip.text = element_text(face = "bold", 
                                  colour = "black", size = 11), 
        strip.background = element_rect(fill = "white", colour = NULL)
  )
ggsave("Figures/5_CSTA_vs_CSTAG.pdf", width = 9, height = 5)
