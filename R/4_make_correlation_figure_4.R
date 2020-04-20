library(reshape2)
library(ggplot2)
source("R/Functions.R")
M <- 100
n <- 100

W <- create.neighbourhood.array(m = c(n,1), sp = 4, type = "rook", torus =TRUE)
library("matrixcalc")
PHI <-  W[,,1] * .5 + W[,,2] * .1
C   <- PHI * 0
for(k in 0:100){
  C <- C + matrix.power(PHI, 2*k)
}

rho <- matrix(NA_real_, ncol = 15, nrow = 4)
for(h in 1:15)
  rho[,h] <- (diag(sqrt(1/diag(C))) %*% matrix.power(PHI,h) %*% diag(sqrt(1/diag(C))))[1:4, 1]
rownames(rho) <- 1:4
colnames(rho) <- 1:15
Mrho<-melt(rho)
names(Mrho)<-c("Space","Time","rho")
Mrho <- Mrho[Mrho$Time <10,]
Mrho$Space <- factor(Mrho$Space-1)
lab <- 1:9
lab[seq(1,9,2)] <- ""
library(grid)
ggplot(Mrho,  aes(x=Time, y = rho)) +
  geom_line(aes(col = Space), size = 1.05, lty = 1) +
  scale_x_continuous(expand = c(0.05,0), breaks = seq(1,9,1), labels = lab) +
  scale_y_continuous(expand = c(0,0.005)) +
  scale_color_manual(name = "Spatial lag",values = c("0"=rgb(0.90, 0.31, 0.28),
                                                     "1"=rgb(0.31, 0.63, 0.72),
                                                     "2"=rgb(0.87, 0.74, 0.30), 
                                                     "3"=rgb(0.56, 0.66, 0.41))) +
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  theme_minimal() + ylab("Autocorrelation") + xlab("Temporal lag") +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key=element_rect(fill=NA, color = "white"),
        axis.ticks = element_line(color = "black", size = .9), 
        axis.line = element_line(color = "black", size = .9),
        panel.background = element_rect(fill ="white", linetype = "blank"),
        plot.background = element_rect(fill = "white", linetype = "blank"))
ggsave("Figures/4_correlation_example.pdf", width = 8, height = 3)
