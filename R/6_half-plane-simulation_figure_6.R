library(ggplot2)
#library(ggpubr)
uibblue <- rgb(0.31, 0.63, 0.72)
uibred <- rgb(0.90, 0.31, 0.28)
uibgreen <- rgb(0.56, 0.66, 0.41)
uibyellow<-rgb(0.87, 0.74, 0.30)
uibwhite <- rgb(0.949, 0.929, 0.918)
greyly <- rgb(0.98, 0.98, 0.98)


data <- data.frame(x1 = c(0,20,150), 
                   x2 = c(500, 500,500), 
                   y1 = c(0,10,50), 
                   y2 = c(160, 150,110),
                   col =factor(c("Initial values","Burnin", "Sample"))
)

ggplot(data=data, aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2))+geom_rect(aes(fill = col))+
  scale_x_discrete(expand = c(0.01,0))+
  scale_y_discrete(expand = c(0.01,0))+
  scale_fill_manual(values = c(uibblue,uibgreen,uibyellow), name = "")+
  xlab("Time")+ylab("Space")+
  theme_minimal()+
  theme(axis.line = element_line(),#arrow = arrow(angle=8, type = "closed")),
        legend.position = "none",
        axis.title = element_text(size =14))
ggsave("Figures/6_halfplane-simulation8x4.png", width = 8, height = 4)
cat("This is just a shell figure. The details has been added in powerpoints.")