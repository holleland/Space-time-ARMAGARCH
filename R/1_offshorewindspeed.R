# -------------------------------------------------
# ---  Script for running wind-speed example:  ----
# -------------------------------------------------
# To install starmagarch package: 
# devtools::install_github("holleland/starmagarch")
library(starmagarch)
library(fpp3)
library(tidyverse)
# Working directory: 
setwd("C:/Users/s15052/Documents/Space-time-ARMAGARCH/")

# Load data and reorder the locations: 
nve <-readRDS("data/offshore_windspeeds_NVE_areas.rds") %>% 
  mutate(
    locID = case_when(
      locID == 12 ~ 13,
      locID == 13 ~ 12,
      locID == 18 ~ 19,
      locID == 19 ~ 18,
      TRUE ~ locID
    ))
nve_points <- nve %>% select(locID, lon, lat) %>% distinct()
# ------------------
# Make MAP figure: 
# ------------------
library(sf)
proj <-  "+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +no_defs +R=6.371e+06"
land <- rnaturalearth::ne_countries(continent = "Europe", 
                                    returnclass = "sf", scale = 10) %>%
  st_transform(proj) %>% 
  select(geometry)
locsSF <- st_as_sf(nve_points, coords = c("lon","lat"), crs= 4326)  %>% st_transform(proj)
locsSF$lon <- st_coordinates(locsSF)[,"X"]
locsSF$lat <- st_coordinates(locsSF)[,"Y"]
locations <- readxl::read_excel("data/NVE_20_locations.xlsx") %>% 
  mutate(name = str_replace(name, "\u00f8","o")) %>% st_as_sf(coords = c("lon","lat"), crs= 4326)  %>% st_transform(proj)
locations$lon <- st_coordinates(locations)[,"X"]
locations$lat <- st_coordinates(locations)[,"Y"]

land %>%
  ggplot() + 
  geom_sf(pch = 15) +
  geom_path(data = locsSF %>% arrange(locID),
            aes(x=lon,y=lat), lty = 2, color = "grey")+
  geom_path(data = locsSF %>%filter(locID%in%c(1,20)),
            aes(x=lon,y=lat), lty = 2, color = "red")+
  geom_polygon(data=locations, aes(x=lon,y=lat, group = name, fill = name), alpha = .2)+
  coord_sf(expand = FALSE, xlim = range(locsSF$lon)+c(-1e5,1e5), ylim = range(locsSF$lat)+c(-1e5,1e5)) + 
  geom_text(data = locsSF,
            aes(x=lon,y=lat, label = locID), vjust=0.5,hjust=.5, angle = 90)+
  
  theme_bw()+
  scale_fill_discrete()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", color = "transparent"),
        strip.text.y.right = element_text(angle = 90),
        strip.text.x = element_text(angle = 90),
        text = element_text(size = 11),
        plot.caption.position = "plot",
        panel.border = element_blank(),
        legend.title = element_text(angle = 90, vjust =.9, size = 18),
        legend.text = element_text(angle = 90, hjust = 0.5),
        plot.background = element_rect(fill = "white", color = "white"),
        legend.background = element_rect(fill = "white", color = "white"))#+
ggsave(file = "figures/map_offshore_areas.pdf", height = 10,width = 5)

# -------------------
# Set up time series: 
origin <- lubridate::as_datetime("1970-01-01", tz = "UTC")
nve$time <- as.integer(nve$time)
nve$datetime <- origin + lubridate::hours(nve$time)
df <- nve %>% 
  select(locID, value) %>% 
  group_by(locID) %>% 
  mutate(time = 1:n()) %>% 
  mutate(datetime = as_datetime("1995-12-31 23:00:00", tz= "UTC")+lubridate::hours(time)) %>% 
  filter(yearmonth(datetime) %in% c("2018 Nov", "2018 Dec",
                                    "2019 Jan", "2019 Feb"))%>% 
  mutate(locID= factor(locID, levels = 1:20))

# Convert to matrix form: 
df_wide <- df %>% pivot_wider(names_from = locID, values_from = value)
y <- df_wide[, as.character(c(1:11,13,12,14:17,19,18,20))] %>%  as.matrix() %>% t()

# Calculate 3-hour averages: 
y <- t(apply(y, 1, function(x){
  colMeans(matrix(x, nrow = 3, byrow=F))}) )

# train-test split: 
train <- which((df_wide$datetime<as.Date("2019-02-01"))[seq(3,nrow(df_wide),3)])
y_train <- y[,train]
y_test <- y[,-train]

# # Set up neighbourhood matrix correcting weights:
W_c <- create.neighbourhood.array(m = 20, sp = 2, type = "rook", torus =TRUE)
diag(W_c[,,2]) <- 0
W_c[,,2] <- W_c[,,2]*2 


# ---------------------------
# ---- SET-UP MODELS: -------
# ---------------------------
# Initial parameters. 
phi <- matrix(0, ncol = 1, nrow =2)
theta <- matrix(0.0, ncol =1, nrow =2)
alpha <- matrix(0, ncol = 1, nrow = 2)
beta <- matrix(0, ncol = 1, nrow=2)
omega <- 1
initial.parameters <- list(
  mu    =  mean(y_train),
  phi   = phi,
  theta = theta,
  omega = omega,
  alpha = alpha,
  beta  = beta
)

# CSTAG model: 
f_c <- CreateLikelihood(y_train, W=W_c,
                        init = apply(y_train,1,var), 
                        parameters=initial.parameters)
fit_c <- fitSTARMAGARCH(f_c, data = y_train, print = TRUE)
fit_c %>% coef() %>% tail(4) %>% sum()

# Check stationarity condtions:
A <- fit_c$coef[2]*W_c[,,1]+fit_c$coef[3]*W_c[,,2]
all(eigen(A)$values < 1)
B <- fit_c$coef[4]*W_c[,,1]+fit_c$coef[5]*W_c[,,2]
eigs <- fft(as.numeric(B[1,]), inverse = FALSE)  # complex eigenvalues
all(Mod(eigs) <1)

B <- (fit_c$coef[7]+fit_c$coef[9])*W_c[,,1]+(fit_c$coef[8]+fit_c$coef[10])*W_c[,,2]
B <- (.6)*W_c[,,1]+(.39)*W_c[,,2]
eigs <- fft(as.numeric(B[1,]), inverse = FALSE)  # complex eigenvalues
all(Mod(eigs)<1)

# Non-spatial model: 
map_n <- list(
  phi = factor(c(1, NA)),
  theta = factor(c(2,NA)),
  alpha = factor(c(3,NA)),
  beta = factor(c(4,NA))
)
f_n <- CreateLikelihood(y_train, W=W_c,
                        init = apply(y_train,1,var), 
                        parameters=initial.parameters,
                        map = map_n)
fit_n <- fitSTARMAGARCH(f_n, data = y_train, print = TRUE)

# Naive model: 
map_naive <- list(
  mu = factor(NA), 
  phi = factor(c(NA, NA)),
  theta = factor(c(NA,NA)),
  alpha = factor(c(NA,NA)),
  beta = factor(c(NA,NA))
)
initial.parameters$phi[1,1] <-1 
f_naive <- CreateLikelihood(y_train, W=W_c,
                            init = apply(y_train,1,var), 
                            parameters=initial.parameters,
                            map = map_naive)
fit_naive <- fitSTARMAGARCH(f_naive, data = y_train, print = TRUE)

# Pure ARMA model:
map_arma <- list(
  alpha = factor(c(NA,NA)),
  beta = factor(c(NA,NA))
)
f_arma <- CreateLikelihood(y_train, W=W_c,
                            init = apply(y_train,1,var), 
                            parameters=initial.parameters,
                            map = map_arma)
fit_arma <- fitSTARMAGARCH(f_arma, data = y_train, print = TRUE)



# Function for predicting on new dataset:
predictions <- function(fit, W, new_y, model.name=" "){
  init <- sigma(fit)[,ncol(sigma(fit))]
  pars <- coef(fit)
  if(model.name =="N"){
    names(pars) <- c("mu","phi1","theta1","omega","alpha1","beta1")
    pars["phi2"]<-pars["theta2"]<-pars["alpha2"]<-pars["beta2"] <-0
  }
  if(length(pars)==1){
    pars["phi1"] <- 1
    pars["theta1"] <- pars["alpha1"]<-pars["beta1"] <- pars["mu"] <-0
    pars["phi2"]<-pars["theta2"]<-pars["alpha2"]<-pars["beta2"] <-0
  }
  if(model.name == "arma"){
    pars["alpha1"]<-pars["alpha2"]<-pars["beta1"]<-pars["beta2"]<-0
  }
  x<- z <- new_sig<- new_yhat <- new_y*0
  new_sig[,1] <- sqrt(pars["omega"] + 
    pars["alpha1"] * W[,,1] %*% fit$garch[,ncol(y_test)]^2 + 
    pars["alpha2"] *W[,,2] %*% fit$garch[,ncol(y_test)]^2 + 
    pars["beta1"] * W[,,1]%*%init^2 +  
    pars["beta2"] * W[,,2]%*%init^2)
  new_yhat[,1] <- pars["mu"] + 
    pars["phi1"]*W[,,1]%*% (y_train[,ncol(y_train)]-pars["mu"])+
    pars["phi2"]*W[,,2]%*% (y_train[,ncol(y_train)]-pars["mu"]) +
    pars["theta1"] * W[,,1]%*% fit$garch[,ncol(y_train)]+ 
    pars["theta2"] * W[,,2] %*% fit$garch[,ncol(y_train)]  
  x[,1]<- new_y[,1]-new_yhat[,1]
  z[,1] <- x[,1]/new_sig[,1]
  for(t in 2:ncol(new_y)){
    new_sig[,t] <- sqrt(pars["omega"] + 
      pars["alpha1"] * W[,,1] %*% x[,t-1]^2 + 
      pars["alpha2"] *W[,,2] %*% x[,t-1]^2 + 
      pars["beta1"] * W[,,1]%*%new_sig[,t-1]^2 +  
      pars["beta2"] * W[,,2]%*%new_sig[,t-1]^2)
    new_yhat[,t] <- pars["mu"] + 
      pars["phi1"]*W[,,1]%*% (new_y[,t-1]-pars["mu"])+
      pars["phi2"]*W[,,2]%*% (new_y[,t-1]-pars["mu"]) +
      pars["theta1"] * W[,,1]%*% x[,t-1]+ 
      pars["theta2"] * W[,,2] %*% x[,t-1]  
    x[,t]<- new_y[,t]-new_yhat[,t]
    z[,t] <- x[,t]/new_sig[,t]
  }
  plot(c(new_y), c(new_yhat))
  
  pred_tbl <- left_join(
    reshape2::melt(x) %>% rename(x=value),
    reshape2::melt(new_y) %>% rename(y=value)) %>% 
    left_join(reshape2::melt(new_yhat) %>% rename(yhat=value)) %>% 
    left_join(reshape2::melt(new_sig) %>% rename(sigma=value)) %>% 
    left_join(reshape2::melt(z) %>% rename(z=value)) %>% 
    rename(u = Var1, t = Var2)
  pred_tbl %>% 
    pivot_longer(cols = -(1:2)) %>% 
    filter(name %in% c("y","yhat")) %>% 
  ggplot(aes(x=t, y =u, fill = value)) + geom_tile()+
    facet_wrap(~name, ncol = 1,strip.position = "left")+
    scale_x_continuous(expand =c(0,0), name = "Space")+
    scale_y_continuous(expand =c(0,0), name = "Time")+
    theme_bw()+
    theme(strip.placement = "outside",
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", color = "transparent"))+
    scale_fill_viridis_c(name = "Wind speed (demeaned)")
  pred_tbl %>% 
    pivot_longer(cols = -(1:2)) %>% 
    filter(!(name %in% c("y","yhat"))) %>% 
    mutate(name = factor(name, levels = c("x","sigma","z"))) %>% 
    ggplot(aes(x=t, y =u, fill = value)) + geom_tile()+
    facet_wrap(~name, ncol = 1,strip.position = "left")+
    scale_x_continuous(expand =c(0,0), name = "Space")+
    scale_y_continuous(expand =c(0,0), name = "Time")+
    theme_bw()+
    theme(strip.placement = "outside",
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", color = "transparent"))+
    scale_fill_viridis_c(name = "Wind speed (demeaned)")
  
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggnewscale)
  # Split by variable, then plot separately in one figure
  p1 <- pred_tbl %>%
    pivot_longer(cols = -(1:2)) %>%
    filter(name == "y") %>%
    ggplot(aes(x = t, y = u, fill = value)) +
    geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                  midpoint = 0.0, name ="Y")+
    scale_y_continuous(expand =c(0,0), name = "Space")+
    scale_x_continuous(expand =c(0,0), name = "Time")+
    theme_bw()+
    theme(axis.title.y = element_blank())+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  p2 <- pred_tbl %>%
    pivot_longer(cols = -(1:2)) %>%
    filter(name == "yhat") %>%
    ggplot(aes(x = t, y = u, fill = value)) +
    geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                  midpoint = 0.0, name ="Yhat")+
    scale_y_continuous(expand =c(0,0), name = "Space")+
    scale_x_continuous(expand =c(0,0), name = "Time")+
    theme_bw()+
    theme(axis.title.y = element_blank())+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  
  p3 <- pred_tbl %>%
    pivot_longer(cols = -(1:2)) %>%
    filter(name == "x") %>%
    ggplot(aes(x = t, y = u, fill = value)) +
    geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                  midpoint = 0.0, name ="X")+
    scale_y_continuous(expand =c(0,0), name = "Space")+
    scale_x_continuous(expand =c(0,0), name = "Time")+
    theme_bw()+
    
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  
  p4 <- pred_tbl %>%
    pivot_longer(cols = -(1:2)) %>%
    filter(name == "sigma") %>%
    ggplot(aes(x = t, y = u, fill = value)) +
    geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                  midpoint = 0.0, name ="Sigma")+
    scale_y_continuous(expand =c(0,0), name = "Space")+
    scale_x_continuous(expand =c(0,0), name = "Time")+
    theme_bw()
  p4_n <- p4 +
    theme(axis.title.y = element_blank())+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  p5 <- pred_tbl %>%
    pivot_longer(cols = -(1:2)) %>%
    filter(name == "z") %>%
    ggplot(aes(x = t, y = u, fill = value)) +
    geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                  midpoint = 0.0, name ="Z")+
    scale_y_continuous(expand =c(0,0), name = "Space")+
    scale_x_continuous(expand =c(0,0), name = "Time")+
    theme_bw() + theme(axis.title.y = element_blank())
  
  library(patchwork)
  print(p1 / p2 / p3 / p4_n / p5)
  ggsave(paste0("figures/windspeed_", model.name, "_model_test_set.pdf"),
         width = 8, height = 10)
  print(p4)
  # Figure 9 (for CSTAG model):
  ggsave(paste0("figures/windspeed_sigma_", model.name, "_model_test_set.pdf"),
         width = 8, height = 4)
  
  return(pred_tbl)
}
pred_c<- predictions(fit_c, W_c, y_test, model.name = "C")
pred_n<- predictions(fit_n, W_c, y_test, model.name = "N")
pred_naive<- predictions(fit_naive, W_c, y_test, model.name = "naive")
pred_arma<- predictions(fit_arma, W_c, y_test, model.name = "arma")

# ------------------------
# ---- MAKING TABLE ------
# ------------------------
rmse <- c(sqrt(mean((pred_c$x)^2)),
          sqrt(mean((pred_arma$x)^2)),
          sqrt(mean((pred_n$x)^2)),
          sqrt(mean((pred_naive$x)^2))
  )
rmse <- round(rmse, 3)
aic <-c( AIC(fit_c), AIC(fit_arma), AIC(fit_n),AIC(fit_naive))
bic <-c( BIC(fit_c), BIC(fit_arma),BIC(fit_n),BIC(fit_naive))

aic <- round(aic-aic[1])

preds <- bind_rows(pred_c %>% mutate(model = "Circular"),
                   pred_n %>% mutate(model = "Non-spatial"),
                   pred_naive %>% mutate(model = "Naive"),
                   pred_arma %>% mutate(model = "ARMA"))
coverage <- preds %>%  group_by(model) %>% 
  mutate(covered = between(y, yhat-1.96*sigma, yhat+1.96*sigma)) %>% 
  summarize(mean = mean(covered))
cover <- paste0(round(100*coverage[c(2,1,4,3),2] %>% pull(mean),1),"\\%")

coefs_csta <- c(paste0(round(coef(fit_arma),3),"(",round(100*fit_arma$matcoef$SD,2),")"),rep("",4))
coefs_cstag <- paste0(round(coef(fit_c),3),"(",round(100*fit_c$matcoef$SD,2),")")
coefs_ns <- c(round(coef(fit_n)[1:2],3),"",
              round(coef(fit_n)[3],3), "",
              round(coef(fit_n)[4:5],3), "",
              round(coef(fit_n)[6],3), "")
se_ns <- c(round(100*fit_n$matcoef$SD[1:2],2),"",
           round(100*fit_n$matcoef$SD[3],2), "",
           round(100*fit_n$matcoef$SD[4:5],2),"",
           round(100*fit_n$matcoef$SD[6],2), "")
coef_ns <- paste0(coefs_ns,"(",se_ns,")")
coef_ns[c(3,5,8,10)] <- ""
coef_naive <- rep("",10)
coef_naive[6]<- paste0(round(coef(fit_naive),3),"(", 100*round(fit_naive$matcoef$SD[1],2),")")
tbl <- tibble(
  "CSTAG" = c(aic[1], rmse[1], cover[1], coefs_cstag),
  "CSTA"= c(aic[2],rmse[2], cover[2],coefs_csta),
  "Non-Spatial"  = c(aic[3], rmse[3], cover[3],
                             coef_ns),
  "Naive" = c(aic[4], rmse[4], cover[4],coef_naive)) %>% 
  mutate(rownames = c( "AIC", "RMSE","Coverage", paste0("$\\", names(coef(fit_c)), "$"))) %>% 
  column_to_rownames(var = "rownames")

# PRINT TABLE: 
library(kableExtra)
knitr::kable(tbl, format = "latex", digits = 4, escape = FALSE, booktabs = TRUE) %>%
      kable_styling(latex_options = c("striped", "hold_position")) %>%
      add_header_above(c(" " = 1, "Models" = 4))

# Not used graphic: 
pred_c %>% 
  ggplot(aes(x=t)) + 
  geom_line(aes(y=yhat), color= "black", lwd = .3, alpha= .9) +
  geom_ribbon(aes(ymin = yhat - 1.96*sigma, ymax = yhat + 1.96*sigma), fill = "blue", alpha = .5)+
  geom_point(aes(y=y), color = "red", size = .002) +
  facet_wrap(~u, scales = "fixed")+
  scale_x_continuous("Time", expand = c(0,0))+
  scale_y_continuous("Log wind speed") +
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
ggsave("figures/one_step_ahead_prediction_testset.pdf", width = 12, height = 8)

# Figure 8: 
preds %>%
  filter(
    u==4,
          model %in% c("Circular","ARMA")
         ) %>% 
  mutate(model = ifelse(model == "Circular", "CSTAG", "CSTA")) %>% 
  ggplot(aes(x=t))+
#  geom_line(aes(y =  1.96*sigma)) +
#  geom_line(aes(y = -1.96*sigma))+
  geom_ribbon(aes(ymin = -1.96*sigma,ymax = 1.96*sigma), fill = "skyblue", alpha = .4)+
  geom_point(aes(y=y-yhat), color = "black")+
 # facet_grid(model~u) + theme_bw()+
  facet_wrap(~model, ncol =1,strip.position = "left") + theme_bw()+
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "transparent"))+
  scale_x_continuous(expand = c(0,0), name = "Time")+
  scale_y_continuous("Prediction error on location 4")
ggsave("Figures/difference_between_circular_A-vs-AG.pdf", 
       width = 8, height = 6)  

# Check p-value statement
all(c(fit_c$matcoef$Pvalue,
    fit_n$matcoef$Pvalue,
    fit_naive$matcoef$Pvalue,
    fit_arma$matcoef$Pvalue)<1e-2)

# Per-location RMSE and coverage:
location_stats <- preds %>%
  group_by(model, u) %>%
  mutate(covered = between(y, yhat - 1.96*sigma, yhat + 1.96*sigma)) %>%
  summarize(
    RMSE = round(sqrt(mean(x^2)), 3),
    Coverage = round(100 * mean(covered), 1),
    .groups = "drop"
  ) %>%
  filter(model %in% c("Circular", "ARMA")) %>%
  mutate(model = ifelse(model == "Circular", "CSTAG", "CSTA"))

# Wide format table for paper:
location_table <- location_stats %>%
  pivot_wider(
    names_from = model,
    values_from = c(RMSE, Coverage)
  ) %>%
  rename(Location = u) %>%
  select(Location, RMSE_CSTAG, RMSE_CSTA, Coverage_CSTAG, Coverage_CSTA)

# Print as latex table:
knitr::kable(location_table, format = "latex", booktabs = TRUE,
             col.names = c("Location", "CSTAG", "CSTA", "CSTAG", "CSTA"),
             escape = FALSE) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  add_header_above(c(" " = 1, "RMSE" = 2, "Coverage (\\%)" = 2))
