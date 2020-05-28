# --
library(devtools)
install_github("holleland/starmagarch")
library(starmagarch)
# RUN
setwd("../")
library(TMB)
sapply(paste("Cpp/",list.files("Cpp/", pattern = "*.cpp"), sep = ""), compile)
rm(list=ls())
tstart <- Sys.time() 
source("R/1_mDia1.R") # Around ~ 25 mins due to Sandwich estimators
tslutt <- Sys.time()     
tslutt-tstart

rm(list=ls())
tstart <- Sys.time()
source("R/2_edgevelocity.R") # Around 25 mins due to Sandwich estimators
tslutt <- Sys.time()
tslutt-tstart



# -- Making figures: --
source("R/4_make_correlation_figure_4.R")
source("R/5_make_CSTA_vs_CSTAG_figure_5.R")
source("R/6_half-plane-simulation_figure_6.R")

# -- run bootstrap --
source("R/7_run_bootstrap_1d.R") # ~ 4 hrs
source("R/8_run_bootstrap_2d.R") # ~ 18 hrs

# Run SSTA example: 
source("R/9_SSTA_example.R") # ~ be patient if you want to run this. The code is not optimized.
