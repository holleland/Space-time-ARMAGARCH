# --
devtools::install_github("holleland/starmagarch")
library(starmagarch)
library(TMB)

# Compile special versions of the tmb source files: 
sapply(paste("Cpp/",list.files("Cpp/", pattern = "*.cpp"), sep = ""), compile)
rm(list=ls())

source("R/1_offshorewindspeed.R") 


# -- Making figures: --
source("R/3_make_correlation_figure_4.R")
source("R/4_make_CSTA_vs_CSTAG_figure_5.R")
source("R/5_half-plane-simulation_figure_6.R")

# -- run bootstrap --
source("R/6_run_bootstrap_1d.R") # ~ 4 hrs
source("R/7_run_bootstrap_2d.R") # ~ 18 hrs


# Check Chen's condition for moments on X. OK if <1; 
source("R/9_chens_condition.R")
