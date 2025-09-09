Stationary spatiotemporal ARMAGARCH models and their boundary issues
================

Contributors: Sondre Hølleland (Norwegian School of Economics), Hans
Arnfinn Karlsen (University of Bergen).

This repository contains the necessary code for reproducing results from
the paper **Stationary spatiotemporal ARMAGARCH models and their
boundary issues**.

## Data

Wind speed data from
[NORA3-WP](https://archive.sigma2.no/dataset/482CC467-9E4F-4377-9E05-CB9822938D07)
by

[Solbrekke, I., Sorteberg, A., University of Bergen (2021). Norwegian
hindcast archive’s wind power data set (NORA3-WP) \[Data set\]. NIRD
RDA.](https://doi.org/10.11582/2021.00068).

The portion of the dataset used in the paper is available in the file
**data/offshore_windspeeds_NVE_areas.rds**, under the same licence as
the original data.

## Code

In the R folder, you will find the necessary code for reproducing
figures and tables in the simulations and wind speed forecasting
application. There is also some additional cpp-files in the Cpp folder
used with the Template Model Builder (TMB). We also rely on the
[starmagarch](https://github.com/holleland/starmagarch) package,
especially for the wind speed application.

## Licence

This project is licensed under [CC BY 4.0 LEGAL
CODE](https://creativecommons.org/licenses/by/4.0/legalcode), same as
the NORA3-WP data.
