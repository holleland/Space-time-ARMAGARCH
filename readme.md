Space-Time ARMA-GARCH models with applications
================

*Contributors: Sondre Hølleland<sup>†</sup>, Hans Arnfinn Karlsen.
University of Bergen, Norway.*

*<sup>†</sup> Responsible for the code.*

*Correspondance to: <sondre.holleland@uib.no>*

*The full paper can be found here (link will come).*

This repository contains the necessary code for reproducing results from
the paper *Space-Time ARMA-GARCH models with applications*.

To redo the analysis, run the separate scripts from *R/0\_main.R*.
Notice the approximate times it takes to execute some of the scripts.
This is mostly due to not treating the neighbourhood matrices as sparse
matrices. This will (hopefully) be implemented later.

## Datasets

The cell data for mDia1 (Data/Chan2\_1L\_imActmap.txt) and edge velocity
(Data/Chan0\_1L\_imActmap.txt) in a skin cell is published by courtesy
of [Lee et al (2015)](https://doi.org/10.1016/j.cels.2015.07.001).
Thanks to Jaewon Huh for providing us with the data. The spatially
differenced sea surface temperature anomalies
(Data/SSTA\_spatially\_differenced\_without\_land.RData) used by
[Hølleland and Karlsen(2020)](https://doi.org/10.1111/jtsa.12498) are
derived from an example dataset for the book by Cressie and Wikle (2011)
and the original data can be downloaded
[here](ftp://ftp.wiley.com/public/sci_tech_med/spatio_temporal_data).

## References

  - Hølleland, S., & Karlsen, H. A. (2020). A Stationary Spatio‐Temporal
    GARCH Model. Journal of Time Series Analysis, 41(2), 177-209.
  - Lee, K., Elliott, H. L., Oak, Y., Zee, C. T., Groisman, A., Tytell,
    J. D., & Danuser, G. (2015). Functional hierarchy of redundant actin
    assembly factors revealed by fine-grained registration of intrinsic
    image fluctuations. Cell systems, 1(1), 37-50.
  - Cressie, N., & Wikle, C. K. (2011). Statistics for spatio-temporal
    data. John Wiley & Sons.

## Author’s github account

**Sondre Hølleland** - [holleland](https://github.com/holleland)

## License

This project is licensed under the GNU GPLv3 License - see
[LICENSE.md](LICENSE.md) for details.
