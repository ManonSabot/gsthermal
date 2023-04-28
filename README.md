
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsthermal

<!-- badges: start -->
<!-- badges: end -->

## Optimize stomatal conductance while accounting for thermal costs

The package gsthermal models stomatal conductance by maximizing the
instantaneous difference between photosynthetic carbon gain and losses
associated with hydraulic failure and thermal damage. This optimization
model is based on the profit maximization approach of Sperry et
al. (2016) but is distinguished by its inclusion of the cost of thermal
damage, which the Sperry et al. model does not account for.

## Installation

You can install the development version of gsthermal from
[GitHub](https://github.com/CamilleSicangco/gsthermal) with:

``` r
# install.packages("devtools")
devtools::install_github("CamilleSicangco/gsthermal")
```

## Contents

- Hydraulics
  - The function `fit_Weibull` solves for the Weibull parameters for a
    plant’s vulnerability curve.
  - From the Weibull parameters, `calc_Pcrit` solves for the critical
    leaf water potential ($\psi_{leaf}$).
  - `Ps_to_Pcrit` generates a vector of $\psi_{leaf}$ values ranging
    from $\psi_{soil}$ to $\psi_{crit}$.
  - `calc_kmax` calculates the whole-plant conductance at a given
    temperature.
  - `vulnerability_curve` calculates the vulnerability curve this
    vector, and `trans_from_vc` calculates the associated transpiration
    supply function.
- Leaf physiology
  - The functions `calc_Tleaf`, `calc_Dleaf`, `calc_gw`, and `calc_A`
    solve for leaf temperature, leaf VPD, stomatal conductance to water
    vapor, and photosynthesis, respectively.
- Atmospheric variables
  - The functions `calc_Rnet` and `calc_VPD` calculate net radiation and
    atmospheric VPD, respectively.
- Thermal damage
  - `F0_func` calculates the minimum fluorescence F<sub>0</sub> versus
    temperature curve.
- Cost and gain functions
  - `hydraulic_cost` and `thermal_cost` calculate the normalized
    hydraulic and thermal costs, respectively.
  - `C_gain` calculates the normalized carbon gain from photosynthesis.
