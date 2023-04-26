#' Hydraulic cost function
#' @description Calculates the normalized hydraulic cost as described in Sperry
#'     et al. 2016
#'
#' @param P Vector of equally spaced water potentials ranging from Ps to Pcrit,
#'     -MPa
#' @param b Weibull scale parameter
#' @param c Weibull shape parameter
#' @param kmax_25 Max plant conductance at 25 deg C, mmol s-1 m-2 MPa-1
#' @param T_air Air temperature, deg C
#' @param ratiocrit Percentage of maximum conductivity at which hydraulic damage
#'     is considered irreversible
#'
#' @return Hydraulic cost, unitless
#' @export
#'
#' @examples
#' # Calculate leaf VPD along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#'
#' # Calculate hydraulic cost
#' hydraulic_cost(P, b, c)
hydraulic_cost = function(P,
                          b = -2.5,
                          c = 2,
                          kmax_25 = 4,
                          T_air = 25,
                          ratiocrit = 0.05
)
{
  kmax = calc_kmax(kmax_25, T_air)
  kcrit = ratiocrit*kmax
  k = kmax * vulnerability_curve(P, b, c)
  kmax_i = max(k)
  cost = (kmax_i - k) / (kmax_i - kcrit)
  return(cost)
}




#' Thermal cost function
#' @description Calculates the normalized thermal cost based on F0-T curve
#'
#' @param P Vector of equally spaced water potentials ranging from Ps to Pcrit, -MPa
#' @param b Weibull scale parameter
#' @param c Weibull shape parameter
#' @param kmax_25 Max plant conductance at 25 deg C, mmol s-1 m-2 MPa-1
#' @param T_air Air temperature, deg C
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param RH Relative humidity, in \%
#' @param Patm Atmospheric pressure, kPa
#' @param u Wind speed above the leaf boundary layer, m s-1
#' @param leaf_width Leaf width, m
#' @param T50 Leaf temperature midway between Tcrit and Tmax, deg C
#' @param F0_max Maximum F0 (achieved at Tmax)
#' @param F0_min Minimum (i.e. undamaged) F0
#' @param r Scaling factor for F0-T curve
#'
#' @return Normalized thermal cost, unitless
#' @export
#'
#' @examples
#' # Calculate leaf VPD along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#'
#' thermal_cost(P, b, c)
thermal_cost = function(P,
                        b = -2.5,
                        c = 2,
                        kmax_25 = 4,
                        T_air = 25,
                        PPFD = 1000,
                        RH = 90,
                        Patm = 101.325,
                        u = 2,
                        leaf_width = 0.01,
                        T50  = 51,
                        F0_max = 1000,
                        F0_min = 500,
                        r = 4
                        )
{
  E = trans_from_vc(P, kmax_25, T_air, b, c)
  T_leaf = calc_Tleaf(T_air, PPFD, RH, E, u, Patm, leaf_width)

  F0 = F0_func(T_leaf, T50, F0_max, F0_min, r)

  cost = (F0 - F0_min) / (F0_max - F0_min)
  return(cost)
}





#' Carbon gain
#' @description Calculates the normalized carbon gain as described in Sperry et
#'     al. 2016
#'
#' @param P Vector of equally spaced water potentials ranging from Ps to Pcrit,
#'     -MPa
#' @param b Weibull scale parameter
#' @param c Weibull shape parameter
#' @param kmax_25 Max plant conductance at 25 deg C, mmol s-1 m-2 MPa-1
#' @param T_air Air temperature, deg C
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param Patm Atmospheric pressure, kPa
#' @param u Wind speed above the leaf boundary layer, m s-1
#' @param leaf_width Leaf width, m
#' @param RH Relative humidity, in \%
#'
#' @return Normalized carbon gain
#' @export
#'
#' @examples
#' # Calculate leaf VPD along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#'
#' C_gain(P, b, c)
C_gain = function(P,
                  b = -2.5,
                  c = 2,
                  kmax_25 = 4,
                  T_air = 25,
                  PPFD = 1000,
                  Patm = 101.325,
                  u = 2,
                  leaf_width = 0.01,
                  RH = 60
                  )
{
  E = trans_from_vc(P, kmax_25, T_air, b, c)
  A = calc_A(T_air, PPFD, Patm, E, u, leaf_width, RH)

  Amax = max(A)

  gain = A/Amax
  return(gain)
}
