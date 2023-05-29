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
#' @param constant_kmax TRUE if the kmax does not vary with temperature for
#'     simulations; else FALSE
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
                          ratiocrit = 0.05,
                          constant_kmax = FALSE
)
{
  kmax = calc_kmax(kmax_25, T_air, constant_kmax)
  kcrit = ratiocrit*kmax
  k = kmax * vulnerability_curve(P, b, c)
  kmax_i = max(k)
  cost = (kmax_i - k) / (kmax_i - kcrit)
  return(cost)
}




#' Thermal cost function
#' @description Calculates the normalized thermal cost based on F0-T curve
#'
#' @inheritParams hydraulic_cost
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param RH Relative humidity, in \%
#' @param Patm Atmospheric pressure, kPa
#' @param u Wind speed above the leaf boundary layer, m s-1
#' @param leaf_width Leaf width, m
#' @param Tcrit Leaf temperature at which F0(T) transitions from slow- to fast- rise, deg C
#' @param T50 Leaf temperature midway between Tcrit and Tmax, deg C
#' @param constant_kmax TRUE if the kmax does not vary with temperature for
#'     simulations; else FALSE
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
                        Tcrit = 50,
                        T50  = 51,
                        constant_kmax = FALSE
                        )
{
  E = trans_from_vc(P, kmax_25, T_air, b, c, constant_kmax)
  T_leaf = calc_Tleaf(T_air, PPFD, RH, E, u, Patm, leaf_width)

  r = 2 / (T50 - Tcrit)

  cost = 1 / (1 + exp(-r * (T_leaf - T50)))
  return(cost)
}


#' Thermal cost function, old version (use `thermal_cost` instead)
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
#' @param constant_kmax TRUE if the kmax does not vary with temperature for
#'     simulations; else FALSE
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
thermal_cost_v0 = function(P,
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
                        r = 4,
                        constant_kmax = FALSE
)
{
  E = trans_from_vc(P, kmax_25, T_air, b, c, constant_kmax)
  T_leaf = calc_Tleaf(T_air, PPFD, RH, E, u, Patm, leaf_width)

  F0 = F0_func(T_leaf, T50, F0_max, F0_min, r)

  cost = (F0 - F0_min) / (F0_max - F0_min)
  return(cost)
}


#' Carbon gain
#' @description Calculates the normalized carbon gain as described in Sperry et
#'     al. 2016
#'
#' @inheritParams thermal_cost
#' @param Ca Atmospheric CO2 concentration (ppm)
#' @param Jmax Maximum rate of electron transport at 25 deg C (mu mol m-2 s-1)
#' @param Vcmax Maximum carboxylation rate at 25 deg C (mu mol m-2 s-1)
#' @param constant_kmax TRUE if the kmax does not vary with temperature for
#'     simulations; else FALSE
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
                  RH = 60,
                  Ca = 420,
                  Jmax = 100,
                  Vcmax = 50,
                  constant_kmax = FALSE
                  )
{
  E = trans_from_vc(P, kmax_25, T_air, b, c, constant_kmax)
  A = calc_A(T_air, PPFD, Patm, E, u, leaf_width, RH, Ca, Jmax, Vcmax)

  Amax = max(A)

  gain = A/Amax
  return(gain)
}

#' Calculate all costs and gains
#' @description Calculate the carbon gain, hydraulic cost, and thermal costs.
#'
#' @inheritParams hydraulic_cost
#' @inheritParams thermal_cost
#' @inheritParams C_gain
#' @return A data frame with columns "P", "ID" and "cost_gain", where "P" denotes
#'     the leaf water potential, "ID" is one of CG, HC, or TC, and "cost_gain" is the
#'     corresponding value.
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
#' calc_costgain(P, b, c)
calc_costgain = function(
    P,
    b = -2.5,
    c = 2,
    kmax_25 = 4,
    T_air = 25,
    ratiocrit = 0.05,
    PPFD = 1000,
    RH = 90,
    Patm = 101.325,
    u = 2,
    leaf_width = 0.01,
    Tcrit = 50,
    T50 = 51,
    Ca = 420,
    Jmax = 100,
    Vcmax = 50,
    constant_kmax = FALSE
)
{
  HC = hydraulic_cost(P, b, c, kmax_25, T_air, ratiocrit, constant_kmax)
  TC = thermal_cost(P, b, c, kmax_25, T_air, PPFD, RH, Patm, u, leaf_width, Tcrit, T50, constant_kmax)
  CG = C_gain(P, b, c, kmax_25, T_air, PPFD, Patm, u, leaf_width, RH, Ca, Jmax, Vcmax, constant_kmax)

  cost_gain = c(HC, TC, CG)
  ID = c(rep("HC", length(HC)),
         rep("TC", length(TC)),
         rep("CG", length(CG)))
  df = data.frame(P, ID, cost_gain)
  return(df)
}

#' Calculate gain minus costs
#' @description Calculates the difference between the carbon gain and either
#'     hydraulic cost or summed hydraulic and thermal costs.
#'
#' @inheritParams calc_costgain
#'
#' @return A data frame with columns "P", "ID" and "CG_min_C", where "P" denotes
#'     the leaf water potential, "ID" is either HC or CC, and "CG_min_C" is the
#'     corresponding value of carbon gain minus the cost identified in "ID".
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
#' gain_min_costs(P, b, c)
gain_min_costs = function(
    P,
    b = -2.5,
    c = 2,
    kmax_25 = 4,
    T_air = 25,
    ratiocrit = 0.05,
    PPFD = 1000,
    RH = 90,
    Patm = 101.325,
    u = 2,
    leaf_width = 0.01,
    Tcrit = 50,
    T50 = 51,
    Ca = 420,
    Jmax = 100,
    Vcmax = 50,
    constant_kmax = FALSE
    )
{
  HC = hydraulic_cost(P, b, c, kmax_25, T_air, ratiocrit, constant_kmax)
  TC = thermal_cost(P, b, c, kmax_25, T_air, PPFD, RH, Patm, u, leaf_width, Tcrit, T50, constant_kmax)
  CG = C_gain(P, b, c, kmax_25, T_air, PPFD, Patm, u, leaf_width, RH, Ca, Jmax, Vcmax, constant_kmax)
  CC = HC + TC

  CG_min_HC = CG - HC
  CG_min_CC = CG - CC

  CG_min_C = c(CG_min_HC, CG_min_CC)
  ID = c(rep("HC", length(CG_min_HC)),
         rep("CC", length(CG_min_CC)))

  df = data.frame(P, ID, CG_min_C)
  return(df)
}



#' Calculate combined costs
#' @description Calculates the summed hydraulic and thermal costs.
#'
#' @inheritParams calc_costgain
#'
#' @return A data frame with columns "P", "ID" and "cost", where "P" denotes
#'     the leaf water potential, "ID" is either HC, TC, or CC, and "cost" is the
#'     corresponding cost value identified in "ID".
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
#' combined_costs(P, b, c)
combined_costs = function(
    P,
    b = -2.5,
    c = 2,
    kmax_25 = 4,
    T_air = 25,
    ratiocrit = 0.05,
    PPFD = 1000,
    RH = 90,
    Patm = 101.325,
    u = 2,
    leaf_width = 0.01,
    Tcrit = 50,
    T50 = 51,
    Ca = 420,
    Jmax = 100,
    Vcmax = 50,
    constant_kmax = FALSE
  )
{
  HC = hydraulic_cost(P, b, c, kmax_25, T_air, ratiocrit, constant_kmax)
  TC = thermal_cost(P, b, c, kmax_25, T_air, PPFD, RH, Patm, u, leaf_width, Tcrit, T50, constant_kmax)

  CC = HC + TC

  cost = c(HC, TC, CC)
  ID = c(rep("HC", length(HC)),
         rep("TC", length(TC)),
         rep("CC", length(CC)))

  df = data.frame(P, ID, cost)
  return(df)
}
