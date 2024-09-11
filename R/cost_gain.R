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
#' # Calculate leaf water potential along transpiration supply stream
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


#' Respiratory cost function
#' @description Calculates the normalized respiratory cost along the
#'     transpiration supply stream.
#'
#' @inheritParams calc_Rd
#' @inheritParams C_gain
#'
#' @return Respiratory cost, unitless
#' @export
#'
#' @examples
#' # Calculate leaf water potential along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#'
#' # Calculate respiratory cost
#' respiratory_cost(P, b, c)
respiratory_cost = function(P,
                            b = -2.5,
                            c = 2,
                            Amax = NULL,
                            kmax_25 = 4,
                            T_air = 25,
                            PPFD = 1000,
                            Patm = 101.325,
                            u = 2,
                            leaf_width = 0.01,
                            RH = 60,
                            constant_kmax = FALSE,
                            Rd0 = 0.92,
                            TrefR = 25)
{
  E = trans_from_vc(P, kmax_25, T_air, b, c, constant_kmax)
  Rd = calc_Rd(T_air, PPFD, RH, E, u, Patm, leaf_width, Rd0, TrefR)
  Rd_min = min(Rd)

  if (is.null(Amax)) {
    Rd_max = max(Rd)
    cost = (Rd - Rd_min) / (Rd_max - Rd_min)
  } else {
    cost = (Rd - Rd_min) / (Amax - Rd_min)
  }
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
  T_leaf = calc_Tleaf(T_air = T_air, PPFD = PPFD, RH = RH, E = E, u = u,
                      Patm = Patm, leaf_width = leaf_width)

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

#' Maximum potential photosynthetic rate
#' @description Calculates the maximum potential photosynthetic rate for the given
#'     hydraulic parameters and range of water potentials.
#'
#' @inheritParams C_gain
#' @param Tair_range Range of air temperatures considered for determination of
#'     Amax (deg C)
#'
#' @return Maximum potential photosynthetic rate used to normalize the carbon
#'     gain function and the temperature at which it occurs.
#' @export
#'
#' @examples
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#'
#' Amax_overT(P, b, c, constant_kmax = TRUE)
Amax_overT = function(P,
                      b = -2.5,
                      c = 2,
                      kmax_25 = 4,
                      PPFD = 1000,
                      Patm = 101.325,
                      u = 2,
                      leaf_width = 0.01,
                      RH = 60,
                      Ca = 420,
                      Jmax = 100,
                      Vcmax = 50,
                      Tair_range= seq(20, 50, by = 0.5),
                      constant_kmax = FALSE,
                      net = FALSE,
                      Rd0 = 0.92,
                      TrefR = 25,
                      netOrig = TRUE
                      )
{
  # Calculate the maximum transpiration rate given the vulnerability curve parameters
  # for the specified temperature range
  VC = vulnerability_curve(P, b, c)
  AUC = pracma::trapz(P, vulnerability_curve(P, b, c))
  kmax = calc_kmax(kmax_25, Tair_range, constant_kmax)
  E = kmax * AUC

  # Calculate A over the temperature range
  calc_A_vect = Vectorize(calc_A, c("T_air", "E"))
  A = calc_A_vect(Tair_range, PPFD, Patm, E, u, leaf_width, RH, Ca, Jmax,
                  Vcmax, net, Rd0, TrefR, netOrig)

  # Select the maximum photosynthetic rate and corresponding temperature
  i = ifelse(all(A <= 0), which.max(abs(A)), which.max(A))
  #i = which.max(A)
  Amax = A[i]
  Tair_opt = Tair_range[i]
  return(data.frame(Tair_opt, Amax))
}

#' Carbon gain
#' @description Calculates the normalized carbon gain as described in Sperry et
#'     al. 2016
#'
#' @inheritParams thermal_cost
#' @inheritParams calc_A
#' @param Amax Maximum potential photosynthetic rate (A) that A values are normalized
#'     by to obtain the carbon gain. If NULL, this is the instantaneous maximum
#'     potential rate described by Sperry et al. (2016).

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
                  Amax = NULL,
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
                  constant_kmax = FALSE,
                  net = FALSE,
                  Rd0 = 0.92,
                  TrefR = 25,
                  netOrig = TRUE
                  )
{
  # Calculate the photosynthesis over the transpiration supply stream
  E = trans_from_vc(P, kmax_25, T_air, b, c, constant_kmax)
  A = calc_A(T_air, PPFD, Patm, E, u, leaf_width, RH, Ca, Jmax, Vcmax, net, Rd0,
             TrefR, netOrig)

  # Calculate Amax if not provided
  #Amax = ifelse(is.null(Amax), abs(max(A)), Amax)
  Amax = if (is.null(Amax) & !(all(A <= 0))){
    max(A)
  } else if (is.null(Amax) & (all(A <= 0))) {
    max(abs(A))
  } else {
    Amax
  }

  gain = A/Amax
  return(gain)
}

#' Calculate all costs and gains
#' @description Calculate the carbon gain, hydraulic cost, and thermal costs.
#'
#' @inheritParams hydraulic_cost
#' @inheritParams thermal_cost
#' @inheritParams C_gain
#' @param Amax_net Maximum potential net photosynthetic rate (A) that A values are normalized
#'     by to obtain the carbon gain. If NULL, this is the instantaneous maximum
#'     potential rate described by Sperry et al. (2016).
#' @param Amax_gross Maximum potential gross photosynthetic rate (A) that A values are normalized
#'     by to obtain the carbon gain. If NULL, this is the instantaneous maximum
#'     potential rate described by Sperry et al. (2016).
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
    Amax_gross = NULL,
    Amax_net = NULL,
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
    constant_kmax = FALSE,
    Rd0 = 0.92,
    TrefR = 25
)
{
  HC = hydraulic_cost(P, b, c, kmax_25, T_air, ratiocrit, constant_kmax)
  TC = thermal_cost(P, b, c, kmax_25, T_air, PPFD, RH, Patm, u, leaf_width, Tcrit, T50, constant_kmax)
  CG_net = C_gain(P, b, c, Amax_net, kmax_25, T_air, PPFD, Patm, u, leaf_width,
              RH, Ca, Jmax, Vcmax, constant_kmax, net = TRUE, Rd0, TrefR, netOrig = FALSE)
  CG_gross = C_gain(P, b, c, Amax_gross, kmax_25, T_air, PPFD, Patm, u, leaf_width,
                    RH, Ca, Jmax, Vcmax, constant_kmax, net = FALSE, Rd0, TrefR)

  cost_gain = c(HC, TC, CG_net, CG_gross)
  ID = c(rep("HC", length(HC)),
         rep("TC", length(TC)),
         rep("CG_net", length(CG_net)),
         rep("CG_gross", length(CG_gross)))
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
    Amax_gross = NULL,
    Amax_net = NULL,
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
    constant_kmax = FALSE,
    Rd0 = 0.92,
    TrefR = 25
    )
{
  HC = hydraulic_cost(P, b, c, kmax_25, T_air, ratiocrit, constant_kmax)
  TC = thermal_cost(P, b, c, kmax_25, T_air, PPFD, RH, Patm, u, leaf_width, Tcrit, T50, constant_kmax)
  CG_net = C_gain(P, b, c, Amax_net, kmax_25, T_air, PPFD, Patm, u, leaf_width,
              RH, Ca, Jmax, Vcmax, constant_kmax, net = TRUE, Rd0, TrefR, netOrig = FALSE)
  CG_gross = C_gain(P, b, c, Amax_gross, kmax_25, T_air, PPFD, Patm, u, leaf_width,
                  RH, Ca, Jmax, Vcmax, constant_kmax, net = FALSE, Rd0, TrefR)
  CC = HC + TC

  CGnet_min_HC = CG_net - HC
  CGnet_min_CC = CG_net - CC
  CGgross_min_HC = CG_gross - HC
  CGgross_min_CC = CG_gross - CC

  CG_min_C = c(CGnet_min_HC, CGnet_min_CC, CGgross_min_HC, CGgross_min_CC)
  cost = c(rep("HC", length(CGnet_min_HC)),
         rep("CC", length(CGnet_min_CC)),
         rep("HC", length(CGgross_min_HC)),
         rep("CC", length(CGgross_min_CC))
         )
  A = c(rep("net", length(CGnet_min_HC)),
           rep("net", length(CGnet_min_CC)),
           rep("gross", length(CGgross_min_HC)),
           rep("gross", length(CGgross_min_CC))
  )

  df = data.frame(P, CG_min_C, cost, A)
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
