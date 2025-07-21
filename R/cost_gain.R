#' Hydraulic cost function
#' @description Calculates the normalized hydraulic cost as described in Sperry
#'     et al. 2016
#'
#' @param P Vector of equally spaced water potentials ranging from Ps to Pcrit,
#'     -MPa
#' @param b Weibull scale parameter
#' @param c Weibull shape parameter
#' @param kmax_25 Max plant conductance at 25 deg C, mmol s-1 m-2 MPa-1
#' @param Tair Air temperature, deg C
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
                          Tair = 25,
                          ratiocrit = 0.05,
                          constant_kmax = FALSE
)
{
  # Adjust kmax based on the temperature response of viscosity
  kmax = calc_kmax(kmax_25, Tair, constant_kmax)

  # Calculate conductivity at Pcrit
  kcrit = ratiocrit*kmax

  # Calculate conductivity over the range of possible water potentials
  k = kmax * vulnerability_curve(P, b, c)

  # Calculate the hydraulic cost
  kmax_i = max(k)
  cost = (kmax_i - k) / (kmax_i - kcrit)
  return(cost)
}


#' Thermal cost function
#' @description Calculates the normalized thermal cost based on F0-T curve
#'
#' @inheritParams hydraulic_cost
#' @param VPD Air vapor pressure deficit, kPa
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param Patm Atmospheric pressure, kPa
#' @param Wind Wind speed above the leaf boundary layer, m s-1
#' @param Wleaf Leaf width, m
#' @param LeafAbs Leaf absorptance of solar radiation (0-1)
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
                        Tair = 25,
                        VPD = 1.5,
                        PPFD = 1000,
                        Patm = 101.325,
                        Wind = 2,
                        Wleaf = 0.025,
                        LeafAbs = 0.5,
                        Tcrit = 50,
                        T50  = 51,
                        constant_kmax = FALSE
                        )
{
  # Calculate Tleaf over the transpiration supply stream
  E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax)
  Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, E = E, Wind = Wind,
                      Patm = Patm, Wleaf = Wleaf, LeafAbs = LeafAbs)

  # Calculate thermal cost
  r = 2 / (T50 - Tcrit)
  cost = 1 / (1 + exp(-r * (Tleaf - T50)))
  return(cost)
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
                  kmax_25 = 0.5,
                  Tair = 25,
                  VPD = 1.5,
                  PPFD = 1000,
                  Patm = 101.325,
                  Wind = 2,
                  Wleaf = 0.025,
                  LeafAbs = 0.5,
                  Ca = 420,
                  Jmax = 100,
                  Vcmax = 50,
                  constant_kmax = FALSE,
                  net = FALSE,
                  Rd0 = 0.92,
                  TrefR = 25,
                  netOrig = TRUE,
                  ...
                  )
{
  E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax)
  A = calc_A(Tair, VPD, PPFD, Patm, E, Wind, Wleaf, LeafAbs,
             Ca, Jmax, Vcmax, net, Rd0, TrefR, netOrig, ...)
  E = ifelse(E == 0, NA, E)
  A = ifelse(E == 0, NA, A)
  Amax = if (is.null(Amax)) {
    max(abs(A[!is.na(A)]))
  } else {
    Amax
  }

  gain =
    if (!is.na(Amax) & Amax != 0) {
      A/Amax
    } else {
      rep(0, length = length(A))
    }
  return(gain)
}

#' Carbon gain
#' @description Calculates the normalized carbon gain with constrained Ci
#'
#' @inheritParams C_gain
#'
#' @return Normalized carbon gain
#' @export
C_gain_alt = function (P,
                       b = -2.5,
                       c = 2,
                       Amax = NULL,
                       kmax_25 = 4,
                       Tair = 25,
                       VPD = 1.5,
                       PPFD = 1000,
                       Patm = 101.325,
                       Wind = 2,
                       Wleaf = 0.01,
                       LeafAbs = 0.5,
                       Ca = 420,
                       Jmax = 100,
                       Vcmax = 50,
                       constant_kmax = FALSE,
                       Rd0 = 0.92,
                       TrefR = 25,
                       net = TRUE,
                       ...)
{
  # Calculate transpiration supply stream
  E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax)

  # Create 2D array of possible Ci's for each value of E
  Cis = array(seq(1, Ca, length.out = 500), dim = c(500, length(E)))

  # Get leaf temperature and stomatal conductance for each value of E
  Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD,
                     E = E, Wind = Wind, Patm = Patm, Wleaf = Wleaf,
                     LeafAbs = LeafAbs)
  g_w = calc_gw(E, Tleaf, Patm, Tair, VPD, PPFD, Wind,
                Wleaf)
  Dleaf = plantecophys::VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)

  # Calculate supply A from Fick's law
  g_ws = t(array(g_w, dim = c(500, length(E))))
  A_supply = g_ws/1.6 * (Ca - Cis)

  Es = t(array(E, dim = c(500, length(E))))

  # Calculate demand A with Farqhuar model
  Tleaves = t(array(Tleaf, dim = c(500, length(E))))
  Dleaves = t(array(Dleaf, dim = c(500, length(E))))
  Photosyn_out = mapply(plantecophys::Photosyn, VPD = Dleaves,
                        Ca = Ca, PPFD = PPFD, Tleaf = Tleaves,
                        Patm = Patm,
                        Ci = Cis, Jmax = Jmax, Vcmax = Vcmax,
                        Rd0 = Rd0, TrefR = TrefR, new_JT = FALSE,
                        ...)

  if (isTRUE(net)) {
    A_demand = as.numeric(Photosyn_out[2,])
  } else {
    Anet = as.numeric(Photosyn_out)[2,]
    Rd = as.numeric(Photosyn_out)[8,]
    A_demand = Anet + Rd
  }
  dim(A_demand) = c(500, 500)

  # Find which values of supply and demand A match closest for each value of E
  idx = apply(abs(A_supply - A_demand), 2, FUN = which.min)

  # Select the corresponding values of Ci and A
  Ci = sapply(1:length(E), function(e) {Cis[idx[e], e]})
  A_P = sapply(1:length(E), function(e) {A_demand[idx[e], e]})

  # Calculate Amax
  Amax = if (is.null(Amax) & !(all(A_P <= 0))) {
    max(abs(A_P))
  } else if (is.null(Amax) & (all(A_P <= 0))) {
    max(abs(A_P))
  } else {
    Amax
  }

  # Normalise A to get CG
  gain = A_P/Amax
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
#' @param EaV,EdVC,delsC Vcmax temperature response parameters
#' @param EaJ,EdVJ,delsJ Jmax temperature response parameters
#' @param constr_Ci If TRUE, constrains Ci between 0 and Ca.
#' @return A data frame with columns "P", "ID" and "cost_gain", where "P" denotes
#'     the leaf water potential, "ID" is one of CG, HC, or TC, and "cost_gain" is the
#'     corresponding value.
#' @export
calc_costgain = function(
    P,
    b = -2.5,
    c = 2,
    Amax_gross = NULL,
    Amax_net = NULL,
    kmax_25 = 4,
    Tair = 25,
    VPD = 1.5,
    ratiocrit = 0.05,
    PPFD = 1000,
    Patm = 101.325,
    Wind = 2,
    Wleaf = 0.025,
    LeafAbs = 0.5,
    Tcrit = 50,
    T50 = 51,
    Ca = 420,
    Rd0 = 0.92,
    TrefR = 25,
    Vcmax=34,
    EaV=62307,
    EdVC=2e5,
    delsC=639,
    Jmax = 60,
    EaJ=33115,
    EdVJ=2e5,
    delsJ=635,
    constr_Ci = FALSE,
    ...
)
{
  if(is.null(P)) {
    return(NULL)
  }

  # Hydraulic costs
  HC_constkmax = hydraulic_cost(P, b, c, kmax_25, Tair, constant_kmax = TRUE)
  HC_varkmax = hydraulic_cost(P, b, c, kmax_25, Tair, constant_kmax = FALSE)

  # Thermal cost
  TC = thermal_cost(P, b, c, kmax_25, Tair, VPD, PPFD, Patm,
                    Wind, Wleaf, LeafAbs, Tcrit, T50, constant_kmax = FALSE)

  # Carbon gains
  if(isFALSE(constr_Ci)) {
    CG_gross_constkmax = C_gain(
      P, b, c, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD,
      VPD = VPD, Tcrit = Tcrit, T50 = T50,
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = TRUE, net = FALSE, new_JT = FALSE,
      ...)
    CG_gross_varkmax = C_gain(
      P, b, c, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD,
      VPD = VPD, Tcrit = Tcrit, T50 = T50,
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, net = FALSE, new_JT = FALSE,
      ...)
    CG_net = C_gain(
      P, b, c, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD,
      VPD = VPD, Tcrit = Tcrit, T50 = T50,
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, net = TRUE, netOrig = TRUE, new_JT = FALSE,
      ...)
    CG_net_newJT = C_gain(
      P, b, c, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD,
      VPD = VPD, Tcrit = Tcrit, T50 = T50,
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, net = TRUE, netOrig = TRUE, new_JT = TRUE,
      ...)
  } else {
    CG_gross_constkmax = C_gain_alt(
      P, b, c, Amax = NULL, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD,
      VPD = VPD, Tcrit = Tcrit, T50 = T50,
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = TRUE, new_JT = FALSE, net = FALSE,
      ...)
    CG_gross_varkmax = C_gain_alt(
      P, b, c, Amax = NULL, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD,
      VPD = VPD, Tcrit = Tcrit, T50 = T50,
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, new_JT = FALSE, net = FALSE,
      ...)
    CG_net = C_gain_alt(
      P, b, c, Amax = NULL, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD,
      VPD = VPD, Tcrit = Tcrit, T50 = T50,
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, new_JT = FALSE, net = TRUE,
      ...)
    CG_net_newJT = C_gain_alt(
      P, b, c, Amax = NULL, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD,
      VPD = VPD, Tcrit = Tcrit, T50 = T50,
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, new_JT = TRUE, net = TRUE,
      ...)
  }

  df = data.frame(P,
                  HC_constkmax, HC_varkmax,
                  TC,
                  CG_gross_constkmax, CG_gross_varkmax,
                  CG_net, CG_net_newJT)
  return(df)
}
