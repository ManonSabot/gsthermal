# Non-vectorized function for calculating leaf temperature
calc_Tleaf_fn = function (
    E = 0.2,
    T_air = 25,
    ...
)
{
  Tleaf <- try(uniroot(leaf_energy_balance,
                       interval = c(T_air - 15, T_air + 15), E = E, T_air = T_air,
                       ...)$root)
  return(Tleaf)
}

#' Leaf temperature
#' @description Calculate leaf temperature from energy balance. Adapted from
#'     '\code{plantecophys::FindTleaf}.
#'
#' @param T_air Air temperature, deg C
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param RH Relative humidity, in \%
#' @param E Transpiration rate, mmol m-2 s-1
#' @param u Wind speed above the leaf boundary layer, m s-1
#' @param Patm Atmospheric pressure, kPa
#' @param leaf_width Leaf width, m
#'
#' @return Leaf temperature, deg C
#'
#' @examples
#' # Calculate leaf temperature for scalar value of E
#' calc_Tleaf(E = 0.2, T_air = 25)
#'
#'
#' # Calculate leaf temperature along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#'
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#' trans = trans_from_vc(P = P, b = b, c = c) # Create vector of transpiration supply stream
#'
#' Tleaves = calc_Tleaf(E = trans) # Get leaf temperatures
#'
#' @export calc_Tleaf
#' @rdname calc_Tleaf
calc_Tleaf = function(
    E = 0.2,
    T_air = 25,
    PPFD = 1000,
    RH = 60,
    u = 2,
    Patm = 101.325,
    leaf_width = 0.01
)
{
  fn = Vectorize(calc_Tleaf_fn)
  out = fn(E = E, T_air = T_air, PPFD = PPFD, RH = RH, u = u, Patm = Patm,
           leaf_width = leaf_width)
  return(out)
}

#' Leaf energy balance
#'
#' @param T_leaf Leaf temperature, deg C
#' @param T_air Air temperature, deg C
#' @param E Transpiration rate, mmol m-2 s-1
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param RH Relative humidity, in \%
#' @param u Wind speed above the leaf boundary layer, m s-1
#' @param Patm Atmospheric pressure, kPa
#' @param leaf_width Leaf width, m
#' @param leaf_abs Leaf absorptance of solar radiation (0-1)
#'
#' @return Leaf temperature, deg C
#' @noRd
leaf_energy_balance = function(
    T_leaf = 25,
    T_air = 25,
    E = 0.2,
    PPFD = 1500,
    RH = 70,
    Patm = 101.325,
    u = 2, # m s-1
    leaf_width = 0.01, # m
    leaf_abs = 0.5,
    returnwhat = c("balance", "fluxes")
    )
{
  returnwhat <- match.arg(returnwhat)
  Boltz <- 5.67 * 10^-8
  Emissivity <- 0.95
  LatEvap <- 2.54
  CPAIR <- 1010
  H2OLV0 <- 2501000
  H2OMW <- 0.018
  AIRMA <- 0.029
  AIRDENS <- 1.204
  UMOLPERJ <- 4.57
  DHEAT <- 2.15e-05
  Tair_k <- T_air + 273.15
  Tleaf_k <- T_leaf + 273.15
  AIRDENS <- Patm * 1000/(287.058 * Tair_k)
  LHV <- (H2OLV0 - 2365 * T_air) * H2OMW
  SLOPE <- (plantecophys::esat(T_air + 0.1) - plantecophys::esat(T_air))/0.1
  Gradiation <- 4 * Boltz * Tair_k^3 * Emissivity/(CPAIR *
                                                     AIRMA)
  CMOLAR <- Patm * 1000/(8.314 * Tair_k)
  Gbhforced <- 0.003 * sqrt(u/leaf_width) * CMOLAR
  GRASHOF <- 1.6e+08 * abs(T_leaf - T_air) * (leaf_width^3)
  Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25)/leaf_width * CMOLAR
  Gbh <- 2 * (Gbhfree + Gbhforced)
  Rsol <- 2 * PPFD/UMOLPERJ
  VPD <- calc_VPD(T_air = T_air, RH = RH)
  ea <- plantecophys::esat(T_air) - 1000 * VPD
  ema <- 0.642 * (ea/Tair_k)^(1/7)
  Rnetiso <- leaf_abs * Rsol - (1 - ema) * Boltz * Tair_k^4
  GAMMA <- CPAIR * AIRMA * Patm * 1000/LHV
  ET <- E/1000
  lambdaET <- LHV * ET
  Y <- 1/(1 + Gradiation/Gbh)
  H2 <- Y * (Rnetiso - lambdaET)
  H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (T_air - T_leaf)
  Tleaf2 <- T_air + H2/(CPAIR * AIRDENS * (Gbh/CMOLAR))
  EnergyBal <- T_leaf - Tleaf2
  if (returnwhat == "balance")
    return(EnergyBal)
  if (returnwhat == "fluxes") {
    l <- data.frame(ELEAFeb = 1000 * ET, Gradiation = Gradiation,
                    Rsol = Rsol, Rnetiso = Rnetiso, H = H, lambdaET = lambdaET,
                    Gbh = Gbh, H2 = H2, Tleaf2 = Tleaf2)
    return(l)
  }
}

#' Leaf temperature, old version
#'
#' @param T_air Air temperature, deg C
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param RH Relative humidity, in \%
#' @param E Transpiration rate, mmol m-2 s-1
#' @param u Wind speed above the leaf boundary layer, m s-1
#' @param Patm Atmospheric pressure, kPa
#' @param leaf_width Leaf width, m
#'
#' @return Leaf temperature, deg C
#' @export
#'
#' @examples
#' # Calculate leaf temperature for scalar value of E
#' calc_Tleaf(E = 10)
#'
#'
#' # Calculate leaf temperature along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#'
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#' trans = trans_from_vc(P = P, b = b, c = c) # Create vector of transpiration supply stream
#'
#' Tleaves = calc_Tleaf0(E = trans) # Get leaf temperatures
calc_Tleaf0 = function(
    T_air = 25,
    PPFD = 1000,
    RH = 60,
    E,
    u = 2,
    Patm = 101.325,
    leaf_width = 0.01
    )
{
  # Constants
  C_p = 29.3 # Specific heat of air at constant pressure, J mol^-1 deg C^-1
  LH2O = 2.501e6  # latent heat H2O (J kg-1)
  MH = 1.00794  # molar mass of hydrogen (g mol-1)
  MO = 15.999  # molar mass of O (g mol-1)
  MH2O = 2. * MH + MO  # H2O molar mass (g mol-1)

  # Radiative conductance, mol m-2 s-1
  g_r = calc_gr(T_air)

  g_Ha = calc_gHa(T_air, u, Patm, leaf_width)

  # Atmospheric VPD (kPa)
  VPD = calc_VPD(T_air, RH)

  # Net radiation (W m-2)
  R_net = calc_Rnet(T_air, PPFD, VPD)

  # Latent heat of vaporization, J mol-1
  lambda = (LH2O - 2.365e3 * T_air) * MH2O * 1.e-3

  # Canopy / leaf sensible heat flux, W m-2
  H = R_net - lambda * E * 1e-3

  Tleaf = T_air + H / (C_p * (g_Ha + g_r))
  return(Tleaf)
}




#' Vapor pressure deficit
#'
#' @param T_leaf Leaf temperature, deg C
#' @param T_air Air temperature, deg C
#' @param RH Relative humidity, in \%
#'
#' @return Leaf vapor pressure deficit, kPa
#' @export
#'
#' @examples
#' # Calculate Tleaf for a given scalar value of leaf temperature
#' calc_Dleaf(T_leaf = 30)
#'
#'
#' # Calculate leaf VPD along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#'
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#' trans = trans_from_vc(P = P, b = b, c = c) # Create vector of transpiration supply stream
#' Tleaves = calc_Tleaf(E = trans) # Get leaf temperatures
#'
#' Dleaves = calc_Dleaf(T_leaf = Tleaves)
calc_Dleaf = function(T_leaf = 25,
                      T_air = 25,
                      RH = 60
)
{
  esat_l = vpsat(T_leaf)
  esat_a = vpsat(T_air)
  D_leaf = esat_l - (RH/100)*(esat_a)
  return(D_leaf)
}




#' Stomatal conductance
#'
#' @param E Transpiration rate, mmol m-2 s-1
#' @param D_leaf Leaf VPD, kPa
#' @param Patm Atmospheric pressure, kPa
#'
#' @return Stomatal conductance to water vapor, mol m-2 s-1
#' @export
#'
#' @examples
#' # Calculate leaf conductance for a scalar value of E
#' calc_gw(E = 3)
#'
#' # Calculate gw along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#'
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#' trans = trans_from_vc(P = P, b = b, c = c) # Create vector of transpiration supply stream
#' Tleaves = calc_Tleaf(E = trans) # Get leaf temperatures
#' Dleaves = calc_Dleaf(T_leaf = Tleaves) # Get leaf VPD
#'
#' calc_gw(E = trans, D_leaf = Dleaves)
calc_gw = function(
    E,
    D_leaf = 1.5,
    Patm = 101.325
)
{
  g_w = Patm * E * 1e-3 / D_leaf
  return(g_w)
}




#' Calculate photosynthesis
#'
#' @param T_air Air temperature, deg C
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param Patm Atmospheric pressure, kPa
#' @param E Transpiration rate, mmol m-2 s-1
#' @param u Wind speed above the leaf boundary layer, m s-1
#' @param leaf_width Leaf width, m
#' @param RH Relative humidity, in \%
#' @param Ca Atmospheric CO2 concentration (ppm)
#' @param Jmax Maximum rate of electron transport at 25 deg C (mu mol m-2 s-1)
#' @param Vcmax Maximum carboxylation rate at 25 deg C (mu mol m-2 s-1)
#' @param net TRUE if desired output is net photosynthesis; FALSE if desired output
#'     is gross photosynthesis
#' @param Rd0 Day respiration rate at reference temperature (TrefR). Must be a positive value.
#' @param TrefR Reference temperature for Rd (deg C)
#' @param netOrig TRUE if net photosynthesis is to be calculated within
#'     plantecophys::Photosyn. FALSE if Anet is to be calculated as Agross - Rd.
#'
#' @return Photosynthetic rate (mol m-2 s-1)
#' @export
#'
#' @examples
#' # Calculate A for a scalar value of E
#' calc_A(E = 10)
#'
#' # Calculate leaf VPD along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#'
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#' trans = trans_from_vc(P = P, b = b, c = c) # Create vector of transpiration supply stream
#'
#' calc_A(E = trans)
calc_A = function(T_air = 25,
                  PPFD = 1000,
                  Patm = 101.325,
                  E,
                  u = 2,
                  leaf_width = 0.01,
                  RH = 60,
                  Ca = 420,
                  Jmax = 100,
                  Vcmax = 50,
                  net = FALSE,
                  Rd0 = 0.92,
                  TrefR = 25,
                  netOrig = TRUE
)
{
  T_leaf = calc_Tleaf(T_air = T_air, PPFD = PPFD, RH = RH, E = E, u = u,
                      Patm = Patm, leaf_width = leaf_width)
  D_leaf = calc_Dleaf(T_leaf, T_air, RH)
  g_w = calc_gw(E, D_leaf, Patm)

  if (net == FALSE){
  Photosyn_out = mapply(plantecophys::Photosyn,
                        VPD = D_leaf, Ca = Ca, PPFD = PPFD, Tleaf = T_leaf,
                        Patm = Patm, GS = g_w, Rd0 = 0, Jmax = Jmax, Vcmax = Vcmax)
  A = as.numeric(Photosyn_out[2,])
   } else {
    Rd = Rd0 * exp(0.1012 * (T_leaf - TrefR) - 0.0005 * (T_leaf**2 - TrefR**2))

    if (isTRUE(netOrig)) {
      Photosyn_out = mapply(plantecophys::Photosyn,
                            VPD = D_leaf, Ca = Ca, PPFD = PPFD, Tleaf = T_leaf,
                            Patm = Patm, GS = g_w, Rd = Rd, Jmax = Jmax, Vcmax = Vcmax)
      A = as.numeric(Photosyn_out[2,])
      } else {
      Photosyn_out = mapply(plantecophys::Photosyn,
                            VPD = D_leaf, Ca = Ca, PPFD = PPFD, Tleaf = T_leaf,
                            Patm = Patm, GS = g_w, Rd0 = 0, Jmax = Jmax, Vcmax = Vcmax)
      A = as.numeric(Photosyn_out[2,]) - Rd
    }
  }
  return(A)
}

#' Radiative conductance
#'
#' @param T_air Air temperature, deg C
#' @param epsilon Leaf emissivity
#' @export
#' @return Radiative conductance, mol m-2 s-1
#' @examples
#' calc_gr(T_air = 30)
calc_gr = function(
    T_air = 25,
    epsilon = 0.97
)
{
  # Constants
  sigma = 5.67e-8 # Stefan-Boltzman constant, W m^-2 K^-4
  C_p = 29.3 # Specific heat of air at constant pressure, J mol^-1 deg C^-1

  TairK = T_air + 273 # air temperature in Kelvin
  g_r = max(0., 4 * epsilon * sigma * TairK ** 3 / C_p)
  return(g_r) # mol m-2 s-1
}



#' Boundary layer conductance to forced convection
#'
#' @param T_air Air temperature, deg C
#' @param u Wind speed above the leaf boundary layer, m s-1
#' @param Patm Atmospheric pressure, kPa
#' @param leaf_width Leaf width, m
#'
#' @export
#'
#' @return Boundary layer conductance to forced convection, mol m-2 s-1
calc_gHa = function(
    T_air = 25,
    u = 2,
    Patm = 101.325,
    leaf_width = 0.01
)
{
  # Constants
  Mair = 28.9644  # molar mass of dry air (g mol-1)
  DH = 21.5e-6  # molecular diffusivity to heat (m2 s-1)
  R = 8.3144598  # molar ideal gas constant (J mol-1 K-1)

  TairK = T_air + 273 # air temperature in Kelvin
  cmolar = Patm * 1e3 / (R * TairK)  # air molar density

  # Sutherland Eq for dynamic viscosity
  mu = 1.458e-6 * TairK ** 1.5 / (TairK + 110.4)  # Pa s

  # Kinematic viscosity
  nu = mu * R * TairK / (Patm * Mair)  # m2 s-1
  prandtl = nu / DH  # unitless

  # Boundary layer cond to forced convect. (Eq. 7.29 Campbell & Norman, 1998)
  d = 0.72 * leaf_width  # leaf width, m
  reynolds = u * d / nu  # unitless
  g_Ha = (0.664 * cmolar * DH * (reynolds ** 0.5) * (prandtl ** (1. / 3.))
          / d)
}


#' Daytime respiration rate
#' @description Calculates the daytime respiration rate (Rd) with the Heskel et
#'     al. 2016 global polynomial model.
#'
#' @inheritParams calc_Tleaf
#' @param Rd0 Day respiration rate (mu mol m-2 s-1) at reference temperature
#'     (TrefR). Must be a positive value.
#' @param TrefR Reference temperature for Rd (deg C)
#'
#' @return Rd (mu mol m-2 s-1)
#' @export
#'
#' @examples
#' # Calculate daytime respiration for scalar value of E
#' calc_Rd(E = 10)
#'
#'
#' # Calculate leaf temperature along transpiration supply stream
#' Weibull = fit_Weibull() # Fit Weibull parameters
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c) # Calculate Pcrit based on Weibull curve
#'
#' P = Ps_to_Pcrit(Pcrit = Pcrit) # Create Ps to Pcrit vector
#' E = trans_from_vc(P = P, b = b, c = c) # Create vector of transpiration supply stream
#'
#' Rd_vect = calc_Rd(E = E) # Get leaf temperatures
calc_Rd = function(T_air = 25, PPFD = 1000, RH = 60, E, u = 2, Patm = 101.325,
                   leaf_width = 0.01, Rd0 = 0.92, TrefR = 25)
{
  T_leaf = calc_Tleaf(T_air = T_air, PPFD = PPFD, RH = RH, E = E, u = u,
                      Patm = Patm, leaf_width = leaf_width)
  Rd = Rd0 * exp(0.1012 * (T_leaf - TrefR) - 5e-04 * (T_leaf^2 - TrefR^2))
  return(Rd)
}
