# Non-vectorized function for calculating leaf temperature
calc_Tleaf_fn <- function(
    Tair = 25,
    E = 0.2,
    VPD = 1.5,
    PPFD = 1000,
    ...,
    initial_width = 30,
    max_attempts = 50,
    width_step = 2
) {
  fun <- function(Tleaf) leaf_energy_balance(
    Tleaf = Tleaf,
    Tair = Tair,
    E = E,
    VPD = VPD,
    PPFD = PPFD,
    ...
  )

  width <- initial_width
  attempt <- 1
  success <- FALSE

  while (attempt <= max_attempts) {
    lower <- Tair - width
    upper <- Tair + width

    f.lower <- fun(lower)
    f.upper <- fun(upper)

    if (!is.finite(f.lower) || !is.finite(f.upper)) {
      # Shrink interval if function gives NA at either bound
      width <- width - width_step
      attempt <- attempt + 1
      next
    }

    if (sign(f.lower) != sign(f.upper)) {
      # Found valid interval
      success <- TRUE
      break
    }

    # Both ends are finite but same sign -> extend the interval
    width <- width + width_step
    attempt <- attempt + 1
  }

  if (!success) {
    warning("Could not find valid interval after ", max_attempts, " attempts (Tair = ", Tair, ", E = ", E, ")")
    return(NA)
  }

  # Now safely apply uniroot
  result <- try(uniroot(fun, interval = c(lower, upper)), silent = TRUE)

  if (inherits(result, "try-error")) {
    warning("uniroot failed: ", conditionMessage(attr(result, "condition")))
    return(NA)
  } else {
    return(result$root)
  }
}


#' Leaf temperature
#' @description Calculate leaf temperature from energy balance. Adapted from
#'     '\code{plantecophys::FindTleaf}.
#'
#' @param Tair Air temperature, deg C
#' @param VPD Air vapor pressure deficit, kPa
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param E Transpiration rate, mmol m-2 s-1
#' @param Wind Wind speed above the leaf boundary layer, m s-1
#' @param Patm Atmospheric pressure, kPa
#' @param Wleaf Leaf width, m
#' @param LeafAbs Leaf absorptance of solar radiation (0-1)
#'
#' @return Leaf temperature, deg C
#'
#' @examples
#' # Calculate leaf temperature for scalar value of E
#' calc_Tleaf(E = 0.2, Tair = 25)
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
    Tair = 25,
    VPD = 1.5,
    E = 0.2,
    PPFD = 1000,
    Wind = 2,
    Patm = 101.325,
    Wleaf = 0.025,
    LeafAbs = 0.5
)
{
  fn = Vectorize(calc_Tleaf_fn)
  out = fn(E = E, Tair = Tair, PPFD = PPFD, VPD = VPD, Wind = Wind, Patm = Patm,
           Wleaf = Wleaf, LeafAbs = LeafAbs)
  return(out)
}

#' Leaf energy balance
#' @description
#' The net leaf energy balance, given that we know Tleaf, E. Modified from
#'     '\code{plantecophys::LeafEnergyBalance} to take E rather than gs as an
#'     input.
#'
#'
#' @param Tleaf Leaf temperature, deg C
#' @param Tair Air temperature, deg C
#' @param VPD Air vapor pressure deficit, kPa
#' @param E Transpiration rate, mmol m-2 s-1
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param Wind Wind speed above the leaf boundary layer, m s-1
#' @param Patm Atmospheric pressure, kPa
#' @param Wleaf Leaf width, m
#' @param LeafAbs Leaf absorptance of solar radiation (0-1)
#'
#' @return Leaf temperature, deg C
#' @noRd
leaf_energy_balance = function(
    Tleaf = 25,
    VPD = 1.5,
    Tair = 25,
    E = 0.2,
    PPFD = 1500,
    Patm = 101.325,
    Wind = 2, # m s-1
    Wleaf = 0.025, # m
    LeafAbs = 0.5,
    returnwhat = c("balance", "fluxes")
    )
{
  returnwhat <- match.arg(returnwhat)

  # Constants
  Boltz <- 5.67 * 10^-8     # w M-2 K-4
  Emissivity <- 0.95        # unitless
  LatEvap <- 2.54           # MJ kg-1
  CPAIR <- 1010             # J kg-1 K-1

  H2OLV0 <- 2501000         # J kg-1
  H2OMW <- 0.018            # J kg-1
  AIRMA <- 0.029            # mol mass air (kg/mol)
  AIRDENS <- 1.204          # kg m-3
  UMOLPERJ <- 4.57
  DHEAT <- 2.15e-05         # molecular diffusivity for heat

  # Get leaf VPD
  Dleaf = plantecophys::VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)

  # Convert temperatures to Kelvin
  Tair_k <- Tair + 273.15
  Tleaf_k <- Tleaf + 273.15

  # Density of dry air
  AIRDENS <- Patm * 1000/(287.058 * Tair_k)

  # Latent heat of water vapour at air temperature (J mol-1)
  LHV <- (H2OLV0 - 2365 * Tair) * H2OMW

  # Const s in Penman-Monteith equation  (Pa K-1)
  SLOPE <- (plantecophys::esat(Tair + 0.1) - plantecophys::esat(Tair))/0.1

  # Radiation conductance (mol m-2 s-1)
  Gradiation <- 4 * Boltz * Tair_k^3 * Emissivity/(CPAIR *
                                                     AIRMA)

  # See Leuning et al (1995) PC&E 18:1183-1200 Appendix E
  # Boundary layer conductance for heat - single sided, forced convection
  CMOLAR <- Patm * 1000/(8.314 * Tair_k)
  Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR

  # Free convection
  GRASHOF <- 1.6e+08 * abs(Tleaf - Tair) * (Wleaf^3) # Grashof number
  Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25)/Wleaf * CMOLAR

  # Total conductance to heat (both leaf sides)
  Gbh <- 2 * (Gbhfree + Gbhforced)

  # Rnet
  Rsol <- 2 * PPFD/UMOLPERJ # W m-2

  # Isothermal net radiation (Leuning et al. 1995, Appendix)
  ea <- plantecophys::esat(Tair) - 1000 * Dleaf
  ema <- 0.642 * (ea/Tair_k)^(1/7)
  Rnetiso <- LeafAbs * Rsol - (1 - ema) * Boltz * Tair_k^4 # isothermal net radiation

  GAMMA <- CPAIR * AIRMA * Patm * 1000/LHV

  # Convert transpiration to mol m-2 s-1
  ET <- E/1000

  # Latent heat loss
  lambdaET <- LHV * ET

  # Heat flux calculated using Gradiation (Leuning 1995, Eq. 11)
  Y <- 1/(1 + Gradiation/Gbh)
  H2 <- Y * (Rnetiso - lambdaET)

  # Heat flux calculated from leaf-air T difference.
  # (positive flux is heat loss from leaf)
  H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf)

  # Leaf-air temperature difference recalculated from energy balance.
  # (same equation as above!)
  Tleaf2 <- Tair + H2/(CPAIR * AIRDENS * (Gbh/CMOLAR))

  # Difference between input Tleaf and calculated, this will be minimized.
  EnergyBal <- Tleaf - Tleaf2
  if (returnwhat == "balance")
    return(EnergyBal)
  if (returnwhat == "fluxes") {
    l <- data.frame(ELEAFeb = 1000 * ET, Gradiation = Gradiation,
                    Rsol = Rsol, Rnetiso = Rnetiso, H = H, lambdaET = lambdaET,
                    Gbh = Gbh, H2 = H2, Tleaf2 = Tleaf2)
    return(l)
  }
}




#' Stomatal conductance
#'
#' @inheritParams calc_Tleaf
#' @param Tleaf Leaf temperature, deg C
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
#' Dleaves = plantecophys::VPDairToLeaf(Tleaf = Tleaves, Tair = 25, VPD = 1.5) # Get leaf VPD
#'
#' calc_gw(E = trans, Tleaf = Tleaves, VPD = Dleaves)
calc_gw = function (
    E,
    Tleaf = 25,
    Patm = 101.325,
    Tair = 25,
    VPD = 1.5,
    PPFD = 1000,
    Wind = 8,
    Wleaf = 0.025) {

  # Constants
  H2OLV0 <- 2.501e6         # J kg-1
  H2OMW <- 18e-3            # J kg-1
  CPAIR <- 1010             # J kg-1 K-1
  AIRMA <- 0.029            # mol mass air (kg/mol)
  UMOLPERJ <- 4.57
  DHEAT <- 2.15e-05         # molecular diffusivity for heat

  # Convert temperature to Kelvin
  Tair_k <- Tair + 273.15

  # Latent heat of water vapour at air temperature (J mol-1)
  LHV <- (H2OLV0 - 2.365E3 * Tair) * H2OMW

  # Psychrometric constant
  GAMMA <- CPAIR * AIRMA * Patm * 1000/LHV

  # Const s in Penman-Monteith equation  (Pa K-1)
  SLOPE <- (plantecophys::esat(Tair + 0.1) - plantecophys::esat(Tair))/0.1

  # See Leuning et al (1995) PC&E 18:1183-1200 Appendix E
  # Boundary layer conductance for heat - single sided, forced convection
  CMOLAR <- Patm * 1000/(8.314 * Tair_k)
  Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR

  # Free convection
  GRASHOF <- 1.6e+08 * abs(Tleaf - Tair) * (Wleaf^3) # Grashof number
  Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25)/Wleaf * CMOLAR

  # Total conductance to heat (both leaf sides)
  Gbh <- 2 * (Gbhfree + Gbhforced)

  # Rnet
  Rsol <- 2 * PPFD/UMOLPERJ # W m-2

  # Get leaf VPD
  Dleaf = plantecophys::VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)

  # diffusive equation
  g_w = E / 1000 * Patm / Dleaf

  # Penman-Monteith equation
  #g_w = GAMMA * Gbh * 1/((SLOPE * Rsol + Dleaf*1000 * Gbh * CPAIR * AIRMA)/(LHV * E/1000) - SLOPE)

  return(g_w)
}




#' Calculate photosynthesis
#'
#' @inheritParams calc_Tleaf
#' @param Ca Atmospheric CO2 concentration (ppm)
#' @param Jmax Maximum rate of electron transport at 25 deg C (mu mol m-2 s-1)
#' @param Vcmax Maximum carboxylation rate at 25 deg C (mu mol m-2 s-1)
#' @param net TRUE if desired output is net photosynthesis; FALSE if desired output
#'     is gross photosynthesis
#' @param Rd0 Day respiration rate at reference temperature (TrefR). Must be a positive value.
#' @param TrefR Reference temperature for Rd (deg C)
#' @param netOrig TRUE if net photosynthesis is to be calculated within
#'     \code{plantecophys::Photosyn}. FALSE if Anet is to be calculated as Agross - Rd.
#' @param g1 Parameter of Ball-Berry type stomatal conductance models.
#' @param g0 Parameter of Ball-Berry type stomatal conductance models.
#' @param g_w Stomatal conductance to water vapor (mol m-2 s-1)
#' @param Tleaf Leaf temperature (deg C)
#' @param \dots Further parameters passed to \code{plantecophys::Photosyn}.
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
calc_A = function(Tair = 25,
                  VPD = 1.5,
                  PPFD = 1000,
                  Patm = 101.325,
                  E = 2,
                  Wind = 2,
                  Wleaf = 0.025,
                  LeafAbs = 0.5,
                  Ca = 420,
                  Jmax = 100,
                  Vcmax = 50,
                  net = FALSE,
                  Rd0 = 0.92,
                  TrefR = 25,
                  netOrig = TRUE,
                  g1 = 2.9,
                  g0 = 0.003,
                  g_w = NULL,
                  Tleaf = NULL,
                  ...
)
{
  if (is.null(g_w) & is.null(Tleaf)) {
    Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD,
                       E = E, Wind = Wind, Patm = Patm, Wleaf = Wleaf,
                       LeafAbs = LeafAbs)
    g_w = calc_gw(E, Tleaf, Patm, Tair, VPD, PPFD, Wind,
                  Wleaf)
  }

  # Get leaf VPD
  Dleaf = plantecophys::VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)

  if (net == FALSE) {
    Photosyn_out = mapply(plantecophys::Photosyn, VPD = Dleaf,
                          Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm,
                          GS = g_w, Rd = 0, Jmax = Jmax, Vcmax = Vcmax, g1 = g1,
                          g0 = g0, ...)
    Anet = as.numeric(Photosyn_out[2, ])
    Ci <- as.numeric(Photosyn_out[1, ])
    Rd = as.numeric(Photosyn_out[8, ])
    A <- Anet + Rd
  }
  else {
    if (isTRUE(netOrig)) {
      Photosyn_out = mapply(plantecophys::Photosyn, VPD = Dleaf,
                            Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm,
                            GS = g_w, Jmax = Jmax, Vcmax = Vcmax, g1 = g1,
                            g0 = g0, ...)
      A <- as.numeric(Photosyn_out[2, ])
    }
    else {
      Rd = Rd0 * exp(0.1012 * (Tleaf - TrefR) - 5e-04 *
                       (Tleaf^2 - TrefR^2))
      Photosyn_out = mapply(plantecophys::Photosyn, VPD = Dleaf,
                            Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm,
                            GS = g_w, Rd = 0, Jmax = Jmax, Vcmax = Vcmax,
                            g1 = g1, g0 = g0, ...)
      A <- as.numeric(Photosyn_out[2, ]) - Rd
    }
  }
  out = list(A=A, Ci=Ci)
  return(out)
}


#' Radiative conductance
#'
#' @param Tair Air temperature, deg C
#' @param epsilon Leaf emissivity
#' @export
#' @return Radiative conductance, mol m-2 s-1
#' @examples
#' calc_gr(Tair = 30)
calc_gr = function(
    Tair = 25,
    epsilon = 0.97
)
{
  # Constants
  sigma = 5.67e-8 # Stefan-Boltzman constant, W m^-2 K^-4
  C_p = 29.3 # Specific heat of air at constant pressure, J mol^-1 deg C^-1

  TairK = Tair + 273 # air temperature in Kelvin
  g_r = max(0., 4 * epsilon * sigma * TairK ** 3 / C_p)
  return(g_r) # mol m-2 s-1
}



#' Boundary layer conductance to forced convection
#'
#' @param Tair Air temperature, deg C
#' @param Wind Wind speed above the leaf boundary layer, m s-1
#' @param Patm Atmospheric pressure, kPa
#' @param Wleaf Leaf width, m
#'
#' @export
#'
#' @return Boundary layer conductance to forced convection, mol m-2 s-1
calc_gHa = function(
    Tair = 25,
    Wind = 2,
    Patm = 101.325,
    Wleaf = 0.025
)
{
  # Constants
  Mair = 28.9644  # molar mass of dry air (g mol-1)
  DH = 21.5e-6  # molecular diffusivity to heat (m2 s-1)
  R = 8.3144598  # molar ideal gas constant (J mol-1 K-1)

  TairK = Tair + 273 # air temperature in Kelvin
  cmolar = Patm * 1e3 / (R * TairK)  # air molar density

  # Sutherland Eq for dynamic viscosity
  mu = 1.458e-6 * TairK ** 1.5 / (TairK + 110.4)  # Pa s

  # Kinematic viscosity
  nu = mu * R * TairK / (Patm * Mair)  # m2 s-1
  prandtl = nu / DH  # unitless

  # Boundary layer cond to forced convect. (Eq. 7.29 Campbell & Norman, 1998)
  d = 0.72 * Wleaf  # leaf width, m
  reynolds = Wind * d / nu  # unitless
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
calc_Rd = function(Tair = 25,
                   VPD = 1.5,
                   PPFD = 1000,
                   E = 2,
                   Wind = 2,
                   Patm = 101.325,
                   Wleaf = 0.025,
                   LeafAbs = 0.5,
                   Rd0 = 0.92,
                   TrefR = 25
                   )
{
  Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, E = E, Wind = Wind,
                      Patm = Patm, Wleaf = Wleaf, LeafAbs = LeafAbs)
  Rd = Rd0 * exp(0.1012 * (Tleaf - TrefR) - 5e-04 * (Tleaf^2 - TrefR^2))
  return(Rd)
}
