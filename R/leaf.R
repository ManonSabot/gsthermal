#' Leaf temperature
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
#' Tleaves = calc_Tleaf(E = trans) # Get leaf temperatures
calc_Tleaf = function(
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
#'
#' @return Photosynthetic rate
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
                  RH = 60
)
{
  T_leaf = calc_Tleaf(T_air, PPFD, RH, E, u, Patm, leaf_width)
  D_leaf = calc_Dleaf(T_leaf, T_air, RH)
  g_w = calc_gw(E, D_leaf, Patm)
  Photosyn_out = mapply(plantecophys::Photosyn,
                        VPD = D_leaf, Tleaf = T_leaf, GS = g_w, Rd0 = 0)
  A = as.numeric(Photosyn_out[2,])
  return(A)
}

#' Radiative conductance
#'
#' @param T_air Air temperature, deg C
#' @param epsilon Leaf emissivity
#' @noRd
#' @return Radiative conductance, mol m-2 s-1
#' @example
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
#' @noRd
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
