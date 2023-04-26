#' Vapor pressure deficit
#' @description Calculate the atmospheric VPD at given relative humidity and air
#'     temperature
#'
#' @param T_air Air temperature, deg C
#' @param RH Relative humidity, in \%
#'
#' @return Vapor pressure deficit, kPa
#' @export
calc_VPD = function(
    T_air = 25, # deg C
    RH = 60 # unitless
)
{
  esat_a = vpsat(T_air)
  VPD = esat_a * (1 - RH / 100) # kPa
  return(VPD)
}




#' Calculate net radiation
#'
#' @param T_air Air temperature, deg C
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param VPD Atmospheric vapor pressure deficit, kPa
#' @param albedo Leaf albedo, unitless
#' @param epsilon Leaf emissivity, unitless
#'
#' @return Net radiation, W m-2
#' @export
calc_Rnet = function(T_air = 25,
                     PPFD = 1000,
                     VPD = 1.5,
                     albedo = 0.15,
                     epsilon = 0.97
)
{
  # Constants
  SW_2_PAR = 4.57 * 0.5  # SW (W m-2) to PAR (mumol m-2 s-1)
  PAR_2_SW = 1 / SW_2_PAR
  sigma = 5.67e-8 # Stefan-Boltzman constant, W m-2 K-4

  # unit conversions
  TairK = T_air + 273 # air temperature in Kelvin
  emissivity = calc_emissivity(T_air, VPD)
  # incoming short and long wave radiation
  Rsw = (1. - albedo) * PPFD * PAR_2_SW  # W m-2
  Rlw = emissivity * sigma * TairK ** 4.  # W m-2

  Rnet = Rsw + Rlw - epsilon * sigma * TairK ** 4 # W m-2
  return(Rnet)
}

#' Saturation vapor pressure at a given temperature
#'
#' @param T Temperature, deg C
#' @return Saturation vapor pressure, kPa
#' @noRd
vpsat = function(
    T = 25
)
{
  # Constants
  a = 0.61078 # kPa
  b = 17.27 # unitless
  c = 237.3 # deg C

  vpsat = a*exp(b*T/(T+c)) #kPa
  return(vpsat)
}

#' Calculate atmospheric emissivity
#'
#' @param T_air Air temperature, deg C
#' @param VPD Atmospheric vapor pressure deficit, kPa
#'
#' @return Apparent atmospheric emissivity, unitless
#' @noRd
calc_emissivity = function(T_air = 25,
                           VPD = 1.5
)
{
  # Constants
  sigma = 5.67e-8 # Stefan-Boltzman constant, W m-2 K-4

  TairK = T_air + 273 # air temperature in Kelvin
  ea = (vpsat(T_air) - 4)*1e3
  emissivity = (0.031 * ea + 2.84 * TairK - 522.5) / (sigma * TairK ** 4.)
  return(emissivity)
}
