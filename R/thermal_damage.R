#' F0-T Curve
#' @description Calculates the value of minimum fluorescence (F0) associated
#'     with a given leaf temperature. F0 values serve as a proxy of thermal
#'     damage to Photsystem II.
#'
#' @param Tleaf Leaf temperature, deg C
#' @param Tcrit Leaf temperature at which F0(T) transitions from slow- to fast- rise, deg C
#' @param T50 Leaf temperature midway between Tcrit and Tmax, deg C
#' @param F0_max Maximum F0 (achieved at Tmax)
#' @param F0_min Minimum (i.e. undamaged) F0
#'
#' @return Minimum fluorescence (F0)
#' @export
#'
#' @examples
#' Tleaf = seq(25, 50, length.out = 500)
#' curve = F0_func(Tleaf = Tleaf)
#' plot(Tleaf, curve, type = "l")
F0_func = function(Tleaf,
                   Tcrit = 50,
                   T50  = 51,
                   F0_max = 1000,
                   F0_min = 500
)
{
  r = 2 / (T50 - Tcrit)

  F0 = (F0_max - F0_min) / (1 + exp(-r * (Tleaf - T50))) + F0_min
  return(F0)
}

#' Modified Jmax-T response
#' @description
#' Modified temperature response of Jmax. Based on the peaked Arrhenius temperature
#' response, but multiplied by (1 - TC) so that Jmax decreases as the TC increases.
#'
#' @param Tleaf Leaf temperature (deg C)
#' @param EaJ,EdVJ,delsJ Jmax temperature response parameters
#' @param Tcrit Leaf temperature at which F0(T) transitions from slow- to fast- rise, deg C
#' @param T50 Leaf temperature midway between Tcrit and Tmax, deg C
#'
#' @return Jmax (normalised to Jmax at 25 deg C)
#' @export
TJmax_updated <- function(Tleaf,
                          EaJ,
                          delsJ,
                          EdVJ,
                          Tcrit = 43.4,
                          T50 = 49.6){
  J1 <- 1+exp((298.15*delsJ-EdVJ)/8.314/298.15)
  J2 <- 1+exp(((Tleaf + 273.15)*delsJ-EdVJ)/8.314/(Tleaf + 273.15))

  r = 2/(T50 - Tcrit)
  functionality = 1 - 1/(1 + exp(-r * (Tleaf - T50)))

  J = functionality * (exp(EaJ/8.314*(1/298.15 - 1/(Tleaf + 273.15)))*J1/J2)

  return(J)
}
