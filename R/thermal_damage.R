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
