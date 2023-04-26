#' F0-T Curve
#' @description Calculates the value of minimum fluorescence (F0) associated
#'     with a given leaf temperature. F0 values serve as a proxy of thermal
#'     damage to Photsystem II.
#'
#' @param T_leaf Leaf temperature, deg C
#' @param T50 Leaf temperature midway between Tcrit and Tmax, deg C
#' @param F0_max Maximum F0 (achieved at Tmax)
#' @param F0_min Minimum (i.e. undamaged) F0
#' @param r Scaling factor
#'
#' @return Minimum fluorescence (F0)
#' @export
#'
#' @examples
#' Tleaf = seq(25, 50, length.out = 500)
#' curve = F0_func(T_leaf = Tleaf)
#' plot(Tleaf, curve, type = "l")
F0_func = function(T_leaf,
                   T50  = 51,
                   F0_max = 1000,
                   F0_min = 500,
                   r = 4 #
)
{
  F0 = (F0_max - F0_min) / (1 + exp(-r * (T_leaf - T50))) + F0_min
  return(F0)
}
