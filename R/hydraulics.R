#' Fit Weibull
#' @description Fit Weibull parameters given values of P50 and P88
#'
#' @param P50 LWP at 50 PLC (-MPa)
#' @param P88 LWP at 88 PLC (-MPa)
#'
#' @return A data frame containing the Weibull parameters b and c
#' @export
fit_Weibull = function(
    P50 = 4.9,
    P88 = 6.41
)
{
  x1 = 0.50
  x2 = 0.88

  c = log(log(1. - x1) / log(1. - x2)) / (log(P50) - log(P88))
  b = - P50 / ((- log(1 - x1)) ** (1. / c))

  return(data.frame(b, c))
}




#' Solve for Pcrit
#' @description Calculate the critical leaf water potential
#'
#' @param b Weibull scale parameter
#' @param c Weibull shape parameter
#' @param ratiocrit Percentage of maximum conductivity at which hydraulic damage
#'     is considered irreversible
#'
#' @return The critical leaf water potential (-MPa)
#' @export
#'
#' @examples
#' Weibull = fit_Weibull()
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#'
#' Pcrit = calc_Pcrit(b, c)
#'
#' # Change the critical ration
#' Pcrit1 = calc_Pcrit(b, c, 0.025)
calc_Pcrit = function(b = -2.5,
                      c = 2,
                      ratiocrit = 0.05)
{
  Pcrit = -b*(-log(ratiocrit))**(1/c)
  return(Pcrit)
}




#' Leaf water potential vector
#' @description Create vector of equally-spaced steps ranging from soil water
#'     potential to the critical water potential
#'
#' @param Ps Soil water potential, -MPa
#' @param Pcrit Critical leaf water potential, -MPa
#' @param pts Number of steps
#'
#' @return A vector of equally-spaced leaf water potentials ranging from the
#'    soil water potential to the critical water potential
#' @export
#'
#' @examples
#' # Calculate Pcrit based on Weibull curve
#' Weibull = fit_Weibull()
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c)
#'
#' # Create Ps to Pcrit vector
#' P = Ps_to_Pcrit(Pcrit = Pcrit)
Ps_to_Pcrit = function(Ps = 0, # -MPa
                       Pcrit = 4, # -MPa
                       pts = 500 # no. of steps
                       )
{seq(Ps, Pcrit, length.out = pts)}




#' Vulnerability curve
#' @description Calculate the vulnerability curve for a given soil and critical
#'     water potential
#'
#' @param P Vector with water potentials, -MPa
#' @param b Weibull scale parameter
#' @param c Weibull shape parameter
#'
#' @return Vulnerability curve
#' @export
#'
#' @examples
#' # Calculate Pcrit based on Weibull curve
#' Weibull = fit_Weibull()
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c)
#'
#' # Create Ps to Pcrit vector
#' P = Ps_to_Pcrit(Pcrit = Pcrit)
#'
#' # Fit vulnerability curve
#' vulnerability_curve(P, b, c)
vulnerability_curve = function(P,
                               b = -2.5,
                               c = 2)
{
  curve = (exp(-(-P/b) ** c))
  return(curve)
}




#' Temperature dependence of conductivity
#' @description Calculate whole plant conductivity as a function of air
#'     temperature
#'
#' @param kmax_25 Max plant conductance at 25 deg C, mmol s-1 m-2 MPa-1
#' @param Tair Air temperature, deg C
#' @param constant_kmax TRUE if the kmax does not vary with temperature for
#'     simulations; else FALSE
#'
#' @return Whole plant conductivity at Tair, mmol s-1 m-2 MPa-1
#' @export
calc_kmax = function(kmax_25 = 0.5,
                     Tair = 25,
                     constant_kmax = FALSE
)
{
  if (isFALSE(constant_kmax)) {
  TairK = Tair+273.15
  kmax = kmax_25 * (TairK**7 / 298.15**7)
  } else {
    kmax = kmax_25
  }
  return(kmax)
}




#' Transpiration supply function
#' @description Calculate transpiration as integral of a vulnerability curve
#'
#' @param P Vector with water potentials, -MPa
#' @param kmax_25 Max plant conductance at 25 deg C, mmol s-1 m-2 MPa-1
#' @param Tair Air temperature, deg C
#' @param b Weibull scale parameter
#' @param c Weibull shape parameter
#' @param constant_kmax TRUE if the kmax does not vary with temperature for
#'     simulations; else FALSE
#'
#' @return Transpiration, mmol s-1 m-2
#' @export
#'
#' @examples
#' # Calculate Pcrit based on Weibull curve
#' Weibull = fit_Weibull()
#' b = Weibull[1,1]
#' c = Weibull[1,2]
#' Pcrit = calc_Pcrit(b, c)
#'
#' # Create Ps to Pcrit vector
#' P = Ps_to_Pcrit(Pcrit = Pcrit)
#'
#' # Create vector of transpiration supply stream
#' trans = trans_from_vc(P = P, b = b, c = c)
trans_from_vc = function(P,
                         kmax_25 = 0.5,
                         Tair = 25,
                         b = -2.5,
                         c = 2,
                         constant_kmax = FALSE
                         )
{
  VC = vulnerability_curve(P, b, c)
  kmax = calc_kmax(kmax_25, Tair, constant_kmax)

  # Create lists of first n water potential and PLC values for all n from 1 to the number of points in P
  P_list = sapply(1:length(P), function(n) P[1:n])
  #VC_list = sapply(1:length(VC), function(n) VC[1:n])

  # Approximate integral as a trapezoidal sum
  #AUC = mapply(pracma::trapz, P_list, VC_list)
  AUC = (b / c) * (expint::gammainc(a=1 / c, x=(-P[1] / b)^c) - expint::gammainc(a=1 / c, x=(-P / b)^c))

  E = kmax*AUC
  return(E)
}
