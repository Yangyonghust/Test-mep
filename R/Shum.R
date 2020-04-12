#'Calculate specific humidity for MEP model input
#'@param TA air temperature(unit: deg C)
#'@param RH Relative humidity(unit: %)
#'@param PA Atmospheric pressure(unit: kPa)
#'@return Specific humidity(unit:kg/kg)
#'@export
Shum <- function(TA, PA, RH){
  Tk <- C2K(TA)
  Es <- SVP(Tk)
  E <- WVP2(RH,Es)
  qs <- SH(E,PA*1000)   #specific humidity in kg/kg
  return(qs)
}
