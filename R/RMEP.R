#' Calculate energy fluxes and actual evapotranspiration
#'
#' Calculate energy fluxes and evapotranspiration based on the maximum entropy production(MEP) model
#' @param Rn net radiation(unit:W/m2)
#' @param RnL net long-wave radiation(unit:W/m2)
#' @param Ts surface temperature(unit:Celsius)
#' @param qs specific humidity(kg/kg)
#' @return A list includes latent heat flux(EMEP),sensible heat flux(HMEP),ground heat flux(GMEP) and evapotranspiration(ETMEP)
#' @export
RMEP<-function(Rn, RnL, qs, Ts, type)   UseMethod ("RMEP")
RMEP.default=function(Rn, RnL, qs, Ts, rhoa=1.18, I=600, z=2.5, type=1){

  if(length(Rn)>1){
    rhoa=Rn * 0 + 1.18; I=Rn * 0 + 600.0; z=Rn * 0 + 2.5 #type=Rn*0+1;
  }
# Parameters
  T0 = Rn*0 + 273.15                              # freezing point
  Ts = Ts + T0                                    # Ts in K
  Tr = Rn*0 + 300                                 # reference temperature
  Cp = Rn*0 + 1006                                # specific heat of air
  Rv = Rn*0 + 465                                 # gas constant of water vapor
  a = Rn*0 + 1                                    # coefficeint of water version
  g = Rn*0 + 9.81                                 # gravitational accelaration
  k = Rn*0 + 0.4                                  # Von Karman constant
  Lv = Rn*0 + 2.5E06                              # vaporization heat of water vapor
  Liv = Rn*0 + 2.83E06                            # vaporization heat of ice
  alpha =Rn*0 + 1                                 # MOST constant
  beta = Rn*0 + 4.7                               # MOST constant
  r1 = Rn*0 + 15                                  # MOST constant
  r2 = Rn*0 + 9                                   # MOST constant
  C1 = matrix(c(sqrt(3)/alpha, 2/(1 + 2 * alpha)),ncol=2,byrow=FALSE)    # MEP constant
  C2 = matrix(c(r2/2, 2 * beta),ncol=2,byrow=FALSE)                    # MEP constant
  e16 = 1.0/6                                     # exponent

# Calculate Variables
  stab = ones(1,length(Rn))                       # Create vector of stability
  stab[which(Rn < 0)] = 0                         # Calculate Variables
  I0 = zeros(1,length(Rn))                        # Create vector of I0
  unsta = which(stab ==1)                         # Unstable condition
  sta = which(stab ==0)                           # Stable condition

  if (1-isempty(unsta)){
    I0[unsta] = rhoa[unsta] * Cp * sqrt(C1[,1] * k * z) * (C2[,1] * k * z * g/(rhoa[unsta] * Cp * Tr))^(e16);
  }

  if (1-isempty(sta)){
    I0[sta] = rhoa[sta] * Cp * sqrt(C1[,1] * k * z) * (C2[,1] * k * z * g/(rhoa[unsta] * Cp * Tr))^(e16);
  }
  ice_num = which(Ts < T0);                      # ice/snow condition (only used in water-snow-ice version type = 3)
  Lv_vec = ones(1,length(Ts)) * Lv;              # create Lv vector
  Lv_vec[ice_num] = Liv;
  sg = sqrt(a) * Lv_vec^2 * qs/(Cp * Rv * Ts^2);       # sigma computation
  B = 6 * (sqrt(1 + 11/36 * sg) - 1);                # B computation
  Bsr = B / sg;                                    # B/sigma
  IdI0 = I / I0;                                   # I/I0 Stable condition
  # grab effective data (by skipping NaN values)
  id_1 = which(is.na(Rn) == 0)
  id_2 = which(is.na(RnL) == 0)
  id_3 = which(is.na(qs) == 0)
  id_4 = which(is.na(Ts) == 0)
  id_5 = which(is.na(rhoa) == 0)
  id_6 = which(is.na(I) == 0)
  id_7= which(is.infinite(IdI0) == 0)
  id_8 = which(is.na(IdI0) == 0)

  ef_id = intersect(id_1, id_2)
  ef_id = intersect(ef_id, id_3)
  ef_id = intersect(ef_id, id_4)
  ef_id = intersect(ef_id, id_5)
  ef_id = intersect(ef_id, id_6)
  ef_id = intersect(ef_id, id_7)
  ef_id = intersect(ef_id, id_8)

  # MEP calculation
  HMEP = zeros(1,length(ef_id))
  EMEP = zeros(1,length(ef_id))
  GMEP = zeros(1,length(ef_id))

  if (type == 1){           # Bare Soil
    for (n in c(1:length(ef_id))) {
      temp = fzero(function(H) abs(H)^(e16) * (B[ef_id[n]]+1) * H +
                     Bsr[ef_id[n]] * IdI0[ef_id[n]] * H - abs(H)^(e16) * Rn[ef_id[n]],0.5 * Rn[ef_id[n]])[1];
      HMEP[ef_id[n]]=temp[[1]];
      EMEP[ef_id[n]] = B[ef_id[n]] * HMEP[[ef_id[n]]];
      GMEP[ef_id[n]] = Rn[ef_id[n]] - EMEP[ef_id[n]] - HMEP[[ef_id[n]]];
    }
  } else if (type == 2) {   # Canopy
    for (n in c(1:length(ef_id))) {
      temp = fzero(function(H) (B[ef_id[n]]+1) * H -
                     Rn[ef_id[n]],0.5 * Rn[ef_id[n]])[1];
      HMEP[ef_id[n]]=temp[[1]];
      EMEP[ef_id[n]]=Rn[ef_id[n]] - HMEP[[ef_id[n]]];
      GMEP[ef_id[n]]=0;
    }
  } else {                  # Water-snow-ice (type == 3)
    for (n in c(1:length(ef_id))) {
      temp = fzero(function(H) abs(H)^(e16) * (B[ef_id[n]]+1) * H +
                     Bsr[ef_id[n]] * IdI0[ef_id[n]] * H - abs(H)^(e16) * Rn[ef_id[n]],0.5 * Rn[ef_id[n]])[1];
      HMEP[ef_id[n]]=temp[[1]];
      EMEP[ef_id[n]] = B[ef_id[n]] * HMEP[ef_id[n]];
      GMEP[ef_id[n]]=RnL[ef_id[n]] - HMEP[ef_id[n]] - EMEP[ef_id[n]];
    }
  }
  return(list(ETMEP = EMEP * 0.0352653,HMEP = HMEP,EMEP = EMEP,GMEP = GMEP,ef_id = ef_id))   #1Wm-2=0.0864mmd-1
}


#'Calculate potential ET
#'
#'Calculate energy fluxes and potential evapotranspiration based on the MEP model
#' @param Rn net radiation(unit:W/m2)
#' @param RnL net long-wave radiation(unit:W/m2)
#' @param Ts surface temperature(unit:Celsius)
#' @return A list concludes latent heat flux(EMEP),sensible heat flux(HMEP),ground heat flux(GMEP) and potential ET(PETMEP)
#' @export
RMEPPET<-function(Rn, RnL, Ts, type)   UseMethod ("RMEPPET")
library(humidity)
RMEPPET.default=function(Rn, RnL,Ts,pres=101325, rhoa=1.18, I=600, z=2.5, type=1){
  t_temp<-C2K(Ts)
  Es<-SVP.ClaCla(t_temp)        #calculating saturation vapor pressure Es
  SShums<-SH(Es * 100,p=pres)  #calculating saturated specific humidity
  output<-RMEP(Rn=Rn,RnL=RnL,qs=SShums,Ts=Ts)
  print("RMEPPET completed!")
  PETMEP<-output$ETMEP;HMEP<-output$HMEP;EMEP<-output$EMEP;GMEP<-output$GMEP;ef_id=output$ef_id;
  return(list(PETMEP = PETMEP,HMEP = HMEP,EMEP = EMEP,GMEP = GMEP,ef_id = ef_id))
}




