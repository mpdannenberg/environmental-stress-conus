# Stress Index Functions

## Based on diameter
StressIndex <- function(rwl, pars, sims=FALSE, p0=NULL, model=NULL){
  G = pars$G
  Dmax = pars$Dmax
  Hmax = pars$Hmax
  
  if(is.na(pars$b2)) b2 = 2 * ((Hmax-137)/Dmax)
  else b2 = pars$b2
  
  if(is.na(pars$b3)) b3 = (Hmax-137)/Dmax^2
  else b3 = pars$b3
  
  inds = which(!is.na(rwl))
  S <- rep(NA, length(rwl))
  #Dopt_out <- rep(NA, length(rwl))
  
  if(sims==FALSE){
		D = 2 * cumsum(rwl[inds]*0.1)
	} else{
		pith = rbinom(1,1,p0)
		if(pith==1){
			D = 2 * cumsum(rwl[inds]*0.1)
		} else{
			bias = rgamma(1, shape=model$estimate[1], rate=model$estimate[2])
			D = 2 * cumsum(rwl[inds]*0.1) + bias
		}
	}
  H = 137 + b2*D - b3*D^2
  Dopt = G*D*(1 - (D*H)/(Dmax*Hmax)) / (274 + 3*b2*D - 4*b3*D^2)
  S[inds] = (2*rwl[inds]*0.1) / Dopt
  #Dopt_out[inds] = Dopt
  return(S)
}


## Based on basal area
StressIndex2 <- function(rwl, pars, sims=FALSE, p0=NULL, model=NULL){
	# rwl: ring widths (mm)
  G = pars$G
  Dmax = pars$Dmax
  Hmax = pars$Hmax
  
  if(is.na(pars$b2)) b2 = 2 * ((Hmax-137)/Dmax)
  else b2 = pars$b2
  
  if(is.na(pars$b3)) b3 = (Hmax-137)/Dmax^2
  else b3 = pars$b3
  
  inds = which(!is.na(rwl))
  S <- rep(NA, length(rwl))
  Dopt_out <- rep(NA, length(rwl))
  
  if(sims==FALSE){
		D = 2 * cumsum(rwl[inds]*0.1)
	} else{
		pith = rbinom(1,1,p0)
		if(pith==1){
			D = 2 * cumsum(rwl[inds]*0.1)
		} else{
			bias = rgamma(1, shape=model$estimate[1], rate=model$estimate[2])
			D = 2 * cumsum(rwl[inds]*0.1) + bias
		}
	}
  H = 137 + b2*D - b3*D^2
  Dopt = G*D*(1 - (D*H)/(Dmax*Hmax)) / (274 + 3*b2*D - 4*b3*D^2)
  
  bai_opt = pi*(Dopt/2 + D/2)^2 - pi*(D/2)^2
  bai_act = pi*(rwl[inds]*0.1 + D/2)^2 - pi*(D/2)^2
  
  S[inds] = bai_act / bai_opt
  Dopt_out[inds] = Dopt
  return(list(S, Dopt_out))
}


D_bias <- function(rwl, model){
	pith = rbinom(1,1,p0)
	if(pith==1){
		bias = 0
	} else{
		bias = rgamma(1, shape=model$estimate[1], rate=model$estimate[2])
	}
	return(bias)
}





## Output Dopt, instead of Stress Index
CalcDopt <- function(rwl, pars, sims=FALSE, p0=NULL, model=NULL){
  G = pars$G
  Dmax = pars$Dmax
  Hmax = pars$Hmax
  
  if(is.na(pars$b2)) b2 = 2 * ((Hmax-137)/Dmax)
  else b2 = pars$b2
  
  if(is.na(pars$b3)) b3 = (Hmax-137)/Dmax^2
  else b3 = pars$b3
  
  inds = which(!is.na(rwl))
  #S <- rep(NA, length(rwl))
  Dopt_out <- rep(NA, length(rwl))
  
  if(sims==FALSE){
		D = 2 * cumsum(rwl[inds]*0.1)
	} else{
		pith = rbinom(1,1,p0)
		if(pith==1){
			D = 2 * cumsum(rwl[inds]*0.1)
		} else{
			bias = rgamma(1, shape=model$estimate[1], rate=model$estimate[2])
			D = 2 * cumsum(rwl[inds]*0.1) + bias
		}
	}
  H = 137 + b2*D - b3*D^2
  Dopt = G*D*(1 - (D*H)/(Dmax*Hmax)) / (274 + 3*b2*D - 4*b3*D^2)
  #S[inds] = (2*rwl[inds]*0.1) / Dopt
  Dopt_out[inds] = Dopt
  return(Dopt_out)
}


TrimAge <- function(rwl, SetAge){
  inds = which(!is.na(rwl))
  n = length(inds)
  S <- rep(NA, length(rwl))
  
  if(n > SetAge){
    inds = inds[(SetAge+1):n]
    S[inds] = rwl[inds]
  } 
  
  return(S)
}

TrimLM <- function(rwl){
  # Running linear models to find significant trends
}


