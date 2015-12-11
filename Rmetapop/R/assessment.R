#' Apply a Kalman filter to a time series using a linear DLM 
kalman_assess <- function(TS, obs_sd, sys_sd, fit = FALSE, full.ts = FALSE, log = TRUE) {
    
    if (fit == TRUE) {
        TSBuild <- function(par) {
          dlmModPoly(1, dV = obs_sd^2, dW = exp(par))
        }
        TSMLE <- dlmMLE(log(TS), log(obs_sd^2), TSBuild)
        TSMod <- TSBuild(TSMLE$par)
    } else {
        TSBuild <- function(obs_sd, sys_sd) {
            dlmModPoly(1, dV = obs_sd^2, dW = sys_sd^2)
        }
        TSMod <- TSBuild(obs_sd=obs_sd^2, sys_sd=sys_sd^2)
    }
    TSFilt <- dlmFilter(TS, TSMod)
    filtered_ts <- exp(TSFilt$m[-1] + 0.5 * obs_sd^2)                             
    if (full.ts == TRUE) {
        return(filtered_ts)
    } else {
        return(filtered_ts[length(filtered_ts)])
    }
} 

kalman_assess_ARMA <- function(TS, obs_sd, sys_sd, fit = FALSE,
                               full.ts = FALSE, log = TRUE) {
  
  if (fit == TRUE) {
    TSBuild <- function(par) {
      dlmModARMA(ma= par[1],dV = obs_sd^2, sigma2= exp(par[2]))
    }
    TSMLE <- dlmMLE(log(B[1, , 2]), obs_sd^2, TSBuild)
    TSMod <- TSBuild(TSMLE$par)
  } else {
    TSBuild <- function(obs_sd, sys_sd) {
      dlmModPoly(1, dV = obs_sd^2, dW = sys_sd^2)
    }
    TSMod <- TSBuild(obs_sd=obs_sd^2, sys_sd=sys_sd^2)
  }
  TSFilt <- dlmFilter(TS, TSMod)
  system_variance <- ifelse(is.null(TSMLE$par),"NA",exp(TSMLE$par))
  filtered_ts <- exp(TSFilt$m[-1] + 0.5 * obs_sd^2)
  returned <- list(filt=filtered_ts,sys_var=system_variance)                                 
  if (full.ts == TRUE) {
    return(returned)
  } else {
    return(returned)
  }
} 



