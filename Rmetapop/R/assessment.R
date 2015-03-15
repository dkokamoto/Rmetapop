#' Apply a Kalman filter to a time series 
kalman_assess <- function(TS, obs_sd, sys_sd, fit = FALSE, full.ts = FALSE, log = TRUE) {
    
    if (fit == TRUE) {
        TSBuild <- function(par) {
            dlmModPoly(1, dV = par[1], dW = par[2])
        }
        TSMLE <- dlmMLE(TS, rep(0, 2), TSBuild)
        TSMod <- TSBuild(TSMLE$par)
    } else {
        TSBuild <- function(obs_sd, sys_sd) {
            dlmModPoly(1, dV = obs_sd, dW = sys_sd)
        }
        TSMod <- TSBuild(obs_sd=obs_sd, sys_sd=sys_sd)
    }
    TSFilt <- dlmFilter(TS, TSMod)
    if (full.ts == TRUE) {
        return(exp(TSFilt$m[-1] - 0.5 * obs_sd^2))
    } else {
        return(exp(TSFilt$m[length(TSFilt$m)] - 0.5 * obs_sd^2))
    }
} 
