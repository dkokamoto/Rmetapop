#' Apply a Kalman filter to a time series 
kalman_assess <- function(TS, sd1, sd2, fit = FALSE, full.ts = FALSE, log = TRUE) {
    
    if (fit == TRUE) {
        TSBuild <- function(par) {
            dlmModPoly(1, dV = par[1], dW = par[2])
        }
        TSMLE <- dlmMLE(TS, rep(0, 2), TSBuild)
        TSMod <- TSBuild(TSMLE$par)
    } else {
        TSBuild <- function(sd1, sd2) {
            dlmModPoly(1, dV = sd1, dW = sd2)
        }
        TSMod <- TSBuild(sd1 = sd1, sd2 = sd2)
    }
    TSFilt <- dlmFilter(TS, TSMod)
    if (full.ts == TRUE) {
        return(exp(TSFilt$m[-1] - 0.5 * sd1^2))
    } else {
        return(exp(TSFilt$m[length(TSFilt$m)] - 0.5 * sd1^2))
    }
} 
