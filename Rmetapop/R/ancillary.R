#' Beverton-Holt stock recruit relationship.
#' @param E Number of eggs. 
#' @param E0 number of eggs produced at B0.
#' @param h steepness.
#' @param R0 recruits produced at E0.
#' @param alpha =E0 * (1 - h))/(4 * h * R0)
#' @param beta =(5 * h - 1)/(4 * h * R0)
#' @description Generates the expected number of recruits given egg production and a set of parameters.  The parameters must either be E0, h, and R0 OR alpha and beta. Uses the following form: 
#' \deqn{R = E/(alpha+beta*E)}
#'\deqn{alpha = [E0*(1-h)]/(4*h*R0)}
#' \deqn{beta = (5*h-1)/(4*h*R0)}
#' \deqn{h = BH(0.2*E0)/BH(E0)}
#' @example /inst/examples/BH_example.R
BH <- function(E, E0=NULL, R0=NULL,h=NULL,  alpha = NULL, beta = NULL) {
    if (!is.null(h)) {
        alpha_h <- (E0 * (1 - h))/(4 * h * R0)
        beta_h <- (5 * h - 1)/(4 * h * R0)
        E/(alpha_h + beta_h * E)
    } else {
        E/(alpha + beta * E)
    }
}

### Allometric functions
#' Ludwig Von-Bertalanfy growth equation.
LVB <- function(age, L_inf = 27, k = 0.48, t0 = 0) {
    L_inf * (1 - exp(-k * (age - t0)))
}

#' Power length-weight relationship.
weight <- function(length, a = 4.5e-06, b = 3.127) {
    a * (length)^b
}

#' Weight at age using the LVB growth equation and length-weight relationship.
LVBweight <- function(age) {
    weight(LVB(age))
}

#' Fecundity at length function.
fecundity <- function(L, F1 = 0.000419, F2 = 3.372) {
    F1 * L^(F2) * 10^F2  ### uses length in cm not mm
}

#' Fecundity at age using LVB growth equation and fecundity at length function
fecundity_age <- function(age) {
    fecundity(LVB(age))
}

#' Define a system of ordinary differential equations from a matrix
linear_odes <- function(t, state, A) {
    dX <- A %*% state
    list(as.vector(dX))
} 
