#' Discrete survival, straying, and recruitment projection.
#' @param stray_mat a square L x L stray matrix. 
#' @param mort_mat an S x L mortality matrix OR a vector of length L OR a single value.
#' @param fec_at_age a vector or matrix of fecundity at age including only mature stages.
#' @param eggs a vector of eggs associated with the next years cohort (fed into the stock recruit relationship).
#' @param E0,h,R0  Beverton-Holt parameters (E0= egg production at B0, h= steepness, R0 = recruitment at B0).
#' @param alpha,beta Direct Beverton-Holt parameters (alpha = per capita recruitment as egg production approaches 0; beta = inverse of recruitment asymptote).
#' @param stage_maturity stage class at maturity.
#' @param errors a vector of length L or single value of log-scale recruitment deviations.
#' @param X an S x L matrix of initial abundances of each stage (S) at each location (L).
#' @param n_stages number of stage classes.
#' @param n_loc number of locations. 
#' @example /inst/examples/survexample.R
#' @description Generates a single projection including survival followed by straying and reproduction using a matrix algebra in discrete time.
#' @seealso surv_stray_recr_ode

SSR_linear <- function(stray_mat, mort_mat, fec_at_age, 
                            eggs, E0, h , R0, 
                            alpha=NULL,beta=NULL,
                            stage_maturity, errors = 0, X0, n_stages, n_loc) {
    
    ### create empty array
    X1 <- array(as.numeric(NA), dim = c(n_stages, n_loc))
    
    if (length(dim(mort_mat)) > 2) {
        ### fill in mortality at age
        for (j in 1:5) {
            X1[(stage_maturity - 1):n_stages, j] <- mort_mat[, , j] %*% X0[, j]
        }
    } else {
        X1[(stage_maturity - 1):n_stages, ] <- mort_mat %*% X0
    }
    X1[(stage_maturity - 1):n_stages, ] <- X1[2:n_stages, ] %*% t(stray_mat)
    X1[1, ] <- fec_at_age %*% X1[stage_maturity:n_stages, ]
    
    # add recruitment 
    if (is.null(alpha)){
      X1[2, ] <- mapply(BH, E = eggs, E0 = E0, h = h, R0 = R0) * exp(errors)
    }
    else {
      X1[2, ] <- mapply(BH, E = eggs, alpha= alpha,beta=beta) * exp(errors)
    }
    return(X1)
}

#' Continuous time survival, straying, and harvest with discrete recruitment
#' @param stray a square L x L stray matrix. 
#' @param Z an S x L instantaneous mortality matrix OR a vector of length L OR a single value.
#' @param fec_at_age a vector or matrix of fecundity at age including only mature stages.
#' @param eggs a vector of eggs associated with the next years cohort (fed into the stock recruit relationship).
#' @param E0,h,R0  Beverton-Holt parameters (E0= egg production at B0, h= steepness, R0 = recruitment at B0).
#' @param alpha,beta Direct Beverton-Holt parameters (alpha = per capita recruitment as egg production approaches 0; beta = inverse of recruitment asymptote).
#' @param stage_maturity stage class at maturity.
#' @param errors a vector of length L or single value of log-scale recruitment deviations.
#' @param n_stages number of stage classes.
#' @param n_loc number of locations. 
#' @param inst_h a matrix vector or single value with instantaneous fishing mortality
#' @param X0 an S x L matrix of initial abundances of each stage (S) at each location (L).
#' @param method the method of numerical integration (defaults to "lsoda").
#' @example /inst/examples/survexample_ode.R
#'@description Generates a single projection including simultaneous survival, straying and harvest using a system of ordinary differential equations followed by discrete reproduction 
SSR_linear_ode <- function(stray, Z, fec_at_age, eggs, 
                                E0, h, R0, 
                                alpha= NULL,beta= NULL,
                                stage_maturity, errors = 0, 
                                n_loc , n_stages,
                                inst_h = 0, X0, method = "lsoda") {
    S <- n_stages - 1
    L <- n_loc
    
    X1 <- array(as.numeric(NA), dim = c(n_stages, n_loc))
    
    col <- rep(1:(S * L), L)
    row <- rep(c(1:S), L * L) + rep(0:(L - 1), each = S * L) * S
    
    Zmat <- diag(as.vector(-Z))
    
    Smat <- matrix(simple_triplet_matrix(row, col, rep(t(stray), each = S)), ncol = 45)
    
    diag(Smat) <- rep(diag(stray), each = S) - 1
    
    Fmat <- diag(as.vector(inst_h), ncol = (L * S), nrow = (L * S))
    
    A <- rbind(cbind(Zmat + Smat, matrix(0, nrow = 45, ncol = 45)), 
               cbind(Fmat, matrix(0, nrow = 45, ncol = 45)))
    
    inits <- matrix(X0)
    
    ### generate the survivors and harvest totals ###
    fin <- array(t(ode(y = c(X0, rep(0, (S * L))), times = c(0, 1), 
                       func = linear_odes, parms = A, method = method)[2, -1]), 
                 dim = c(S, L, 2))
    
    ### split survivors and harvests ###
    harvest <- fin[, , 2]
    survivors <- fin[, , 1]
    survivors[S-1, ] <- survivors[S-1, ] + survivors[S, ]    
    
    X1[3:10,] <- survivors[1:8, ]
    X1[1, ] <- fec_at_age %*% X1[stage_maturity:n_stages, ]
    if (is.null(alpha)){
      X1[2, ] <- mapply(BH, E = eggs, E0 = E0, h = h, R0 = R0) * exp(errors)
    }
    else {
      X1[2, ] <- mapply(BH, E = eggs, alpha= alpha,beta=beta) * exp(errors)
    }
    return(list(X = X1, harvest = harvest))
} 
