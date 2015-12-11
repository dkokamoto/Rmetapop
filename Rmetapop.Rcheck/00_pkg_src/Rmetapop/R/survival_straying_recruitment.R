#' Discrete survival, straying, and recruitment projection.
#' @param stray_mat a square L x L stray matrix. 
#' @param surv_array an array of S-1 x S-1 survival matrices of length L (i.e. an array of off diagonal matrices) OR a single S-1 x S-1 survival matrix if constant in space OR a single survival value if constant in space and across ages.
#' @param fec_at_age a vector or matrix of fecundity at age including only mature stages.
#' @param eggs a vector of eggs associated with the next years cohort (fed into the stock recruit relationship).
#' @param E0,h,R0  Beverton-Holt parameters (E0= egg production at B0, h= steepness, R0 = recruitment at B0).
#' @param alpha,beta Direct Beverton-Holt parameters (alpha = per capita recruitment as egg production approaches 0; beta = inverse of recruitment asymptote).
#' @param stage_maturity stage class at maturity.
#' @param errors a vector of length L or single value of log-scale recruitment deviations.
#' @param X an S x L matrix of initial abundances of each stage (S) at each location (L).
#' @param n_stages number of stage classes.
#' @param n_loc number of locations. 
#' @param spat_alloc function of spatial allocation of harvest within a stock (currently only the extremes of 1="diffuse", 2="max", or 3="IDF")
#' @param group_dd groupings for grouped density dependence in the SR function.
#' @example /inst/examples/survexample.R
#' @description Generates a single projection including survival followed by straying and reproduction using a matrix algebra in discrete time.
#' @seealso surv_stray_recr_ode

ssr_linear <- function(stray_mat, surv_array, fec_at_age, 
                            eggs, E0, h , R0, 
                            tons_at_age,harvest,stage_selectivity,
                            stage_maturity,X0, n_stages, n_loc,
                            N_s=NULL,s_ID= NULL, 
                            alpha=NULL,beta=NULL,
                            group_dd=1,
                            spat_alloc= 1,
                            errors = 0,
                            stocklet_list= stocklet_list,
                            const_harvest= FALSE,
                            Fmort= NULL) {
    
    ### create empty array
    X1 <- array(as.numeric(NA), dim = c(n_stages, n_loc))
    
    if (length(dim(surv_array)) > 2) {
    
    ### fill in mortality at age
        for (j in 1:n_loc) {
            X1[(stage_maturity - 1):n_stages, j] <- surv_array[, , j] %*% X0[, j]
        }
    } else {
        X1[(stage_maturity - 1):n_stages, ] <- surv_array %*% X0
    }
    
    X1[(stage_maturity - 1):n_stages, ] <- X1[2:n_stages, ] %*% t(stray_mat)
    
    if(sum(harvest[1:n_stocks])>0&!is.na(sum(harvest[1:n_stocks]))){
      ### HARVEST with density independent selectivity among stages ###
      stage_selectivity= rep(1,length(stage_mat:n_stages))
      
      ### harvestable biomass
      h_biomass <-  apply(X1[stage_maturity:n_stages, ],2,function(x)x*tons_at_age)
      
      ### local harvestable biomass 
      loc_biomass <- colSums(h_biomass)
      stock_biomass <- sapply(stocklet_list,function(x) sum(loc_biomass[x]))
      if(const_harvest == FALSE){
        harvest <- mapply(min, x= harvest,y=stock_biomass)
      ### generate local harvest given quota, local biomass and fleet allocation 
        if(spat_alloc == 1){
          effort_alloc <- as.vector(sapply(stocklet_list,function(x) as.numeric(loc_biomass[x])/sum(as.numeric(loc_biomass[x]))))
        }
        if(spat_alloc == 2){
          effort_alloc <- as.vector(sapply(stocklet_list,
               function(x) as.numeric(x==x[which(loc_biomass[x]==max(loc_biomass[x]))])))
        }
        if(spat_alloc == 3){
          effort_alloc <- rep(0,n_loc)
          for (i in 1:n_stocks) {
            if (harvest[[i]]>0){
              Harvest <- set_effort(quota = harvest[[i]],biomass= loc_biomass[stocklet_list[[i]]],n_times= 2)$Harvest
              effort_alloc[stocklet_list[[i]]] <- Harvest/sum(Harvest)
            } else {
              effort_alloc[stocklet_list[[i]]] <- 0
            }
          }
        }
        
        loc_harvest <- as.vector(t(apply(matrix(effort_alloc,ncol= n_stocks),1,function(x)x*harvest)))
        post_harvested <- (h_biomass-t(apply(h_biomass,1,
                                             function(x) x/colSums(h_biomass)))*
                             (stage_selectivity%*%t(loc_harvest)))/tons_at_age  
      } else {
        loc_harvest <- colSums(Fmort*h_biomass)
        post_harvested <- (h_biomass-Fmort*h_biomass)/tons_at_age  
      }
      X1[stage_maturity:n_stages, ]<- mapply(max,post_harvested,0)
    } else {
      loc_harvest <- rep(0,n_stocks)
    }
    
    X1[1, ] <- mapply(max, 0.5*(fec_at_age %*% X1[stage_maturity:n_stages, ]),0)
    
    # add recruitment 
    if (is.null(alpha)){
      if(!is.list(group_dd)){
        X1[2, ] <- mapply(BH, E = eggs, E0 = E0, h = h, R0 = R0) * exp(errors)
      } else {
        X1[2, ] <-  mapply(BH,E= eggs,group_E= as.vector(sapply(group_dd,function(x) rep(sum(eggs[x]),each = length(x)))),E0 = E0, h = h, R0 = R0/length(group_dd[[1]]))  * exp(errors)
      } 
    } else {
      if(!is.list(group_dd)){
        X1[2, ] <- mapply(BH, E = eggs, alpha= alpha,beta=beta) * exp(errors)
      } else {
        X1[2, ] <-   mapply(BH,E= eggs,group_E= as.vector(sapply(group_dd,function(x) rep(sum(eggs[x]),each = length(x)))),alpha= alpha, beta= beta/length(group_dd[[1]]))* exp(errors)
      }
    }
  return(list(ages=X1, harvest= loc_harvest))
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
ssr_linear_ode <- function(stray, Z, fec_at_age, eggs, 
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
    
    Smat <- matrix(simple_triplet_matrix(row, col, rep(t(stray), each = S)), ncol = S*L)
    
    diag(Smat) <- rep(diag(stray), each = S) - 1
    
    Fmat <- diag(as.vector(inst_h), ncol = (L * S), nrow = (L * S))
    
    A <- rbind(cbind(Zmat + Smat, matrix(0, nrow = (L * S), ncol = (L * S))), 
               cbind(Fmat, matrix(0, nrow = (L * S), ncol = (L * S))))
    
    inits <- matrix(X0)
    
    ### generate the survivors and harvest totals ###
    fin <- array(t(ode(y = c(X0, rep(0, (S * L))), times = c(0, 1), 
                       func = linear_odes, parms = A, method = method)[2, -1]), 
                 dim = c(S, L, 2))
  
    ### split survivors and harvests ###
    harvest <- fin[, , 2]
    survivors <- fin[, , 1]
    survivors[S-1, ] <- survivors[S-1, ] + survivors[S, ]    

    X1[stage_maturity:n_stages,] <- survivors[1:8, ]
    X1[1, ] <- fec_at_age %*% X1[stage_maturity:n_stages, ]
    if (is.null(alpha)){
      X1[2, ] <- mapply(BH, E = eggs, E0 = E0, h = h, R0 = R0) * exp(errors)
    } else {
      X1[2, ] <- mapply(BH, E = eggs, alpha= alpha,beta=beta) * exp(errors)
    }
    return(list(X = X1, harvest = harvest))
} 
