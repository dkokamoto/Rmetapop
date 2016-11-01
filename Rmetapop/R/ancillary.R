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
BH <- function(E, E0=NULL, R0=NULL,h=NULL,  alpha = NULL, beta = NULL,group_E=NULL) {
    if(is.null(group_E)){
      group_E <- E
    } else {
      group_E <- group_E
    }
    if (!is.null(h)) {
        alpha_h <- (E0 * (1 - h))/(4 * h * R0)
        beta_h <- (5 * h - 1)/(4 * h * R0)
        E/(alpha_h + beta_h * group_E)
    } else {
        E/(alpha + beta * group_E)
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


#' A wrapper to simulate a linear stage structued model with a stock assessment and harvest
#' @param n_loc  number of locations at which a stock assessment is implemented
# #n_subloc <- rep(10,5) ### vector of number of sublocations
#' @param n_iter  number of years to simulate, including the warmup period
#' @param n_stages number of stages including eggs and recruits
#' @param stage_mat stage # that indicates maturity 
#' @param warmup  number of warmup iterations prior to starting a fishery  
#' @param M natural mortality (either a vector, matrix, or single value)
#' @param point.estimate logical for whether to estimate a point estimate or a fully Bayesian posterior
#' @param obs_sd observation standard deviation 
#' @param spat_scale spatial scale of stray connectivity
#' @param C degree of stochasticity in the stray matrix where higher numbers are less stochastic
# ### 
#' @param alpha BH stock-recruit parameter
#' @param beta BH stock-recruit parameters
#' @param phi first order recruitment autocorrelation 
#' @param spat_sd gaussian spatial correlation parameter
#' @param site_sd overall log-scale recruitment standard deviation 
#' @example /inst/examples/assessment_example.R

fishery_simulate <- function(n_loc,stock_IDs,
                             a_bh, b_bh, phi, 
                             spat_scale,spat_sd,
                             site_sd,C,obs_sd,
                             M=0.47,
                             n_iter=53,warmup=10,
                             warmup2= 13,
                             spat_alloc= 1,
                             surv_rho=0.05,
                             cor_mat=NULL,
                             collective_dd=FALSE,
                             h_crit=0.25,Fmort=0.2,
                             assessment=FALSE,
                             point.estimate= FALSE,
                             n_stages= 10,stage_mat= 3,
                             fit.compile= NULL,
                             repetitions = 1,
                             const_harvest=FALSE,
                             ret_ts= FALSE){

  n_stocks <- length(unique(stock_IDs))
  
  ### list of unique stocklets ###
  stocklet_list <- split(1:n_loc, stock_IDs)

  
  ### 9 x L mortality matrix (can also be a single number or a single row with L
  ### columns)
  
  mort <- matrix(M, ncol = n_loc, nrow = length((stage_mat - 1):(n_stages)))
  
  surv_array <- array(0, dim = c(n_stages - 1, n_stages - 1, n_loc))
  
  ### create identifier of survival array entries that correspond to each survival rate for quick
  ### replacement in the array
  surv_array_id <- array(matrix(1:n_stages^2, ncol = n_stages) %in% 
                           c(diag(matrix(1:n_stages^2,ncol = n_stages)[-c(1), ]), n_stages^2), 
                         dim = c(n_stages, n_stages, n_loc))[-1,-1, ]
  
  ### fecundity at age
  fec_at_age <- fecundity_age(stage_mat:(n_stages))
  
  ### biomass at age in tons (by the ten thousand)
  tons_at_age <- LVBweight(stage_mat:(n_stages))/1000/10000
  
  # matrix of mean stray probabilities, controled by a Cauchy distribution with scale spat_scale
  stray_probs <- spat_cor_mat(n_loc, spat_scale = spat_scale, sumto1 = TRUE)
  #round(stray_probs,2)# columns are probability of coming from, and row is probability going to (i.e. column 2 row 1 is going from location 2 to location 1)
  #colSums(stray_probs) # check to make sure it sums to 1
  
  ### stochastic realizations of the stray matrix for each year using the Dirichlet
  ### distribution with scale = 1000 (higher = less stochastic probabiliies)
  Crand <- ran_stray_prob(stray_mat=stray_probs,n_iter=n_iter,scale= C)
  #Crand[,,1:3]  #examine the first three stray matrices in the array
  mean_stray <- apply(Crand,c(1,2),mean)
  ### create matrices and vectors from parameters ### create spatiotemporal
  ### recruitment errors
  if(site_sd!=0){
    errors <- spat_temp_ts(n_iter = n_iter, n_loc = n_loc,
                         site_sd = site_sd, spat_sd =spat_sd, 
                         phi = phi,cor_mat=NULL)
    error_mat <- errors$ts  ### save the ts for use in the simulation 
  } else {
    error_mat <- matrix(0,ncol= n_loc,nrow= n_iter)
  }
  
  if(assessment== FALSE){
    if( obs_sd !=0){
          obs_error_mat <- matrix(spat_temp_ts(n_iter = n_iter, n_loc = 2,
                               site_sd = obs_sd, spat_sd =0.00001, 
                               phi = 0,cor_mat=NULL)$ts,ncol= n_stocks)
    } else{
      obs_error_mat <- matrix(0,ncol= n_stocks,nrow= n_iter)
    }
  }
  
  ### initial conditions ### create the data arrays
  d1 <- c(0, 2:n_stages)
  d2 <- paste("Location", 1:n_loc, sep = "_")
  d2b <- paste("Stock",1:n_stocks,sep= "_")
  d3 <- 1:n_iter
  
  ## set up the abundance array
  X <- array(as.numeric(NA), dim = c(n_stages, n_loc, n_iter), dimnames = list(d1, d2, d3))
  X[, , 1:2] <- matrix(rep(1e+08, each = n_loc), ncol = n_loc, byrow = T)
  
  ## set up the adult biomass array
  B <- array(NA, dim = c(n_loc, n_iter), dimnames = list(d2, 1:(n_iter)))
  B_s <- array(NA, dim = c(n_stocks, n_iter,2), dimnames = list(d2b, 1:(n_iter),c("actual","observed")))
  B[, 1] <- tons_at_age %*% X[3:10,, 1]
  B[, 2] <- tons_at_age %*% X[3:10,, 2]
  
  for(i in 1:n_stocks){
    B_s[i, 1,] <- sum(B[stock_IDs==i, 1])
    B_s[i, 2,] <- sum(B[stock_IDs==i, 2])
  }
  
  ## set up the Forecast adult biomass array
  BF <- array(NA, dim = c(n_stocks, n_iter+1, n_iter,6), dimnames = list(d2b, 1:(n_iter+1), d3,c("Estimated Stock Size","Prediction","ESS_UP","ESS_LOW","P_UP","P_LOW")))
  
  ## set up the age frequency arrays
  S_freq <- array(NA, dim = c(n_stages-2,n_loc, n_iter), dimnames = list(d1[-c(1,2)],d2,d3))
  S_freq_s <- array(NA, dim = c(n_stages-2,n_stocks, n_iter), dimnames = list(d1[-c(1,2)],d2b,d3))
  
  for (i in 1:2){
    S_freq[,,i]  <- t(t(X[-c(1,2),,i])/colSums(X[-c(1,2),,i]))
    SN <- mapply(function(x) rowSums(X[-c(1,2),stock_IDs==x,i]),1:n_stocks)
    S_freq_s[,,i]  <- t(t(SN)/colSums(SN))
  }
  
  H <- array(NA, dim = c(n_iter+1,n_stocks), dimnames = list(1:(n_iter+1),d2b))
  H_loc <- array(NA, dim = c(n_iter+1,n_loc), dimnames = list(1:(n_iter+1),d2))
  
  for (i in 1:n_iter){
    H[i,]  <- 0
    H_loc[i,] <- 0
  }
  
  if(collective_dd==0){
    group_dd <- stocklet_list
  }else{
    group_dd <- FALSE
  }
  
  for (i in 3:n_iter) {
    surv_array[surv_array_id] <- exp(-mort)
    X[, , i]  <-  ssr_linear(alpha=a_bh,beta=b_bh, fec_at_age = fec_at_age, 
                             n_loc = n_loc, 
                             n_stages = n_stages,
                             harvest=rep(0,n_stocks),
                             tons_at_age=tons_at_age,
                             stage_selectivity= rep(1,length(stage_mat:n_stages)),
                             stage_maturity = stage_mat,
                             surv_array = surv_array,              
                             group_dd= group_dd,
                             stocklet_list=stocklet_list,
                             spat_alloc= spat_alloc,
                             s_ID=stock_IDs,N_s=n_stocks,
                             eggs = X[1, , i - 2], X0 = X[(stage_mat - 1):n_stages, , i - 1],
                             stray_mat = stray_probs, errors =  0)$ages
    
    B[, i] <- tons_at_age %*% X[3:10, , i]
    B_s[, i,] <- mapply(function(x) sum(B[stock_IDs==x,i]),1:n_stocks) 
    ### age frequency ###
    SN <- mapply(function(x) rowSums(X[-c(1,2),stock_IDs==x,i]),1:n_stocks)
    S_freq[,,i]  <- t(t(X[-c(1,2),,i])/colSums(X[-c(1,2),,i]))
    S_freq_s[,,i]  <- t(t(SN)/colSums(SN))
  } 
  
  ### estimate B0  for each stock and stocklet
  B0 <- sapply(stocklet_list, function(x) sum(B[x,n_iter],na.rm= T))
  B0_loc <- B[,n_iter]
  
  ### estimate beta for the BH parameters
  model_beta <- 1/mean(X[2,,n_iter-1])/(n_loc/n_stocks)
  
  ### reset initial values as last two values of deterministic simulation 
  
  B[, 1] <- B[, n_iter-1]
  B[, 2] <- B[, n_iter]
  B_s[, 1,] <- B_s[, n_iter-1,]
  B_s[, 2,] <- B_s[, n_iter,]
  X[, , 1] <- X[,,n_iter-1]
  X[, , 2] <- X[,,n_iter]
  for (i in 1:2){
    S_freq[,,i]  <- S_freq[,,n_iter-2+i]
    S_freq_s[,,i]  <- S_freq_s[,,n_iter-2+i]
  }
  ### project the population with stochastic recruitment and harvesting ###repea?
  for (i in 3:n_iter) {  
    ptm <-proc.time()
    #surv_array[surv_array_id] <- ran_surv_prob(mean=exp(-mort),corr=surv_rho)
    surv_array[surv_array_id] <- exp(-mort)
    projection <- ssr_linear(alpha=a_bh,beta=b_bh, 
                             fec_at_age = fec_at_age, 
                             n_loc = n_loc,
                             n_stages = n_stages,
                             N_s=n_stocks,
                             s_ID=stock_IDs,
                             spat_alloc=spat_alloc,
                             group_dd= group_dd,
                             stocklet_list=stocklet_list,
                             tons_at_age=tons_at_age,
                             harvest=H[i,],
                             stage_maturity = stage_mat,
                             surv_array = surv_array,
                             stage_selectivity= rep(1,length(stage_mat:n_stages)),
                             eggs = X[1, , i - 2], 
                             X0 = X[(stage_mat - 1):n_stages, , i - 1],
                             stray_mat = Crand[,, i - 1], errors = error_mat[i - 1, ],
                             const_harvest =const_harvest,Fmort= Fmort)

    X[, , i] <- projection$ages            
    H_loc[i,] <- projection$harvest
    
    B[, i] <- tons_at_age %*% X[3:10, , i]
    B_s[, i,1] <- sapply(stocklet_list,function(x) sum(B[x,i]))
    B_s[, i,2] <- B_s[, i,1]* exp(rnorm(n_stocks, 0, obs_sd) - 0.5 * obs_sd^2)
    
    ### age frequency ###
    SN <- mapply(function(x) rowSums(X[-c(1,2),stock_IDs==x,i]),1:n_stocks)
    S_freq[,,i]  <- t(t(X[-c(1,2),,i])/colSums(X[-c(1,2),,i]))
    S_freq_s[,,i]  <- t(t(SN)/colSums(SN))
    
    ### harvest allocation
    if(i>warmup){  
      if (assessment== TRUE){
        if(i==warmup+1){
          start <- max(1,i-20)
          assess <- nlss_assess(b_obs=B_s[, start:i,2],b_true=B_s[, start:i,1],sd_R=site_sd,
                                age_freq=S_freq_s[,,start:i],obs_sd=obs_sd,
                                alpha_BH= a_bh,beta_BH=model_beta,harvest= H[start:i,],
                                selectivity=rep(1,length(stage_mat:n_stages)),n_stocks=n_stocks,
                                fec_at_age=fec_at_age,weight_at_age=tons_at_age,
                                mort=rep(M,n_stocks),stage_mat=stage_mat,
                                n_stages=n_stages,phi=phi,n_warmup=3,compile= FALSE,
                                optimize=point.estimate,compiled_model= fit.compile)
        } else {
          assess <- nlss_assess(b_obs=B_s[, start:i,2],b_true=B_s[, start:i,1],sd_R=site_sd,
                                age_freq=S_freq_s[,,start:i],obs_sd=obs_sd,
                                alpha_BH= a_bh,beta_BH=model_beta,harvest= H[start:i,],
                                selectivity=rep(1,length(stage_mat:n_stages)),n_stocks=n_stocks,
                                fec_at_age=fec_at_age,weight_at_age=tons_at_age,
                                mort=rep(M,n_stocks),stage_mat=stage_mat,
                                n_stages=n_stages,phi=phi,n_warmup=3,compile= FALSE,
                                optimize=point.estimate,compiled_model= fit.compile,
                                sd_pro = sd_process)
        }
        if(point.estimate==TRUE){
          BF[, start:(i), i,1] <- t(assess$assess$mean)
          BF[, start:(i+1), i,2] <- t(assess$assess$pred)
        }
        else{
          BF[, start:(i), i,1] <- t(assess$assess$mean)
          BF[, start:(i+1), i,2] <- t(assess$assess$pred)
          BF[, start:(i), i,3] <- t(assess$assess$upper)
          BF[, start:(i), i,4] <- t(assess$assess$lower)
          BF[, start:(i+1), i,5] <- t(assess$assess$pred_upper)
          BF[, start:(i+1), i,6] <- t(assess$assess$pred_lower)
        }
        
        sd_process <- assess$sd_pro
        cat(noquote(paste("iteration",i, "of", n_iter, "completed", signif(proc.time()[3]-ptm[3],digits=4),"sec\n")))
        ### forecast array 
      } else {          
        stock_surv <- lapply(stocklet_list,function(x) apply(surv_array[,,x],c(1,2),mean))
        for (j in 1:n_stocks) {   
          BF[j, i+1, i,2] <- tons_at_age%*%(stock_surv[[j]]%*%sapply(stocklet_list,function(x) rowSums(X[(stage_mat - 1):n_stages,x , i]))[,j])[2:9]*exp(obs_error_mat[i,j])
        }
      }
      
      to_fish <- ifelse(h_crit*B0<BF[,i+1,i,2],1,0)
      H[i+1,] <- Fmort*BF[,i+1,i,2]*to_fish
    } else {
      if(i==warmup){
        if(assessment== TRUE){
          cat(noquote(paste("iterations 1:",i," of ",n_iter," complete (warmup)\n",sep= "")))
        }
      }
      H[i+1,] <- rep(0,n_stocks) 
    }
  }  
  
  ### generate summaries from the simulation 
  
  stock_catch =matrix(sapply(stocklet_list,function(x) rowSums(H_loc[,x])),ncol= n_stocks)
  g_mu_catch =ifelse(n_stocks>1,mean(rowSums(matrix(stock_catch[warmup2:n_iter,]),na.rm= T)),NA)
  s_mu_catch = mean(stock_catch[warmup2:n_iter,])
  
  metrics_df <- data.frame(list(
    p_cv =mean(apply(B[,warmup2:n_iter],1,sd)/apply(B[,warmup2:n_iter],1,mean)),
    s_cv =mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),1,sd)/apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks) ,1,mean)),
    g_cv =ifelse(n_stocks>1,
                 sd(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,mean))/mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,mean)),
                 NA),
    
    ### mean and equilibrium biomass
    p_mu =mean(apply(B[,warmup2:n_iter],1,mean)),
    s_mu =mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),1,mean))),
    g_mu =ifelse(n_stocks>1,mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,mean)),NA),
    
    s_b0 =mean(B0_loc),
    p_b0 =mean(B0),
    g_b0 =ifelse(n_stocks>1,sum(B0),NA),
    
    ### overfishing metrics
    p_of =mean(apply(B[,warmup2:n_iter],2,function(x) x<=h_crit*B0_loc)),
    s_of =mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,function(x) x<=h_crit*B0)),
    g_of =ifelse(n_stocks>1,
                 mean(!(colMeans(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,function(x) x<=h_crit*B0))<1)),
                 NA),
    
    ### mean depletion length
    g_dep_l =dep_length(colMeans(apply(matrix(B_s[,warmup2:n_iter,1]),2,function(x) x<=h_crit*B0)!=0)),
    s_dep_l =mean(apply(apply(matrix(B_s[,warmup2:n_iter,1]),2,function(x) x<=h_crit*B0),1,dep_length),na.rm= T),
    p_dep_l =ifelse(n_stocks>1,
                    mean(apply(apply(B[,warmup2:n_iter],2,function(x) x<=h_crit*B0_loc),1,dep_length),na.rm= T),
                    NA),
    
    ### mean annual catch 
    g_mu_catch = ifelse(n_stocks>1,g_mu_catch,NA),
    s_mu_catch = s_mu_catch,

    ### average annual variation in catch 
    g_aav =ifelse(n_stocks>1,mean(abs(diff(rowSums(H_loc[warmup2:n_iter,]))),na.rm= T)/g_mu_catch,NA),
    s_aav =mean(apply(matrix(stock_catch[warmup2:n_iter,]),2,function(x) sd(diff(x),na.rm= T)))/(s_mu_catch),
    p_aav =mean(apply(H_loc[warmup2:n_iter,],2,function(x) sd(diff(x),na.rm= T)))/(s_mu_catch/n_loc),
    
    ### average variability in catch 
    g_catch_sd =ifelse(n_stocks>1,sd(rowSums(H_loc[warmup2:n_iter,]),na.rm=T)/g_mu_catch,NA),
    s_catch_sd =mean(apply(matrix(stock_catch[warmup2:n_iter,]),2,function(x) sd(x,na.rm= T)))/(s_mu_catch),
    p_catch_sd =mean(apply(H_loc[warmup2:n_iter,],2,sd,na.rm= T))/(s_mu_catch/n_loc),
    
    ### autocorrelaton in biomass 
    g_phi =ifelse(n_stocks>1,pacf(colMeans(matrix(B_s[,,1],nrow= n_stocks),na.rm=T),na.action=na.pass,plot= FALSE)$acf[1],NA),
    s_phi =mean(apply(matrix(B_s[,,1],nrow= n_stocks),1,function(x) pacf(x,na.action=na.pass,plot= FALSE)$acf[1])),
    p_phi =mean(apply(B,1,function(x) pacf(x,na.action=na.pass,plot= FALSE)$acf[1])),
    
    ### spatial variance 
    spat_cv =mean(apply(B[,warmup2:n_iter],2,function(x) sd(x)/mean(x)))
  )
  
  ### spectral density at different frequencies in biomass
  if(mean(B[,warmup2:n_iter])==0){
    b <- NA; s<- NA; g <- NA;
  } else {
    p <- spec.pgram(log(t((B[,warmup2:n_iter]))),plot= FALSE)
    freq <- p$freq
    p <- rowMeans(p$spec)
    s <- rowMeans(apply(log(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks)),1,function(x) spec.pgram(x,plot= FALSE)$spec))
    g <-ifelse(n_stocks>1, spec.pgram(log(colMeans(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks))),plot= FALSE)$spec,NA)
  }
  spec_df <- data.frame(pop_spec=p,stock_spec=s,glob_spec= g,freq= freq)
  
  probs =seq(0.05,0.95,0.05)
  ### biomass quantiles
  quantile_df <- data.frame(list(
    p_b =rowMeans(apply(B[,warmup2:n_iter],1,quantile, probs)),
    s_b = rowMeans(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),1,quantile, probs)),
    g_b = ifelse(n_stocks>1,quantile(colSums(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks)), probs),NA),
    p_c = rowMeans(apply(H_loc,2,quantile, probs,na.rm= T)),
    s_catch = rowMeans(apply(stock_catch,2,quantile, probs,na.rm=T)),
    g_c = ifelse(n_stocks>1,quantile(rowSums(stock_catch,na.rm=T), probs,na.rm=T),NA)
    )
  )
  
  r_sigma <- sd(log(apply(X[2,,warmup2:n_iter],2,sum)))
  
  quantile_df$quants <- row.names(quantile_df)
  if(ret_ts == TRUE){
    TS <- list(B=B,assess= BF,B_stocks= B_s,ages= X,B0=B0, Harvest= H, forecast= BF)
  } else {
    TS <- NULL
  }
 
  return(list(metrics= metrics_df,spec= spec_df,quantiles= quantile_df,r_sigma= r_sigma, TS=TS))
  
}




