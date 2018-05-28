
rows.to.list <- function( df ) {
  ll<-apply(df,1,list)
  ll<- lapply(ll,function(x) as.data.frame(t(unlist(x))))
}

dep_length <- function(x) {
  a <- rle(as.numeric(x))
  ifelse(length(a$lengths[a$values==1])==0,NA,mean(a$lengths[a$values==1]))
}

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

### function to read in ADMB output
#' @param fn ADMB output file 
reptoRlist=function(fn)
{
  ifile=scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  idx=sapply(suppressWarnings(as.double(ifile)),is.na)
  vnam=ifile[idx] #list names
  nv=length(vnam) #number of objects
  A=list()
  ir=0
  for(i in 1:nv)
  { 
    ir=match(vnam[i],ifile)
    if(i!=nv)
      irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
      dum=NA
      if(irr-ir==2)
        dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
      if(irr-ir>2)
        dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
      if(is.numeric(dum))#Logical test to ensure dealingwith  numbers
      {
        A[[ vnam[i]]]=dum
      }
  }
  return(A)
}


### Allometric functions
#' Ludwig Von-Bertalanfy growth equation.
LVB <- function(age, L_inf = 27, k = 0.48, t0 = 0) {
  L_inf * (1 - exp(-k * (age - t0)))
}

#' Power length-weight relationship (from Martell et al. 2012)
weight <- function(length, a = 4.5e-06, b = 3.127) {
  a * (length)^b
}

#' Weight at age using the LVB growth equation and length-weight relationship.
LVBweight <- function(age) {
  weight(LVB(age))
}

#' Fecundity at weight function with base parameters from Tanasichuk et al. (1993)  
fecundity <- function(w, F1 = exp(4.69),F2 = 1.13) {
  F1 * (w*1000)^(F2) 
}

#' Fecundity at age using LVB growth equation and fecundity at length function
fecundity_age <- function(age) {
  fecundity(LVBweight(age))
}

#' Define a system of ordinary differential equations from a matrix
linear_odes <- function(t, state, A) {
  dX <- A %*% state
  list(as.vector(dX))
} 


#' A wrapper to simulate a linear stage structued model with a stock assessment and harvest
#' @param n_loc  number of locations (populations) simulated
#' @param stock_IDs vector of integers identifying which populations belong to which stock
#' @param harvest_IDs vector of integers identifying which populations are harvested within a stock (NA iS not harvested)
# #n_subloc <- rep(10,5) ### vector of number of sublocations
#' @param n_iter  number of years to simulate, including the warmup period
#' @param n_stages number of stages including eggs and recruits
#' @param stage_mat stage # that indicates maturity 
#' @param warmup  number of warmup iterations prior to starting a fishery  
#' @param M natural mortality (either a vector, matrix, or single value)
#' @param point.estimate logical for whether to estimate a point estimate or a fully Bayesian posterior
#' @param obs_sd observation standard deviation 
#' @param spat_scale spatial scale of stray connectivity
#' @param juv_spat_scale spatial scale of juvenile stray connectivity
#' @param C degree of stochasticity in the stray matrix where higher numbers are less stochastic
# ### 
#' @param alpha BH stock-recruit parameter
#' @param beta BH stock-recruit parameters
#' @param phi first order recruitment autocorrelation 
#' @param spat_sd gaussian spatial correlation parameter
#' @param site_sd overall log-scale recruitment standard deviation 
#' @param follow do young fish follow older fish NULL or TRUE
#' @example /inst/examples/assessment_example.R

fishery_simulate <- function(n_loc,
                             stock_IDs,
                             harvest_IDs,
                             a_bh, b_bh,
                             phi, 
                             spat_scale,
                             juv_spat_scale,
                             spat_sd,
                             site_sd,
                             C,
                             obs_sd,
                             M=0.47,
                             n_iter=53,warmup=10,
                             warmup2= 13,
                             spat_alloc= 1,
                             surv_rho=0.05,
                             cor_mat=NULL,
                             collective_dd=FALSE,
                             h_crit=0.25,
                             of_crit= 0.2,
                             Fmort=0.2,
                             assessment=FALSE,
                             point.estimate= FALSE,
                             stage_forecast= 2,
                             n_stages= 10,
                             stage_mat= 2,
                             maturity = c(0.1,0.9,rep(1,7)),
                             fit.compile= NULL,
                             repetitions = 1,
                             const_harvest=FALSE,
                             ret_ts= FALSE,
                             follow= NULL,
                             dd_IDs=NULL,
                             quota_adj= FALSE,
                             ns_samps=200,
                             nf_samps=200){
  
  n_stocks <- length(unique(stock_IDs[!is.na(stock_IDs)]))
  
  ### list of unique stocklets by stock###
  stocklet_list <- split(1:n_loc, stock_IDs,drop= TRUE)
  
  ### list of harvestable stocklets by stock###
  harvest_list <- split(1:n_loc, harvest_IDs)
  
  ### list of stocklets that compete ###
  if(is.null(dd_IDs)){
    dd_list <- split(rep(1:n_loc,rep(1,n_loc)),rep(1,n_loc))
  } else {
    dd_list <- split(1:n_loc, dd_IDs)
  }
  
  ### 9 x L mortality matrix (can also be a single number or a single row with L
  ### columns)
  
  mort <- matrix(M, ncol = n_loc, nrow = n_stages-1)
  
  surv_array <- array(0, dim = c(n_stages - 1, n_stages - 1, n_loc))
  
  ### create identifier of survival array entries that correspond to each survival rate for quick
  ### replacement in the array
  surv_array_id <- array(matrix(1:n_stages^2, ncol = n_stages) %in% 
                           c(diag(matrix(1:n_stages^2,ncol = n_stages)[-c(1), ]), n_stages^2), 
                         dim = c(n_stages, n_stages, n_loc))[-1,-1, ]
  
  ### fecundity at age
  fec_at_age <- fecundity_age(2:(n_stages))
  
  ### biomass at age in tons (by the ten thousand)
  tons_at_age <- LVBweight(2:(n_stages))/1000/1000
  
  # matrix of mean stray probabilities, controled by a Cauchy distribution with scale spat_scale
  stray_probs <- spat_cor_mat(n_loc, spat_scale = spat_scale, sumto1 = TRUE,offset= 0)
  juv_stray_probs <- spat_cor_mat(n_loc, spat_scale = juv_spat_scale, sumto1 = TRUE,offset= 0)
  
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
      obs_error_mat <- matrix(spat_temp_ts(n_iter = n_iter, n_loc = n_stocks,
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
  d4 <- c("pre-spawn","catch")
  ## set up the abundance array
  X <- array(as.numeric(NA), dim = c(n_stages, n_loc, n_iter), dimnames = list(d1, d2, d3))
  X[, , 1:2] <- matrix(rep(1e+08, each = n_loc), ncol = n_loc, byrow = T)
  
  mat <- array(as.numeric(NA), dim = c(n_stages-1, n_loc, n_iter), dimnames = list(d1[-1], d2, d3))
  
  ## set up the adult biomass index array
  B <- array(NA, dim = c(n_loc, n_iter,2), dimnames = list(d2, 1:(n_iter),c("actual","observed")))
  B_s <- array(NA, dim = c(n_stocks, n_iter,2), dimnames = list(d2b, 1:(n_iter),c("actual","observed")))
  B[, 1,1] <- (tons_at_age*maturity) %*% X[2:n_stages,, 1]
  B[, 2,1] <- (tons_at_age*maturity) %*% X[2:n_stages,, 2]
  
  for(i in 1:n_stocks){
    B_s[i, 1,] <- sum(B[stock_IDs==i, 1,1])
    B_s[i, 2,] <- sum(B[stock_IDs==i, 2,1])
  }
  
  ## set up the Forecast adult biomass array
  BF <- array(NA, dim = c(n_stocks, n_iter+1, n_iter,6), dimnames = list(d2b, 1:(n_iter+1), d3,c("Estimated Stock Size","Prediction","ESS_UP","ESS_LOW","P_UP","P_LOW")))
  
  ## set up the age frequency arrays
  S_freq <- array(NA, dim = c(n_stages-1,n_loc, n_iter), dimnames = list(d1[-1],d2,d3))
  S_freq_f <- array(NA, dim = c(n_stages-1,n_stocks, n_iter,2), dimnames = list(d1[-1],d2b,d3,c("actual","observed")))
  S_freq_s <- array(NA, dim = c(n_stages-1,n_stocks, n_iter,2), dimnames = list(d1[-1],d2b,d3,c("actual","observed")))
  
  ## set up the social outcome array
  SBI <- factor(rep(NA,n_iter), levels= 1:10)
  
  ### add first two years mature age frequencies
  for (i in 1:2){
    S_freq[,,i]  <- t(t(maturity*X[-1,,i])/colSums(maturity*X[-1,,i]))
    SN <- mapply(function(x) rowSums(matrix(maturity*X[-1,stock_IDs==x,i],nrow= n_stages-1),na.rm=TRUE),1:n_stocks)
    S_freq_f[,,i,1]  <- t(t(SN)/colSums(SN))
    S_freq_s[,,i,1]  <- t(t(SN)/colSums(SN))
    for(j in 1:n_stocks){
      S_freq_f[,j,i,2]  <-  rmultinom(1,ns_samps,S_freq_s[,j,i,1])
      S_freq_s[,j,i,2]  <-  rmultinom(1,ns_samps,S_freq_s[,j,i,1])
    }
  }
  
  H <- array(NA, dim = c(n_iter+1,n_stocks), dimnames = list(1:(n_iter+1),d2b))
  H_loc <- array(NA, dim = c(n_iter+1,n_loc), dimnames = list(1:(n_iter+1),d2))
  EH_loc <- array(NA, dim = c(n_iter+1,n_loc), dimnames = list(1:(n_iter+1),d2))
  
  red <- rep(NA,n_iter+1)
  
  for (i in 1:n_iter){
    H[i,]  <- 0
    H_loc[i,] <- 0
  }
  
  if(collective_dd==0){
    group_dd <- dd_list
  }else{
    group_dd <- FALSE
  }
  
  for (i in 3:n_iter) {
    surv_array[surv_array_id] <- exp(-mort)
    projection <- ssr_linear(alpha=a_bh,beta=b_bh, fec_at_age = fec_at_age, 
                             n_loc = n_loc, 
                             n_stages = n_stages,
                             harvest=rep(0,n_stocks),
                             tons_at_age=tons_at_age,
                             maturity= maturity,
                             stage_maturity = stage_mat,
                             surv_array = surv_array,              
                             group_dd= group_dd,
                             harvest_list=harvest_list,
                             spat_alloc= spat_alloc,
                             s_ID=stock_IDs,N_s=n_stocks,
                             eggs = X[1, , i - 2], 
                             X0 = X[2:n_stages, , i - 1],
                             stray_mat = stray_probs, errors =  0,
                             GWOF=follow,
                             juv_rand= 100000,
                             juv_stray_mat = juv_stray_probs)
    
    X[, , i]  <-  projection$ages
    mat[, , i] <- projection$mat
    
    SN <- mapply(function(x) rowSums(projection$ages_harvested),1:n_stocks)
    S_freq_f[,,i,1]  <-  t(t(SN)/colSums(SN))
    
    B[, i,1] <-  tons_at_age%*% projection$mat
    B_s[, i,] <- mapply(function(x) sum(B[stock_IDs==x,i,1]),1:n_stocks) 
    
    ### age frequency ###
    SN <- maturity*mapply(function(x) rowSums(X[-1,stock_IDs==x,i]),1:n_stocks)
    S_freq[,,i]  <-  t(t(maturity*X[-1,,i])/colSums(maturity*X[-1,,i]))
    S_freq_s[,,i,1]  <- t(t(SN)/colSums(SN))
  } 
  
  ### estimate B0  for each stock and stocklet
  B0 <- sapply(stocklet_list, function(x) sum(B[x,n_iter,1],na.rm= T))
  B0_est <- B0
  B0_loc <- B[,n_iter,1]
  
  ### estimate Age  for each stock and stocklet
  ages0 <- sapply(stocklet_list,function(x) rowSums(X[-1,x,n_iter]*maturity))
  ages0_loc <- X[-1,,n_iter]*maturity
  
  ### reset initial values as last two values of deterministic simulation 
  B[, 1,] <- B[, n_iter-1,1]
  B[, 2,] <- B[, n_iter,1]
  B_s[, 1,] <- B_s[, n_iter-1,]
  B_s[, 2,] <- B_s[, n_iter,]
  X[, , 1] <- X[,,n_iter-1]
  X[, , 2] <- X[,,n_iter]
  mat[, , 1] <- mat[,,n_iter-1]
  mat[, , 2] <- mat[,,n_iter]
  
  for (i in 1:2){
    S_freq[,,i]  <- S_freq[,,n_iter-2+i]
    S_freq_s[,,i,1]  <- S_freq_s[,,n_iter-2+i,1]
    
  }
  
  ### create sampling dataframes ###
  catch_comp <- cbind(expand.grid(stock = 1:n_stocks, fleet = c(1,2),year= 1:(n_iter-2),size= nf_samps,a=0),matrix(NA,ncol= n_stages-1))%>%
    arrange(stock,fleet,year)
  age_comp <- cbind(expand.grid(stock = 1:n_stocks, year= 1:(n_iter-2),size= ns_samps,a=0),matrix(NA,ncol= n_stages-1))
  index <- expand.grid(stock = 1:n_stocks, year= 1:(n_iter-2),index =NA,cv= sqrt(exp(obs_sd^2)-1))
  index.true <- expand.grid(stock = 1:n_stocks, year= 1:(n_iter-2),index.true=NA)
  catches <- expand.grid(stock = 1:n_stocks, fleet= c(1,2),year= 1:(n_iter-2),catch=NA)%>%
    arrange(stock,fleet,year)
  
  ### project the population with stochastic recruitment and harvesting ###repea?
  for (i in 3:n_iter) {  
    ptm <-proc.time()
    surv_array[surv_array_id] <- ran_surv_prob(mean=exp(-mort),corr=surv_rho)
    projection <- ssr_linear(alpha=a_bh,beta=b_bh, 
                             fec_at_age = fec_at_age, 
                             n_loc = n_loc,
                             n_stages = n_stages,
                             N_s=n_stocks,
                             s_ID=stock_IDs,
                             spat_alloc=spat_alloc,
                             group_dd= group_dd,
                             harvest_list=harvest_list,
                             tons_at_age=tons_at_age,
                             harvest=H[i,],
                             stage_maturity = stage_mat,
                             surv_array = surv_array,
                             maturity= maturity,
                             eggs = X[1, , i - 2], 
                             X0 = X[2:n_stages, , i - 1],
                             stray_mat = Crand[,, i - 1], 
                             errors = error_mat[i - 1, ],
                             const_harvest =const_harvest,
                             Fmort= Fmort,
                             GWOF=follow,
                             juv_rand= 100000,
                             juv_stray_mat = juv_stray_probs)
    
    ### record number of fish at each age
    X[, , i] <- projection$ages
    mat[, , i] <- projection$mat
    
    ### record frequency of fish harvested by the fleet at each age
    SN <- mapply(function(x) rowSums(projection$ages_harvested),1:n_stocks)
    S_freq_f[,,i,1]  <-  t(t(SN)/colSums(SN))
    
    ### record total biomass of fish harvested by the fleet
    H_loc[i,] <- projection$harvest
    
    ### record total biomass of eggs harvested in space
    EH_loc[i,] <- projection$eggs_harvest
    
    ### record total spawning biomass 
    B[, i,1] <- (tons_at_age %*% projection$mat)
    B[, i,2] <- B[, i,1]*exp(rnorm(n_stocks, 0, obs_sd) - 0.5 * obs_sd^2)
    B_s[, i,1] <- sapply(stocklet_list,function(x) sum(B[x,i,1]))
    B_s[, i,2] <- sapply(stocklet_list,function(x) sum(B[x,i,2]))
    
    ### spawning age frequency ###
    SN <- mapply(function(x) rowSums(projection$mat[,stock_IDs==x]),1:n_stocks)
    S_freq[,,i]  <-  t(t(projection$mat)/colSums(projection$mat))
    S_freq_s[,,i,1]  <- t(t(SN)/colSums(SN))
    
    if(i> 2){
      ### random multinomial sample 
      for(j in 1:n_stocks){
        age_comp[age_comp$stock==j&age_comp$year==i-2,-c(1:4)] <-  rmultinom(1,ns_samps,S_freq_s[,j,i,1]) 
        for (k in 1:2){
          if(i>(warmup)&!is.na(sum(S_freq_f[,j,i,1]))){
            catch_comp[catch_comp$stock==j&catch_comp$fleet==k&catch_comp$year==i-2,-c(1:5)] <- rmultinom(1,nf_samps,S_freq_f[,j,i,1])
          } else {
            catch_comp[catch_comp$stock==j&catch_comp$fleet==k&catch_comp$year==i-2,-c(1:5)] <- rmultinom(1,nf_samps,S_freq_s[,j,i,1])
          }
        }
      }
      
      index[index$year==i-2,"index"] <- B_s[, i,2]
      index.true[index.true$year==i-2,"index.true"] <- B_s[, i,1]
      for (k in 1:2){
        catches[catches$year==i-2&catches$fleet==k,"catch"] <- sapply(harvest_list,function(x) sum(H_loc[i,x]))/2
      }
    }
    
    ### harvest allocation
    if(i>start){  
      if(i>warmup&assessment== TRUE){
        for (j in 1:n_stocks){
          n_iter2= i-2
          if(i==(start+1)) {
            test_pin <- list(Dummy=0,
                             LogRbar= 5,
                             SigmaJ= 0.5,
                             SigmaA= 0.6,
                             Alpha= 5,
                             SigmaR= 0.8,
                             Rec_devs = matrix(0,ncol= n_iter2,nrow= n_stocks))
          } else {
            test_pin <- list(Dummy=0,
                             LogRbar= log(temp_out$"#mfexpLogRbar"),
                             SigmaJ= 0.5,
                             SigmaA= 0.6,
                             Alpha= temp_out$"#Alpha",
                             SigmaR=0.7,
                             Rec_devs = matrix(c(temp_out$"#Rec_dev",0),ncol= n_iter2,nrow= n_stocks))
          }
          data <- list(assessed_pops=n_stocks,
                       number_of_ages=9,
                       #initial_year = -40,
                       number_of_years= n_iter2,
                       Nproj= 1,
                       number_of_fleets= 2,
                       SR_type= 3,
                       phase_adult_diff = -1,
                       phase_juv_diff= -1,
                       PhaseSigmaR= -1,
                       last_rec = n_iter2,
                       catches_pop_fleet_year_catch=as.matrix(subset(catches,year<=n_iter2)),
                       index_pop_year_index_cv=  as.matrix(subset(index,year<=n_iter2)),
                       spawning_age_com_pop_year_SS= as.matrix(subset(age_comp,year<=n_iter2)),
                       catches_age_com_pop_year_SS=as.matrix(subset(catch_comp,year<=n_iter2)),
                       M= M,
                       maturity_at_age = c(0,maturity),
                       #proportion_maturing= ifelse(c(0,maturity)>=1,1,diff(c(0,maturity))/c(1-c(0,maturity[-9]))),
                       weight_at_age = c(LVBweight(1:(n_stages)))*1000,
                       fec_at_age= fecundity_age(1:n_stages)/100,
                       selectivity_at_age_fleet1 = rbind(c(0,maturity),c(0,maturity)),
                       Initial_numbers= 10,
                       #The_Year= n_iter2,
                       #PhaseAlpha= 1,
                       TestCodeDat= 123456)
          
          write_dat("Herring",data)
          write_pin("Herring",test_pin)
          system.time(run_admb("Herring",extra.args= "-est",
                               mcmc= FALSE,verbose= FALSE))
          temp_out <- reptoRlist("OutputFile1.Out")
          temp_par <- reptoRlist("Herring.par")
          
          #  } else {
          
          
          # ### for testing 
          # temp_out2 <- data.frame(temp_out$"#For.All.Output.Ages")
          # temp_out3 <- melt(temp_out2[,c(1,2)],id.vars= 1,variable.name= "pop",value.name= "est")
          # names(temp_out2) <- c("year","1","1_N")
          # temp_out3 <- melt(temp_out2[,c(1,2)],id.vars= 1,variable.name= "pop",value.name= "Andre_est")
          # names(temp_out3) <- c("year","stock","est")
          # test4 <- join_all(list(catches,index,index_real,temp_out3))
          #  ggplot(aes(year,index_real),data= test4)+
          #    geom_line(aes(colour= "true biomass"),size= 2)+
          #    geom_point(aes(y= index))+
          #    geom_line(aes(x= year,y=est))
          # ####
          
          BF[1, 3:i, i,1] <- subset(data.frame(temp_out$"#For.All.Output.Ages"),V1>0)$V2
          BF[1, i+1, i,2] <- temp_out$"#Quota.Calculation"[6]+temp_out$"#Quota.Calculation"[5]
          B0_est[i] <- temp_out$"#Virgin.conditions"[3]
          #}
          cat(noquote(paste("iteration",i, "of", n_iter, "completed", signif(proc.time()[3]-ptm[3],digits=4),"sec\n")))
        }   ### forecast array 
      } else {          
        for (j in 1:n_stocks) {   
          surv_array[surv_array_id] <- exp(-mort)
          BF[j, i, i,1] <- sum(B[stocklet_list[[j]],i,1])
          BF[j, i+1, i,2]<-(((tons_at_age*maturity)[(stage_forecast:n_stages)-1])%*%
                              sapply(stocklet_list,function(x) 
                                rowSums(ssr_linear(alpha=a_bh,beta=b_bh, fec_at_age = fec_at_age, 
                                                   n_loc = n_loc, 
                                                   n_stages = n_stages,
                                                   harvest=rep(0,n_stocks),
                                                   tons_at_age=tons_at_age,
                                                   maturity= maturity,
                                                   stage_maturity = stage_mat,
                                                   surv_array = surv_array,              
                                                   group_dd= group_dd,
                                                   harvest_list=harvest_list,
                                                   spat_alloc= spat_alloc,
                                                   s_ID=stock_IDs,N_s=n_stocks,
                                                   eggs = X[1, , i - 1], 
                                                   X0 = X[2:n_stages, , i],
                                                   stray_mat = stray_probs, 
                                                   errors =  0,
                                                   GWOF=follow,
                                                   juv_rand= 100000,
                                                   juv_stray_mat = juv_stray_probs)$ages[stage_forecast:n_stages,x])))[j]
          B0_est[i] <- B0
        }
      }
      ### set the quota(s) 
      to_fish <- ifelse((h_crit*B0_est[i]<BF[,i+1,i,2])&(h_crit*B0_est[i]<BF[, i, i,1]),1,0)
      H[i+1,] <- ifelse(h_crit*B0_est[i]<((1-Fmort)*BF[,i+1,i,2]),
                        Fmort*BF[,i+1,i,2],
                        BF[,i+1,i,2]-h_crit*B0_est[i])*to_fish
      if(quota_adj==1){
        red[i] <- sapply(harvest_list,function(x) sum(B[x,i,2]))/sapply(stocklet_list,function(x) sum(B[x,i,2]))
        H[i+1,] <- H[i+1,]*red[i]
      } else {
        red[i] <- 1
      }
    } else {
      if(i==warmup){
        if(assessment== TRUE){
          cat(noquote(paste("iterations 1:",i," of ",n_iter," complete (warmup)\n",sep= "")))
        }
      } 
      H[i+1,] <- rep(0,n_stocks) 
    }
    SBI[i] <- ifelse(H[i]>0,
                     ### if fishery is open and biomass is enough in local area
                     ifelse(B[n_loc, i,1]>0.4*B0/n_loc,1,
                            ### if fishery is open and biomass is enough in neighboring area
                            ifelse((B[n_loc-1, i,1]>0.4*B0/n_loc)|(B[1, i,1]>0.4*B0/n_loc),4,5)),
                     ### if fishery is closed and local biomass is crashed
                     ifelse(B[n_loc, i,1]<0.2*B0/n_loc,
                            ### if fishery is closed and local and far biomass is crashed
                            ifelse((B[n_loc-1, i,1]<0.2*B0/n_loc)&(B[1, i,1]<0.2*B0/n_loc),7,6),
                            ### if fishery is closed, local crash but far biomass is low
                            2))
  }  
  
  ### generate summaries from the simulation 
  
  stock_catch =matrix(sapply(stocklet_list,function(x) rowSums(H_loc[,x])),ncol= n_stocks)
  g_mu_catch =ifelse(n_stocks>1,mean(rowSums(matrix(stock_catch[warmup2:n_iter,]),na.rm= T)),NA)
  s_mu_catch = mean(stock_catch[warmup2:n_iter,],na.rm=T)
  stock_egg_catch =matrix(sapply(stocklet_list,function(x) rowSums(EH_loc[,x])),ncol= n_stocks)
  
  B_A <- rowMeans(apply(tons_at_age*mat,c(2,3),function(x) sum(x[1:3])/sum(x[4:9])))
  
  metrics_df <- data.frame(list(
    pv_cv =sd(B[n_loc,warmup2:n_iter,1])/mean(B[n_loc,warmup2:n_iter,1]),
    po_cv =mean(apply(B[1:(n_loc-1),warmup2:n_iter,1],1,sd)/apply(B[1:(n_loc-1),warmup2:n_iter,1],1,mean)),
    s_cv =mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),1,sd)/apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks) ,1,mean)),
    g_cv =ifelse(n_stocks>1,
                 sd(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,mean))/mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,mean)),
                 NA),
    
    ### mean and equilibrium biomass
    pv_mu =mean(B[n_loc,warmup2:n_iter,1]),
    po_mu =mean(apply(B[1:(n_loc-1),warmup2:n_iter,1],1,mean)),
    s_mu =mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),1,mean))),
    g_mu =ifelse(n_stocks>1,mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,mean)),NA),
    
    s_b0 =mean(B0),
    s_b0_est=mean(B0_est[warmup2:n_iter],na.rm=T),
    p_b0 =mean(B0_loc),
    g_b0 =ifelse(n_stocks>1,sum(B0),NA),
    
    ### overfishing metrics
    pv_of =mean(B[n_loc,warmup2:n_iter,1]<=of_crit*B0/n_loc),
    po_of =mean(apply(B[1:(n_loc-1),warmup2:n_iter,1],2,function(x) x<=of_crit*B0/n_loc)),
    s_of =mean(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,function(x) x<=of_crit*B0)),
    g_of =ifelse(n_stocks>1,
                 mean(!(colMeans(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),2,function(x) x<=of_crit*B0))<1)),
                 NA),
    
    ### closure metrics
    s_closure = mean(stock_catch[warmup2:n_iter,]==0),
    s_closure_l = dep_length(stock_catch[warmup2:n_iter,]==0),
    
    ### mean depletion length
    g_dep_l =dep_length(colMeans(apply(matrix(B_s[,warmup2:n_iter,1]),2,function(x) x<=of_crit*B0)!=0)),
    s_dep_l =mean(apply(apply(matrix(B_s[,warmup2:n_iter,1]),2,function(x) x<=of_crit*B0),1,dep_length),na.rm= T),
    pv_dep_l =dep_length(B[n_loc,warmup2:n_iter,1]<=of_crit*B0_loc[n_loc]),
    po_dep_l =mean(apply(apply(B[1:(n_loc-1),warmup2:n_iter,1],2,function(x) x<=of_crit*B0_loc[1:(n_loc-1)]),1,dep_length),na.rm= T),
    
    ### mean annual catch 
    g_mu_catch = ifelse(n_stocks>1,g_mu_catch,NA),
    s_mu_catch = s_mu_catch,
    s_mu_egg_catch = s_mu_catch,
    s_mu_red = mean(red,na.rm=T),
    ### average annual variation in catch 
    g_aav =ifelse(n_stocks>1,mean(abs(diff(rowSums(H_loc[warmup2:n_iter,])))/rowSums(H_loc[warmup2:n_iter-1,]),na.rm= T)/g_mu_catch,NA),
    s_aav =mean(apply(matrix(stock_catch[warmup2:n_iter,]),2,function(x) abs(diff(x))))/(s_mu_catch),
    p_aav =mean(apply(H_loc[warmup2:n_iter,],2,function(x) diff(x)))/(s_mu_catch/n_loc),
    
    ### average variability in catch 
    g_catch_sd =ifelse(n_stocks>1,sd(rowSums(H_loc[warmup2:n_iter,]),na.rm=T)/g_mu_catch,NA),
    s_catch_sd =mean(apply(matrix(stock_catch[warmup2:n_iter,]),2,function(x) sd(x,na.rm= T)))/(s_mu_catch),
    p_catch_sd =mean(apply(H_loc[warmup2:n_iter,],2,sd,na.rm= T))/(s_mu_catch/n_loc),
    
    ### autocorrelaton in biomass 
    g_phi =ifelse(n_stocks>1,pacf(colMeans(matrix(B_s[,,1],nrow= n_stocks),na.rm=T),na.action=na.pass,plot= FALSE)$acf[1],NA),
    s_phi =mean(apply(matrix(B_s[,,1],nrow= n_stocks),1,function(x) pacf(x,na.action=na.pass,plot= FALSE)$acf[1])),
    pv_phi =pacf(B[n_loc,,1],na.action=na.pass,plot= FALSE)$acf[1],
    po_phi =mean(apply(B[1:(n_loc-1),,1],1,function(x) pacf(x,na.action=na.pass,plot= FALSE)$acf[1])),
    
    ### spatial variance 
    spat_cv =mean(apply(B[,warmup2:n_iter,1],2,function(x) sd(x)/mean(x))),
    
    ### age trunctation
    po_age_trunc = mean(B_A[1:(n_loc-1)],na.rm=T),
    pv_age_trunc = mean(B_A[n_loc],na.rm=T)
  )
  
  ### biomass quantiles
  probs =seq(0.05,0.95,0.05)
  quantile_df <- data.frame(list(
    po_b =rowMeans(apply(B[1:(n_loc-1),warmup2:n_iter,1],1,quantile, probs)),
    pv_b =quantile(B[n_loc,warmup2:n_iter,1],probs),
    s_b = rowMeans(apply(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks),1,quantile, probs)),
    g_b = ifelse(n_stocks>1,quantile(colSums(matrix(B_s[,warmup2:n_iter,1],nrow= n_stocks)), probs),NA),
    p_c = rowMeans(apply(H_loc,2,quantile, probs,na.rm= T)),
    s_catch = rowMeans(apply(stock_catch,2,quantile, probs,na.rm=T)),
    g_c = ifelse(n_stocks>1,quantile(rowSums(stock_catch,na.rm=T), probs,na.rm=T),NA)
  )
  )
  
  ### age depletion quantiles
  ages_dep <- aaply(X[-1,,warmup2:n_iter]*maturity,3,function(x)x/ages0_loc)
  age_trunc <- aaply(X[-1,,warmup2:n_iter]*maturity,3,function(x)x/ages0_loc)
  ages_quants <- apply(apply(ages_dep,c(2,3),quantile,probs),c(1,2),mean)
  ages_quants <- data.frame(rbind(ages_quants,mean=rowMeans(apply(ages_dep,c(2,3),mean))))
  ages_quants$quant <- row.names(ages_quants)
  
  ### estimate Age  for each stock and stocklet
  
  r_sigma <- sd(log(apply(X[2,,warmup2:n_iter],2,sum)))
  
  quantile_df$quants <- row.names(quantile_df)
  if(ret_ts == TRUE){
    TS <- list(B=B,assess= BF,B_stocks= B_s,ages= X,mat=tons_at_age*mat,s_ages = S_freq_s,h_ages =S_freq_f,
               B0=B0,B0_est= B0_est, Harvest= H, forecast= BF, Loc_Harvest = H_loc, SBI=SBI[warmup2:n_iter],red= red,
               sampling= list(catches= catches,catch_comp=catch_comp,age_comp=age_comp,index= index,index.true=index.true))
  } else {
    TS <- NULL
  }
  
  SB_df = as.data.frame.matrix(t(table(SBI[warmup2:n_iter])))
  SB_df = SB_df/rowSums(SB_df)
  
  return(list(metrics= cbind(metrics_df,SB_df),quantiles= quantile_df,age_quantiles = ages_quants,r_sigma= r_sigma, TS=TS))
  
}

