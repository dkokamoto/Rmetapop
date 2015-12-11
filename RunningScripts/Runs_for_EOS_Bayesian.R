
source("http://mc-stan.org/rstan/stan.R")
### stock assessment parameters (Table E7- 2014/15 assessment)
E0_assess <- 3355331342
R0_assess <- 527.91*1e6
h_assess <- 0.83
M_assess <- 0.47
stage_mat <- 3
a_bh <- (E0_assess  * (1 - h_assess))/(4 * h_assess  * R0_assess )
b_bh <- (5 * h_assess - 1)/(4 * h_assess  * R0_assess )

stray_s_scale <- c(0.01,0.05,0.1,0.5,1,2,5,10)
rec_corr_sd <- c(0.01,0.05,0.1,0.5,1,2,5,10)
surv_rho  <- c(0.05)
rec_sd  <- c(0.6)
phi  <- c(0.5)
Fmort <- c(0,0.2,0.4)

### warmup and number of simulations###
warmup1 <- 15
warmup2 <-32
totiter <- 72

### dataframe of all parameter combinations for the simulations ###
param.df <- expand.grid(stray_s_scale= stray_s_scale,rec_corr_sd=rec_corr_sd,
                        surv_rho= surv_rho, rec_sd=rec_sd,phi=phi,Fmort= Fmort)

rows.to.list <- function( df ) {
  ll<-apply(df,1,list)
  ll<-lapply(ll,unlist)
}

### list of parameters for parallel runs
param.list <- rows.to.list(param.df)  

### complile the model so it need not be directly compiled each time ###
compiled <-fishery_simulate(n_loc=50,n_stages=10,stage_mat=3,
                         a_bh=a_bh, b_bh=b_bh, phi=0.5,M=M_assess, spat_scale=0.5,
                         spat_sd=0.0001,surv_rho = 0.05,
                         site_sd=0.2,C=1000,obs_sd=0.3,stock_IDs= rep(1:5,each= 10), 
                         n_iter=3,warmup= 2,point.estimate=FALSE,
                         h_crit=0.2, Fmort=0.2)
             
parallel.sim <- function(x){  
    x <- as.data.frame(t(x))
    
    assess <- fishery_simulate(n_loc=50,                     ### total number of subpopulations
                               n_stages=10,                  ### number of age/stage classes
                               stage_mat=3,                  ### stage/age at maturity
                               a_bh=a_bh,                    ### BH alpha parameter
                               b_bh=b_bh,                    ### BH beta parameter
                               phi=x$phi,                    ### lag 1 recruitment autocorrelation 
                               M=M_assess,                   ### natural mortality (a constant, vector or matrix)
                               spat_scale=x$stray_s_scale,   ### spatial Cauchy straying scale
                               spat_sd=x$rec_corr_sd,        ### spatial SD of recruitment synchrony
                               surv_rho=x$surv_rho,          ### degree of stochasticity in survival 
                               site_sd=x$rec_sd,             ### degree of stochasticity in alpha
                               C=1000,                       ### degree of stochasticity in straying (lower = more stochastic)
                               obs_sd=0.3,                   ### observation error on surveys of total biomass
                               stock_IDs= rep(1:5,each= 10), ### ids representing to which stock the populations belong
                               n_iter=totiter,               ### total number of iterations
                               warmup= warmup1,              ### warmup prior to harvest
                               point.estimate=FALSE,         ### assessment point estimate vs full posterior
                               h_crit=0.25,                  ### lower harvest cutoff
                               Fmort=x$Fmort                 ### proportional fishing rate
                               )
    
    return(list(B = assess$B, B_stocks = assess$B_stocks,Eggs= assess$ages[1,,],
                assessment = assess$assess, Harvest= assess$Harvest, 
                cors= assess$corrs, B0= assess$B0))
}

ptm <- proc.time()
for(i in 1:10){
  fit1 <- mclapply(param.list,parallel.sim,mc.cores= 8)  
  save(fit1,file= paste0("Bayes_results","_",Sys.Date(),"_",i,".RData"))
}
proc.time() - ptm
