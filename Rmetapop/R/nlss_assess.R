#' Nonlinear state-space model stock assessment.
#' @param b_obs
#' @param b_true
#' @param age_freq
#' @param obs_sd
#' @param alpha_BH
#' @param sd_R,
#' @param beta_BH
#' @param fec_at_age
#' @param weight_at_age
#' @param mort
#' @param stage_mat
#' @param n_stages
#' @param harvest
#' @param selectivity
#' @param n_loc
#' @param n_iter
#' @param phi
#' @param plot
#' @param n_warmup
#' @param MCMC.iter
#' @param forecast_limit
#' @param optimize
#' @param compile
#' @param compiled_model
#' @param init_states
#' @param sd_pro

### Bayesian State-Space Model using Stan
nlss_assess <- function(b_obs,b_true,age_freq,obs_sd,alpha_BH,sd_R,
                       beta_BH,fec_at_age,weight_at_age,mort,
                       stage_mat,n_stages,
                       harvest,selectivity,
                       n_loc,n_iter,phi,
                       plot= TRUE,
                       n_warmup=3,MCMC.iter=200,forecast_limit=0.95,
                       optimize=TRUE,compile=TRUE,compiled_model=NULL,
                       init_states= NULL,sd_pro= NULL){
  
  ### get biomass to abundance at age conversion given size frequency and size at age

  ### generate a list of data for the 
  model_data <- list(B_obs= t(b_obs),
             sigma_obs=obs_sd,
             alpha=alpha_BH,
             beta=beta_BH,
             fec_aa=fec_at_age,
             w_aa=weight_at_age,
             BpN=aperm(apply(age_freq,c(2,3),function(x) 
               (1/x%*%((weight_at_age)))*x), c(3,2,1)),
             surv=exp(-mort),
             sigma_R=sd_R,
             n_ages = length(stage_mat:n_stages),
             n_i= n_iter,
             n_l= n_loc,
             H=harvest,
             S=selectivity,
             n_w = 3,
             phi= phi,
             offset = sd_R^2/(1-phi^2)/2)

  ### list of parameters to save ###
  params<- c("sigma_pro","B_hat","B")
  
  ### initial values ###
  inits <- function(){list(
    B_hat=  rbind(init_states,init_states[nrow(init_states),]), 
    sigma_pro = 0.4,
    epsilon = matrix(rep(0,(n_iter+1-3)*n_loc),ncol= n_loc))}
  
  
  ### MCMC details
  MCMC.iter <- 200
  n.burnin <- 50
  set.seed=134
  n.chains <- 1

model <- '
  data {
    int n_l; // number of locations
    int n_ages; // number of ages
    int n_i; // number of iterations
    int n_w;// number of warmup iterations
    row_vector<lower=0.000001>[n_ages] fec_aa; // fecundity at age
    row_vector[n_ages] w_aa; // weight at age
    vector<lower=0>[n_l] B_obs[n_i]; // observed biomass
    real sigma_obs; // observation error for biomass
    vector[n_ages] BpN[n_i,n_l]; //biomass per number
    vector[n_l] H[n_i]; //harvest
    vector[n_ages] S; //selectivity
    real alpha; // alpha
    real beta; // beta
    vector[n_l] surv; //survival
    real phi;
    real sigma_R;
    real offset;
  }

  parameters {
    real<lower= 0.000001,upper=5> sigma_pro; // process error variance
    vector<lower= 0.000001,upper=300>[n_l] B_hat[n_i];
    real<lower=-4,upper=4> epsilon[n_i+1-n_w,n_l];
  }

  transformed parameters {
    vector<lower= 0.0000000000001>[n_ages] N[n_i+1,n_l];
    vector<lower= 0.0000000000001>[n_l] B[n_i+1];
    vector[n_l] R[n_i+1];
    vector[n_l*(n_i)] ln_B;
    vector[n_l*n_i] ln_Bhat;
    vector[n_l*n_i] ln_Bobs;
    #vector[n_l] surv;
    #real offset;
    real eps_phi[n_i+1-n_w,n_l];  
    
    eps_phi[1] <- epsilon[1];
    for (e in 2:(n_i+1-n_w)){
      for(j in 1:n_l){
       eps_phi[e,j] <- epsilon[e,j]+phi*eps_phi[e-1,j];
      }
    }
 
  ### initial 
  for(i in 1:n_w){  
    for (j in 1:n_l){
      B[i,j] <- B_obs[i,j];
      N[i,j] <- BpN[i,j]*B[i,j];
      R[i,j] <-  N[i,j,1]; 
      ln_Bhat[i+(j-1)*n_i] <- log(B_hat[i,j]);
      ln_Bobs[i+(j-1)*n_i] <- log(B_obs[i,j]);
    }
  }
    for(i in (n_w+1):n_i){
      for (j in 1:n_l){
        ### number of three year olds ###
        R[i,j] <- (fec_aa*(BpN[i-3,j]*B_hat[i-3,j]))/(alpha + beta * (fec_aa*(BpN[i-3,j]*B_hat[i-3,j])));
        N[i,j,1] <- R[i,j]*exp(eps_phi[i-n_w+1,j]-offset)*surv[j]; 
        ### number of mature fish ###
        for (k in 2:(n_ages-1)){
          N[i,j,k] <- BpN[i-1,j,k-1]*B_hat[i-1,j]*surv[j];
        }
        ### terminal age class
        N[i,j,n_ages] <- BpN[i-1,j,n_ages-1]*B_hat[i-1,j]*surv[j]+
                            BpN[i-1,j,n_ages]*B_hat[i-1,j]*surv[j];
        ### final biomass
        if  (w_aa*N[i,j]-H[i,j] >0)
          B[i,j] <- w_aa*N[i,j]-H[i,j];
        else
           B[i,j] <- 0.00001;
      }
    }
   ### log-transformed values
   for(i in 1:n_i){  
    for (j in 1:n_l){
        ln_B[i+(j-1)*n_i] <- log(B[i,j]);
        ln_Bhat[i+(j-1)*n_i] <- log(B_hat[i,j]);
        ln_Bobs[i+(j-1)*n_i] <- log(B_obs[i,j]);
    }
  }

  ### one-year ahead forecast given the model 
   for (j in 1:n_l){
        ### number of three year olds ###
        R[n_i+1,j] <- (fec_aa*(BpN[n_i-2,j]*B_hat[n_i-2,j]))/(alpha + beta * (fec_aa*(BpN[n_i-2,j]*B_hat[n_i-2,j])))*surv[j];;
        N[n_i+1,j,1] <- R[n_i+1,j]*exp(eps_phi[n_i+1-n_w,j]-offset)*surv[j];     
        for (k in 2:(n_ages-1)){
          N[n_i+1,j,k] <- BpN[n_i,j,k-1]*B_hat[n_i,j]*surv[j];
        }
        N[n_i+1,j,n_ages] <- BpN[n_i,j,n_ages-1]*B_hat[n_i,j]*surv[j]+
                            BpN[n_i,j,n_ages]*B_hat[n_i,j]*surv[j];        
        B[n_i+1,j] <- w_aa*N[n_i+1,j];
  }
}

  model { 
        for (e in 1:(n_i+1-n_w)){
          for (j in 1:n_l){
            epsilon[e,j] ~ normal(0,sigma_R);//stochastic recruitment variation
          }
        }
        ln_Bhat ~ normal(ln_B,sigma_pro); // process equation
        ln_Bobs ~ normal(ln_Bhat,sigma_obs); // observation equation
  }
'
if(optimize==TRUE){
  if(compile==TRUE){
    opt_model<-stan_model(model_code=model)
    return(opt_model)
  }
  else{
    fit.assess <- invisible(optimizing(compiled_model,data= model_data,iter=5000)) 
    Bhat_m<- matrix(fit.assess$par[grep("B_hat",names(fit.assess$par))],,nrow= n_iter)
    B_m <-matrix(fit.assess$par[grep("B\\[([[:digit:]]+),([[:digit:]]+)]",names(fit.assess$par))],nrow= n_iter+1)
    assessment = list(mean=Bhat_m,observed=t(b_obs),
                      actual = t(b_true),pred=B_m)
    return(list(fit=fit.assess,assess=assessment))
  }
}
else{
 if(compile==TRUE){
  fit.assess <- invisible(stan(model_code= model,data= model_data,pars= params,
                   iter= 2,warmup= 1,chains=1,
                   init= "random",seed= set.seed,refresh=0))
 }
 else{
   if(is.null(init_states)){
    fit.assess <- invisible(stan(model_code= model,data= model_data,pars= params,
                    iter= MCMC.iter,warmup= n.burnin,chains=n.chains,
                    init= "random",seed= set.seed,refresh=0))
  }
  else{
    fit.assess <- invisible(stan(model_code= model,data= model_data,pars= params,
                                 iter= MCMC.iter,warmup= n.burnin,chains=n.chains,
                                 init= inits,seed= set.seed,refresh=0))
  }
  MCMC <- extract(fit.assess)
  sd_pro <- mean(MCMC$sigma_pro)
  Bhat_m <- apply(MCMC$B_hat,c(2,3),mean)
  Bhat_HPD <- apply(MCMC$B_hat,c(2,3),
                    function(x) HPDinterval(as.mcmc(x),prob= forecast_limit))
  
  B_m <- apply(MCMC$B,c(2,3),mean)
  B_HPD <- apply(MCMC$B,c(2,3),
                 function(x) HPDinterval(as.mcmc(x),prob= forecast_limit))
  
  assessment = list(mean=Bhat_m,upper=Bhat_HPD[2,,],lower=Bhat_HPD[1,,],observed=t(b_obs),actual = t(b_true),pred=B_m,pred_upper= B_HPD[2,,],pred_lower= B_HPD[1,,],sd_pro= sd_pro)
      return(list(fit=fit.assess,assess=assessment))
    }
  }
}
