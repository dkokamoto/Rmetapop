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


fishery_simulate <- function(n_loc,stock_IDs,n_stages,stage_mat,#n_subloc,
         alpha, beta, phi,M, spat_scale,spat_sd,surv_rho,
         site_sd,C,obs_sd,cor_mat=NULL,
         n_iter,point.estimate,warmup,h_crit){

Nstocks <- length(unique(stock_IDs))
  
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
round(stray_probs,2)# columns are probability of coming from, and row is probability going to (i.e. column 2 row 1 is going from location 2 to location 1)
colSums(stray_probs) # check to make sure it sums to 1

### stochastic realizations of the stray matrix for each year using the Dirichlet
### distribution with scale = 1000 (higher = less stochastic probabiliies)
Crand <- ran_stray_prob(stray_mat=stray_probs,n_iter=n_iter,scale= C)
Crand[,,1:3]  #examine the first three stray matrices in the array


### create matrices and vectors from parameters ### create spatiotemporal
### recruitment errors
errors <- spat_temp_ts(n_iter = n_iter, n_loc = n_loc, site_sd = site_sd, spat_sd = .1, 
                       phi = phi,cor_mat)

errors$cor_mat  ### correlation matrix among sites in log-scale recruitment error 
errors$pacf  ### within site partial-autocorrelation function  
error_mat <- errors$ts  ### save the ts for use in the simulation 

### initial conditions ### create the data arrays
d1 <- c(0, 2:n_stages)
d2 <- paste("Location", 1:n_loc, sep = "_")
d2b <- paste("Stock",1:Nstocks,sep= "_")
d3 <- 1:n_iter

## set up the abundance array
X <- array(as.numeric(NA), dim = c(n_stages, n_loc, n_iter), dimnames = list(d1, d2, d3))
X[, , 1:2] <- matrix(rep(1e+08, each = n_loc), ncol = n_loc, byrow = T)

## set up the adult biomass array
B <- array(NA, dim = c(n_loc, n_iter+1), dimnames = list(d2, 1:(n_iter+1)))
B_s <- array(NA, dim = c(Nstocks, n_iter+1,2), dimnames = list(d2b, 1:(n_iter+1),c("actual","observed")))
B[, 1] <- tons_at_age %*% X[3:10,, 1]
B[, 2] <- tons_at_age %*% X[3:10,, 2]

for(i in 1:Nstocks){
  B_s[i, 1,] <- sum(B[stock_IDs==i, 1])
  B_s[i, 2,] <- sum(B[stock_IDs==i, 2])
}

## set up the Forecast adult biomass array
BF <- array(NA, dim = c(Nstocks, n_iter+1, n_iter,6), dimnames = list(d2b, 1:(n_iter+1), d3,c("Estimated Stock Size","Prediction","ESS_UP","ESS_LOW","P_UP","P_LOW")))

## set up the age frequnecy arrays
S_freq <- array(NA, dim = c(n_stages-2,n_loc, n_iter), dimnames = list(d1[-c(1,2)],d2,d3))
S_freq_s <- array(NA, dim = c(n_stages-2,Nstocks, n_iter), dimnames = list(d1[-c(1,2)],d2b,d3))

for (i in 1:2){
  S_freq[,,i]  <- t(t(X[-c(1,2),,i])/colSums(X[-c(1,2),,i]))
  SN <- mapply(function(x) rowSums(X[-c(1,2),stock_IDs==x,i]),1:Nstocks)
  S_freq_s[,,i]  <- t(t(SN)/colSums(SN))
}

H <- array(NA, dim = c(n_iter+1,Nstocks), dimnames = list(1:(n_iter+1),d2b))

for (i in 1:n_iter){
  H[i,]  <- 0
}
  for (i in 3:n_iter) {
    surv_array[surv_array_id] <- ran_surv_prob(mean=exp(-mort),corr=0.05)
    X[, , i] <- ssr_linear(alpha=alpha,beta=beta, fec_at_age = fec_at_age, 
                         n_loc = n_loc, n_stages = n_stages,harvest=rep(0,Nstocks),
                         tons_at_age=tons_at_age,DI_selectivity= rep(1,length(stage_mat:n_stages)),
                         stage_mat = stage_mat,surv_array = surv_array,
                         s_ID=stock_IDs,N_s=Nstocks,
                         eggs = X[1, , i - 2], X0 = X[(stage_mat - 1):n_stages, , i - 1],
                         stray_mat = stray_probs, errors = 0)
    
    B[, i] <- tons_at_age %*% X[3:10, , i]
    B_s[, i,] <- mapply(function(x) sum(B[stock_IDs==x,i]),1:Nstocks) 
    ### age frequency ###
    SN <- mapply(function(x) rowSums(X[-c(1,2),stock_IDs==x,i]),1:Nstocks)
    S_freq[,,i]  <- t(t(X[-c(1,2),,i])/colSums(X[-c(1,2),,i]))
    S_freq_s[,,i]  <- t(t(SN)/colSums(SN))
  } 

  B0 <- tons_at_age %*% mapply(function(x) rowSums(X[-c(1,2),stock_IDs==x,n_iter]),1:Nstocks)
  B0_loc <- tons_at_age %*% X[-c(1,2),,n_iter]

fit.compile <- nlss_assess(b_obs=array(B0,dim= c(Nstocks,n_iter)),
                          b_true=array(B0,dim= c(Nstocks,n_iter)),
                          age_freq=S_freq_s[,,1:n_iter],obs_sd=obs_sd,sd_R=site_sd,
                          alpha_BH= alpha,beta_BH=beta/(n_loc/Nstocks),harvest= H[1:n_iter,],
                          selectivity=rep(1,length(stage_mat:n_stages)),n_loc=Nstocks,
                          fec_at_age=fec_at_age,weight_at_age=tons_at_age,
                          mort=rep(M,Nstocks),stage_mat=stage_mat,
                          n_stages=n_stages,phi=phi,n_warmup=3,
                          compile=TRUE,optimize=point.estimate,n_iter=n_iter)

### reset initial values as last two values of deterministic simulation 

B[, 1] <- B[, n_iter-1]
B[, 2] <- B[, n_iter]
B_s[, 1,] <- B_s[, n_iter-1,]
B_s[, 2,] <- B_s[, n_iter,]
for (i in 1:2){
  S_freq[,,i]  <- S_freq[,,n_iter-2+i]
  S_freq_s[,,i]  <- S_freq_s[,,n_iter-2+i]
}

### project the population with stochastic recruitment and harvesting ###
for (i in 3:(n_iter)) {
  ptm <- proc.time()
  surv_array[surv_array_id] <- ran_surv_prob(mean=exp(-mort),corr=0.05)
  X[, , i] <- ssr_linear(alpha=alpha,beta=beta, fec_at_age = fec_at_age, 
                         n_loc = n_loc, n_stages = n_stages,
                         N_s=Nstocks,s_ID=stock_IDs,
                         tons_at_age=tons_at_age,harvest=H[i,],
                         stage_mat = stage_mat,
                         surv_array = surv_array,
                         DI_selectivity= rep(1,length(stage_mat:n_stages)),
                         eggs = X[1, , i - 2], 
                         X0 = X[(stage_mat - 1):n_stages, , i - 1],
                         stray_mat = Crand[,, i - 1], errors = errors$ts[i - 1, ])

  B[, i] <- tons_at_age %*% X[3:10, , i]
  B_s[, i,1] <-  mapply(function(x) sum(B[stock_IDs==x,i]),1:Nstocks) 
  B_s[, i,2] <- B_s[, i,1]* exp(rnorm(Nstocks, 0, obs_sd) - 0.5 * obs_sd^2)
  
  ### age frequency ###
  SN <- mapply(function(x) rowSums(X[-c(1,2),stock_IDs==x,i]),1:Nstocks)
  S_freq[,,i]  <- t(t(X[-c(1,2),,i])/colSums(X[-c(1,2),,i]))
  S_freq_s[,,i]  <- t(t(SN)/colSums(SN))
 
  ### stock assessment
    if(i>warmup){  
    assess <- nlss_assess(b_obs=B_s[, 1:i,2],b_true=B_s[, 1:i,1],sd_R=site_sd,
                  age_freq=S_freq_s[,,1:i],obs_sd=obs_sd,
                  alpha_BH= alpha,beta_BH=I(beta/(n_loc/Nstocks)),harvest= H[1:i,],
                  selectivity=rep(1,length(stage_mat:n_stages)),n_loc=Nstocks,
                  fec_at_age=fec_at_age,weight_at_age=tons_at_age,
                  mort=rep(M,Nstocks),stage_mat=stage_mat,
                  n_stages=n_stages,phi=phi,n_warmup=3,compile= FALSE,
                  n_iter=i,optimize=point.estimate,compiled_model= fit.compile)
     
    if(point.estimate==TRUE){
       BF[, 1:(i), i,1] <- t(assess$assess$mean)
       BF[, 1:(i+1), i,2] <- t(assess$assess$pred)
     }
    else{
      BF[, 1:(i), i,1] <- t(assess$assess$mean)
      BF[, 1:(i+1), i,2] <- t(assess$assess$pred)
      BF[, 1:(i), i,3] <- t(assess$assess$upper)
      BF[, 1:(i), i,4] <- t(assess$assess$lower)
      BF[, 1:(i+1), i,5] <- t(assess$assess$pred_upper)
      BF[, 1:(i+1), i,6] <- t(assess$assess$pred_lower)
    }
     ### forecast array 
     Hcrit <- ifelse(h_crit*B0<BF[,i+1,i,2],1,0)
     H[i+1,] <- mapply(min,BF[,i+1,i,2]-h_crit*B0,0.2*BF[,i+1,i,2])*Hcrit 
   cat(noquote(paste("iteration",i, "of", n_iter, "completed in", signif(proc.time()[3]-ptm[3],digits=4),"sec\n")))
  }
  else{
    H[i+1,] <- rep(0,Nstocks)
    if(i==warmup){
      cat(noquote(paste("iterations 1:",i," of ",n_iter," complete (warmup)\n",sep= "")))
    }
  }
}

B0.df <- melt(B0,varnames=c("time","location"),value.name= "B0") 

pred <- melt(apply(BF[,-1,,2],1,diag),varnames=c("time","site"),value.name= "pred")
mean <- melt(apply(BF[,-n_iter,,1],1,diag),varnames=c("time","site"),value.name= "mean")
forecast_array <- join(pred,mean,by =c("time","site"))

assess2 <- reshape(melt(assess$assess,varnames= c("time","location")), timevar = "L1", idvar = c("time","location"), direction = "wide")
names(assess2) <- gsub("value.","",names(assess2))

if(point.estimate==TRUE){
plot1 <- ggplot(aes(time,observed),data= assess2)+
    geom_line(aes(y=mean),col= "black")+
    geom_line(aes(y=pred),col= "blue")+
    geom_point(shape= 21, fill= "black",size= 0.5)+
    geom_point(aes(y=actual),shape= 21, fill= "pink",col= "red",size =0.5)+
    facet_wrap(~location)+
    geom_hline(aes(yintercept= B0),data= B0.df,linetype="dotted")+
    geom_hline(aes(yintercept= 0.25*B0),data= B0.df,col= "red")+
    theme_bw()+geom_vline(xintercept= warmup+1)
}
else{
plot1 <- ggplot(aes(time,observed),data= assess2)+
    geom_ribbon(aes(ymin= lower,ymax= upper),col= "grey70",fill= "grey70")+
    geom_line(aes(y=mean),col= "black")+
    #geom_ribbon(aes(ymin= pred_lower,ymax= pred_upper),col= "blue",fill= NA)+
    geom_line(aes(y=pred),col= "blue")+
    geom_point(shape= 21, fill= "black",size= 0.5)+
    geom_point(aes(y=actual),shape= 21, fill= "pink",col= "red",size =0.5)+
    facet_wrap(~location)+
    geom_hline(aes(yintercept= B0),data= B0.df,linetype="dotted")+
    geom_hline(aes(yintercept= 0.25*B0),data= B0.df,col= "red")+
    theme_bw()+geom_vline(xintercept= warmup+1)
}

### create a dataframe for plotting
df_ts <- as.data.frame.table(X)
names(df_ts) <- c("Age", "Location", "Time", "Number")
df_ts$Time <- as.numeric(as.character(df_ts$Time))
df_ts$Age <- as.numeric(as.character(df_ts$Age))
df_ts$Stage <- with(df_ts, ifelse(Age == 2, "Recruits", ifelse(Age == 0, "Eggs", "Adults")))
df_ts$Age_Class <- with(df_ts, ifelse(Age == 2, "Recruits", ifelse(Age == 0, "Eggs",ifelse(Age>9,"9+", Age))))


### plot the time series of adults plotting options
plot.options <- theme(strip.background = element_rect(fill = NA, colour = NA), 
                      panel.grid.minor = element_line(colour = NA), 
                      panel.grid.major = element_line(colour = NA), 
                      legend.key = element_blank(), 
                      legend.key.size = unit(0.6, "cm"), legend.position = "right", 
                      legend.direction = "vertical", 
                      panel.background = element_rect(fill = NA, colour = "black"),
                      legend.key.width = unit(1.5, "lines"), 
                      panel.margin = unit(1,  "lines"),
                      legend.background = element_rect(fill = NA, colour = NA),
                      plot.background = element_rect(fill = NA), 
                      plot.margin = unit(c(0.075, 0.075, 0.075, 0.2), "inches"))

### generate the plot ###
plot2 <- ggplot(aes(Time, Number/1e+06), data = df_ts) + 
  geom_line(aes(colour = Age_Class)) + 
  facet_grid(Stage ~ Location,scales= "free_y") + 
  plot.options + 
  ylab("Number (Millions)")+
  scale_x_continuous(expand = c(0, 0))

### correlation in recruitment among sites
rec_cor <- cor(t(X[2, ,warmup:n_iter ]))

### correlation in biomass among sites
B_cor <- cor(t(B[, warmup:n_iter]))



# stats <- data.frame(
#   list(
#     "B_hat_lt_h_crit"=c(mean(apply(B_s[,(warmup+2):n_iter,1],2,function(x)h_crit*t(B0)/x)>1),
#                         mean(colMeans(apply(B_s[,(warmup+2):n_iter,1],2,function(x)h_crit*t(B0)/x)>1)==1)),
#     "B_lt_h_crit"=c(mean(H[(warmup+2):n_iter,]==0),
#                     mean(rowMeans(H[(warmup+2):n_iter,]==0)==1)),
#     "Bloc_lt_h_crit"=c(mean(c(test$B[,(warmup+2):n_iter]<=0.25*B0_loc)),
#                        mean(rowMeans(test$H[(warmup+2):n_iter,]==0)==1)),
#     "CV" = c(mean(apply(H[(warmup+2):n_iter,],2,function(x)sd(x)/mean(x))),
#              sd(rowSums(H[(warmup+2):n_iter,]))/mean(rowSums(H[(warmup+2):n_iter,]))),
#     source= c("Local","Global")))


return(list(B=B,asess=BF,B_stocks=B_s,ages=X,B0=B0.df,Harvest=H,forcast=forecast_array,plot_biomass=plot1,plot_nums=plot2,corrs= list(recruitment=rec_cor,biomass=B_cor,B0_loc)))
}
