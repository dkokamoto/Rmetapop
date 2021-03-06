
### set up parameters for the simulation test survival & migration function with two
### locations ### general parameters ###
E0 <- 3355331342
R0 <- 527.91*1e6
h <- 0.83
M <- 0.47

### convert steepness, R0 and E0 to Beverton-Holt alpha and beta parameters
a_bh <- (E0  * (1 - h))/(4 * h  * R0)
b_bh <- (5 * h - 1)/(4 * h  * R0)

### life history characteristics
stage_mat <- 3          ### Stage at maturity
n_stages <- 10          ### Plus group

### spatial characteristics 
n_loc <- 5
C <- 1000

### straying and stochastic info 
spat_scale <- 0.5  ### Cauchy straying scale
site_sd <-0.5       ### spatial standard deviation in recruitment correlations
surv_rho  <- 0.05       ### scale of stochasticity in survival
rec_sd  <- 0.6          ### lognormal recruitment sd
phi  <- 0.5             ### AR(1) autoregressive coefficient in recruitment 
obs_sd <- 0.3           ### Lognormal survey error 

stock_IDs= rep(1,each= 5)
Nstocks <- 1

n_iter <- 100
### 9 x L mortality matrix (can also be a single number or a single row with L
### columns) in this case randomly drawn from mean M and sd 0.25

mort <- matrix(rnorm(n_loc,M,0.25), ncol = n_loc, nrow = length((stage_mat - 1):(n_stages)),byrow= T)

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

### create matrices and vectors from parameters ### create spatiotemporal
### recruitment errors
errors <- spat_temp_ts(n_iter = n_iter, n_loc = n_loc, site_sd = site_sd, spat_sd = .1, 
                       phi = phi,cor_mat=NULL)

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

### set up a harvest array (zero harvest)
H <- array(NA, dim = c(n_iter+1,Nstocks), dimnames = list(1:(n_iter+1),d2b))

for (i in 1:n_iter){
  H[i,]  <- 0
}

### project the population with stochastic recruitment and harvesting ###
for (i in 3:(n_iter)) {
  surv_array[surv_array_id] <- ran_surv_prob(mean=exp(-mort),corr=surv_rho)
  X[, , i] <- ssr_linear(alpha=a_bh,beta=b_bh, fec_at_age = fec_at_age, 
                         n_loc = n_loc, n_stages = n_stages,
                         N_s=Nstocks,s_ID=stock_IDs,
                         tons_at_age=tons_at_age,harvest=H[i,],
                         stage_maturity = stage_mat,
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
ggplot(aes(Time, Number/1e+06), data = df_ts) + 
  geom_line(aes(colour = Age_Class)) + 
  facet_grid(Stage ~ Location,scales= "free_y") + 
  plot.options + 
  ylab("Number (Millions)")+
  scale_x_continuous(expand = c(0, 0))

### correlation in recruitment among sites
levelplot(cor(t(X[2, ,1:(n_iter-1)])))

### correlation in biomass among sites
levelplot(cor(t(B[,1:(n_iter-1)])))


