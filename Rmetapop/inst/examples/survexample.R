
### set up parameters for the simulation test survival & migration function with two
### locations ### general parameters ###
n_loc <- 5  # number of locations
n_iter <- 100  # number of years
n_stages <- 10  # number of stages including eggs and recruits
stage_mat <- 3  # stage # that indicates maturity 

### 9 x L mortality matrix (can also be a single number or a single row with L
### columns)
mort <- matrix(0.334, ncol = n_loc, nrow = length((stage_mat - 1):(n_stages)))

mort_mat <- array(0, dim = c(n_stages - 1, n_stages - 1, n_loc))

### create identifier of matrix entries that correspond to mortality for quick
### replacement in the array
mort_mat_id <- array(matrix(1:n_stages^2, ncol = n_stages) %in% c(diag(matrix(1:n_stages^2, 
                                                                              ncol = n_stages)[-c(1), ]), n_stages^2), dim = c(n_stages, n_stages, n_loc))[-1, 
                                                                                                                                                           -1, ]

### fecundity at age 3-10
fec_at_age <- fecundity_age(stage_mat:(n_stages))

### biomass at age 3-10 in tons
tons_at_age <- LVBweight(stage_mat:(n_stages))/1000

# matrix of mean stray probabilities, controled by a Cauchy distribution with scale = 0.3
spat_scale <- 0.3
stray_probs <- spat_cor_mat(n_loc, spat_scale = spat_scale, sumto1 = TRUE)
stray_probs # columns are probability of coming from, and row is probability going to (i.e. column 2 row 1 is going from location 2 to location 1)
colSums(stray_probs) # check to make sure it sums to 1

### stochastic realizations of the stray matrix for each year using the Dirichlet
### distribution with scale = 1000 (higher = less stochastic probabiliies)
C <- 1000
Crand <- ran_stray_prob(stray_mat=stray_probs,n_iter=n_iter,scale= C)
Crand[,,1:3]  #examine the first three stray matrices in the array

### BH stock-recruit parameters-must be either 1) a single value or 2) a vector of
### length L each
#h = 0.82
# B0 = 60000 
# S0 = B0 * 1000/mean(LVBweight(9))  # convert B0 to mean number of fish
# E0 = S0 * fecundity_age(9)  ### eggs at 6 years old (scaled in millions)
# R0 = 5e+08  # in millions
# alpha <- (E0 * (1 - h))/(4 * h * R0)
# beta <- (5 * h - 1)/(4 * h * R0)
alpha <- 3227
beta <- 1/529032258

### spatiotemporal recruitment deviation parameters
phi = 0.5  ### first order autocorrelation 
spat_sd = 2  ### gaussian spatial correlation paramete
site_sd = 0.2  ### overall log-scale recruitment standard deviation 

### create matrices and vectors from parameters ### create spatiotemporal
### recruitment errors
errors <- spat_temp_ts(n_iter = n_iter, n_loc = n_loc, site_sd = site_sd, spat_sd = spat_sd, 
                       phi = phi)
errors$cor_mat  ### correlation matrix among sites in log-scale recruitment error 
errors$pacf  ### within site partial-autocorrelation function  
error_mat <- errors$ts  ### save the ts for use in the simulation 


### initial conditions ### create the data arrays
d1 <- c(0, 2:10)
d2 <- c(paste("Location", 1:n_loc, sep = "_"))
d3 <- 1:n_iter

## set up the abundance array
X <- array(as.numeric(NA), dim = c(10, n_loc, n_iter), dimnames = list(d1, d2, d3))
X[, , 1:2] <- matrix(rep(1e+08, each = n_loc), ncol = n_loc, byrow = T)

## set up the adult biomass array
B <- array(NA, dim = c(n_loc, n_iter, 3), dimnames = list(d2, d3, c("Actual", "Observed", 
                                                                     "Forecast")))
B[, 1, ] <- tons_at_age %*% X[3:10, , 1]/1000
B[, 2, ] <- tons_at_age %*% X[3:10, , 2]/1000

## set up the Forecast adult biomass array
BF <- array(NA, dim = c(n_loc, n_iter, 3), dimnames = list(d2, d3, c("True", "Observed", 
                                                                      "Estimate")))
BF[, 1, ] <- tons_at_age %*% X[3:10, , 1]/1000
BF[, 2, ] <- tons_at_age %*% X[3:10, , 2]/1000

obs_sd <- 0.3

### project the population with stochastic recruitment ###
system.time(for (i in 3:n_iter) {
  
  mort_mat[mort_mat_id] <- exp(-mort)
  
  X[, , i] <- SSR_linear(alpha=alpha,beta=beta, fec_at_age = fec_at_age, 
                         n_loc = n_loc, n_stages = n_stages,
                         stage_mat = stage_mat,mort_mat = mort_mat,
                         eggs = X[1, , i - 2], X0 = X[(stage_mat - 1):n_stages, , i - 1],
                         stray_mat = Crand[,, i - 1], errors = errors$ts[i - 1, ])

  B[, i, 1] <- tons_at_age %*% X[3:10, , i]/1000
  B[, i, 2] <- B[, i, 1] * exp(rnorm(n_loc, 0, obs_sd) - 0.5 * obs_sd^2)
  B[, i, 3] <- aaply(log(B[, 2:i, 2]), 1, sd1 = obs_sd, sd2 = 0.2, fit = FALSE, 
                     kalman_assess)
})

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

df_ts2 <- as.data.frame.table(B)
names(df_ts2) <- c("Location", "Time", "Source", "Tons")
df_ts2$Time <- as.numeric(as.character(df_ts2$Time))
### generate the plot ###

ggplot(aes(Time, Tons), data = df_ts2) + 
  geom_line(aes(colour = Source, linetype = Source),  size = 1) + 
  scale_linetype_manual(values = c(1, 3, 1)) + 
  plot.options + facet_grid(Location ~ .) + 
  scale_colour_manual(values = c("red", "black", "blue")) + 
  ylab("Tons (Thousands)") + 
  scale_x_continuous(expand = c(0, 0))
### correlation in recruitment among sites
cor(t(X[2, , ]))

### correlation in biomass among sites
cor(t(B[, , 1]))



