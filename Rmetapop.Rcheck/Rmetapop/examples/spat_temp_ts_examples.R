
### generate the time series
errors <- spat_temp_ts(n_iter = 100, n_loc = 5, 
                       site_sd = 0.5, spat_sd = 2, 
                       phi = 0.5)

# examine correlation matrix of time series
errors$cor_mat

# examine pacf of each time series
errors$pacf

# examine sd of each time series
apply(errors$ts,2,sd)

# plot the time series
matplot(errors$ts,type= "l")
