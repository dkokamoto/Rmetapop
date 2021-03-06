
fish_time <- function (time, state, pars, N) {
  state_array <- array(state, dim = c(N,3))
  Fish    <-  state_array[,1]
  Effort  <-  state_array[,2]
  Harvest <-  state_array[,3]
  
  with (as.list(pars), {
    ## fishing equations
    qt <- ifelse(sum(Harvest)<h_max,q,0)
    dFish    <- - Effort * Fish * qt
    dEffort  <- 0
    dHarvest <- Effort * Fish * qt
    
    ## flux 
    FluxEffort    <-  Effort*De *rowSums(outer(Fish,Fish,"-"))
    
    ## Add flux gradient to rate of change
    dEffort       <- dEffort + FluxEffort
    
    return(list(c(as.vector(dFish), as.vector(dEffort),as.vector(dHarvest))))
  })
}


N  <- 10
max_time <- 10
n_times <- 100

set_effort <- function(quota, biomass,catchability=0.05,diffusion= 10000000000,N=10,n_times= 100,max_time=2){
    
    ## ===================
    ## Apply the model 
    ## =================== 
    scalar <- sum(biomass)/N
    pars    <- c(h_max  = quota/scalar,q= catchability,De=diffusion)   
    ## array of initial conditions   
    yini_array <- array(0,dim=c(N,3))
    yini_array[,1] <- biomass/scalar
    yini_array[,2] <- 1
    yini_array[,3] <-0
    yini <- as.vector(yini_array)
    
    ## solve model (use Cash-Karp Runge-Kutta method
    times   <- seq(from=0, to=max_time, length.out= n_times)
    
    solution <- ode.1D(y = yini, times = times, func = fish_time, parms = pars,
                   dimens = N, names = c("Fish", "Effort","Harvest"),
                   N = N, method = rkMethod("rk4"))
    a <- melt(solution[,12:21])
    ggplot(aes(Var1,value),data= a)+geom_line(aes(colour= factor(Var2)))
    
    return(list(Harvest= solution[n_times,(2*N+2):(3*N+1)]*scalar, CPUE = solution[n_times,(2*N+2):(3*N+1)]/solution[n_times,(N+2):(2*N+1)]))
}

