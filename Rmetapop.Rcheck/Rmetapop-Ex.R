pkgname <- "Rmetapop"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Rmetapop')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("BH")
### * BH

flush(stderr()); flush(stdout())

### Name: BH
### Title: Beverton-Holt stock recruit relationship.
### Aliases: BH

### ** Examples

h = 0.82                              # steepness of BC the central coast herring fishery
B0 = 60000                            # tons in central coast in 2013/14
S0 = B0 * 1000/mean(LVBweight(3:10))  # convert B0 to mean number of fish
E0 = S0 * fecundity_age(9)            # eggs at 6 years old (scaled in millions)
R0 = 5e+08                            # in millions
alpha <- (E0 * (1 - h))/(4 * h * R0)
beta <- (5 * h - 1)/(4 * h * R0)

sim_eggs <- sort(c(seq(from = 1, to = E0*1.02, length.out = 1000),E0/5,E0))   # project 
R <- BH(sim_eggs, alpha=alpha,beta= beta)
df1 <- data.frame(list(Recruits=R/1e6,Eggs=sim_eggs/1e12))
df2 <- subset(df1,Eggs%in%(c(E0/5,E0)/1e12))

ggplot(aes(Eggs,Recruits),data= df1)+
  geom_line()+
  theme_bw()+
  ylab("Recruits (Millions)")+
  xlab("Eggs (Trillions)")+
  geom_vline(aes(xintercept = Eggs),data= df2,linetype= "dotted")+
  geom_hline(aes(yintercept = Recruits),data= df2,linetype= "dotted")+
  geom_text(aes(x=Eggs+(E0/1e12)*.02, y= 1/beta/1e6*.25,label = c("20% E0", "E0")),data= df2,hjust= 0,vjust= -1,angle=90)+
  geom_line(yintercept= 1/beta/1e6)+
  geom_errorbar(aes(x=E0/1e12,ymax=max(Recruits),ymin= min(Recruits)),data= subset(df1,Eggs%in%(c(E0/5,E0)/1e12)),colour= "red")+
  geom_errorbar(aes(x=E0/1e12,ymax=min(Recruits),ymin= 0),data= subset(df1,Eggs%in%(c(E0/5,E0)/1e12)),colour= "blue")+
  geom_text(x= E0/1e12*.7,y= 1/beta/1e6*.2,label= "h (steepness)= blue/(red+blue)",data= df1[1,])+
  scale_x_continuous(limits=c(0,E0/1e12*1.02),expand= c(0,0))+
  scale_y_continuous(limits=c(0,1/beta/1e6*1.02),expand= c(0,0))




cleanEx()
nameEx("Rmetapop-package")
### * Rmetapop-package

flush(stderr()); flush(stdout())

### Name: Rmetapop-package
### Title: What the package does (short line) ~~ package title ~~
### Aliases: Rmetapop-package Rmetapop
### Keywords: package

### ** Examples

~~ simple examples of the most important functions ~~



cleanEx()
nameEx("fishery_simulate")
### * fishery_simulate

flush(stderr()); flush(stdout())

### Name: fishery_simulate
### Title: A wrapper to simulate a linear stage structued model with a
###   stock assessment and harvest
### Aliases: fishery_simulate

### ** Examples

### simulate a population with no fishery 

source("http://mc-stan.org/rstan/stan.R")
### stock assessment parameters for BC Central Coast Herring (Table E7- 2014/15 assessment)
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
n_loc <- 50
C <- 1000

### straying and stochastic info 
stray_s_scale <- c(0.1,10)  ### Cauchy straying scale
rec_corr_sd <-c(0.1)       ### spatial standard deviation in recruitment correlations
surv_rho  <- 0.05       ### scale of stochasticity in survival
rec_sd  <- 0.6          ### lognormal recruitment sd
phi  <- 0.5             ### AR(1) autoregressive coefficient in recruitment 
obs_sd <- 0.3           ### Lognormal survey error 

### fishery info 
Fmort <- 0              ### Maximum target fishing mortality (none in this case)
h_crit <- 0.25          ### Lower biomass threshold

### warmup and number of simulations.  Warmup = ###
warmup1 <- 21           ### Warmup prior to inducing fishing 
warmup2 <-21            ### Warmup prior to measuring statistics 
totiter <- 22           ### Total number of iterations 

### function to convert rows of a dataframe to a list for parallel processing
rows.to.list <- function( df ) {
  ll<-apply(df,1,list)
  ll<-lapply(ll,unlist)
}

stock_IDs= rep(1:5,each= 10)

### dataframe of all parameter combinations for the simulations ###
param.df <- expand.grid(stray_s_scale= stray_s_scale,rec_corr_sd=rec_corr_sd,
                        surv_rho= surv_rho, rec_sd=rec_sd,phi=phi,Fmort= Fmort)

### list of parameters for parallel runs
param.list <- rows.to.list(param.df)  

### complile the assessment model so it need not be directly compiled each time ###
compiled <-fishery_simulate(n_loc=n_loc,n_stages=n_stages,stage_mat=stage_mat,
                         a_bh=a_bh, b_bh=b_bh, phi=0.5,M=M, spat_scale=0.5,
                         spat_sd=0.0001,surv_rho = 0.05,
                         site_sd=0.2,C=1000,obs_sd=0.3,stock_IDs= rep(1:5,each= 10), 
                         n_iter=22,warmup= 21,point.estimate=FALSE,
                         h_crit=0.2, Fmort=0.2)
             
parallel.sim <- function(x){  
    x <- as.data.frame(param.list[[1]])
    
    assess <- fishery_simulate(n_loc=n_loc,                  ### total number of subpopulations
                               n_stages=n_stages,            ### number of age/stage classes
                               stage_mat=stage_mat,          ### stage/age at maturity
                               a_bh=a_bh,                    ### BH alpha parameter
                               b_bh=b_bh,                    ### BH beta parameter
                               phi=x$phi,                    ### lag 1 recruitment autocorrelation 
                               M=M,                          ### natural mortality (a constant, vector or matrix)
                               spat_scale=x$stray_s_scale,   ### spatial Cauchy straying scale
                               spat_sd=x$rec_corr_sd,        ### spatial SD of recruitment synchrony
                               surv_rho=x$surv_rho,          ### degree of stochasticity in survival 
                               site_sd=x$rec_sd,             ### degree of stochasticity in alpha
                               C=C,                       ### degree of stochasticity in straying (lower = more stochastic)
                               obs_sd=obs_sd,                   ### observation error on surveys of total biomass
                               stock_IDs= rep(1:5,each= 10), ### ids representing to which stock the populations belong
                               n_iter=totiter,               ### total number of iterations
                               warmup= warmup1,              ### warmup prior to harvest
                               point.estimate=FALSE,         ### assessment point estimate vs full posterior
                               h_crit=h_crit,                  ### lower harvest cutoff
                               Fmort=x$Fmort                 ### proportional fishing rate
                               )
    
    return(list(B = assess$B, B_stocks = assess$B_stocks,Eggs= assess$ages[1,,],
                assessment = assess$assess, Harvest= assess$Harvest, 
                cors= assess$corrs, B0= assess$B0,
                ages = assess$ages))
}

### run the simulation for two examples
ptm <- proc.time()
fit1 <- mclapply(param.list,parallel.sim,mc.cores= 2) 
proc.time() - ptm

### fuction to plot the simulation
plot_assessment <- function(fit,warmup1=71,warmup2=72,totiter=72,stock=3){
  rgb2 <- tim.colors(10)
  rgb.palette <- colorRampPalette(rgb2, space = "rgb",bias= 2)
  B0.df <- fit$B0
  pops <- melt(fit$B)
  names(pops) <- c("population","time","Biomass")
  
  pops <- join(pops,data.frame(list(population=unique(pops$population),location=paste("Stock",rep(1:(dim(fit$B_stocks)[1]),each = dim(fit$B)[1]/dim(fit$B_stocks)[1]),sep= "_"))))
  
  B0.df$location <- paste("Stock",B0.df$location,sep= "_")
  pred <- melt(apply(fit$assessment[,-1,,c(2)],1,diag),
               varnames=c("time","site"),value.name= "pred")
  
  CSUL <- melt(apply(fit$assessment[,-1,,c(5)],1,diag),
               varnames=c("time","site"),value.name= "CSUL")
  
  CSLL <- melt(apply(fit$assessment[,-1,,c(6)],1,diag),
               varnames=c("time","site"),value.name= "CSLL")
  forecast <- join(pred,CSUL)
  forecast <- join(forecast,CSLL)
  names(forecast)[2] <- "location"
  mean <- melt(apply(fit$assessment[,-totiter,,1],1,diag),
               varnames=c("time","site"),value.name= "mean")
  
  Harvest <- melt(fit$Harvest)
  names(Harvest) <- c("time","location","Harvest")
  forecast_array <- join(pred,mean,by =c("time","site"))
  names(forecast_array)[2] <- "location"
  
  surveys <- melt(fit$B_stocks)
  names(surveys) <- c("location","time","Source","Biomass")
  
  assess2 <- melt(fit$assessment)
  names(assess2) <- c("location","time","iteration","Variable","Biomass")
  levels(surveys$Source) <- c("Aggregate Stock Biomass","Survey Index")
  assess3 <- reshape(assess2, timevar = "Variable", idvar = c("iteration","time","location"), direction = "wide")
  names(assess3)[4] <- "ESS"
  
  assess4 <- subset(assess3,iteration ==72)
  assess4$Posterior <- factor("95% Credible Interval \n(HPD from assessment posterior ) ")
  assess4$PosteriorMean <- factor("Assessment Posterior Mean")
  levels(pops$population) <-gsub("Location_","Population_",levels(pops$population))
  plot.options <- theme(strip.background = element_rect(fill = NA, colour = NA), 
                        panel.grid.minor = element_line(colour = NA), 
                        panel.grid.major = element_line(colour = NA), 
                        legend.key = element_blank(), 
                        axis.text.y = element_text(angle= 90,hjust=0.5),
                        legend.key.size = unit(0.6, "cm"), 
                        legend.direction = "horizontal", 
                        legend.position= "top",
                        legend.box= "horizontal",
                        panel.background = element_rect(fill = NA, colour = "black"),
                        legend.key.width = unit(1.5, "lines"), 
                        panel.margin = unit(1,  "lines"),
                        legend.background = element_blank(),
                        plot.background = element_rect(fill = NA,colour= "black"), 
                        plot.margin = unit(c(0.075, 0.075, 0.075, 0.2), "inches"))
  
  plot1 <-  ggplot(data= subset(assess4,location==paste0("Stock_",stock)))+
    geom_rect(xmin = warmup2,xmax= totiter,ymin= 0,ymax= subset(B0.df,location==paste0("Stock_",stock))$B0,
              fill= alpha("grey90",0.1),colour= "black")+
    geom_ribbon(aes(x=time,ymax= Biomass.ESS_UP,ymin= Biomass.ESS_LOW,
                    fill= Posterior))+  
    geom_pointrange(aes(y= pred,ymin= CSLL,ymax= CSUL,x= time+1),data=subset(forecast,location==paste0("Stock_",stock)),col= "grey45",size= 1)+
    geom_line(aes(y= ESS,x= time,linetype= PosteriorMean,colour=PosteriorMean))+       
    geom_point(aes(y= ESS,x= time,size= PosteriorMean,colour=PosteriorMean,size = PosteriorMean))+       
    geom_point(aes(time,Biomass,colour= Source,size= Source),
               data= subset(surveys,location==paste0("Stock_",stock)))+
    geom_line(aes(time,Biomass,colour= Source,size= Source,linetype= Source),
              data=  subset(surveys,location==paste0("Stock_",stock)))+
    # facet_wrap(~location,ncol= 2)+
    theme_bw()+
    coord_cartesian(xlim= c(0,totiter+2))+
    ylab("biomass")+
    scale_linetype_manual(values=c(1,1,0),name= "")+
    scale_colour_manual(values= c("blue","black","black"),name= "")+
    scale_size_manual(values= c(2,0,2),name= "")+
    scale_shape_manual(values= c(18,18,18),name= "")+
    scale_y_continuous(breaks= c(0,0.25*subset(B0.df,location==paste0("Stock_",stock))$B0,
                                 subset(B0.df,location==paste0("Stock_",stock))$B0),
                       labels= c(0,expression(paste(0.25,B[0])),expression(B[0])))+
    geom_hline(aes(yintercept= B0),
               data= subset(B0.df,location==paste0("Stock_",stock)),linetype="dotted")+
    geom_hline(aes(yintercept= 0.25*B0),colour= "red",data=  subset(B0.df,location==paste0("Stock_",stock)))+
    scale_fill_manual(values= c("grey70"),name= "")+
    geom_vline(xintercept = warmup1)+
    theme_bw()+ plot.options+theme(legend.position= "top",legend.box= "horizontal")
  
  plot2 <-  ggplot(data= subset(surveys,location==paste0("Stock_",stock)&Source=="Aggregate Stock Biomass"))+
    geom_rect(xmin = warmup2,xmax= totiter,ymin= 0,ymax= subset(B0.df,location==paste0("Stock_",stock))$B0/(dim(fit$B)[1]/dim(fit$B_stocks)[1]),
              fill= alpha("grey90",0.1),colour= "black")+
    #facet_wrap(~location,ncol= 2)+
    theme_bw()+
    geom_line(aes(time,Biomass,colour= population),data= subset(pops,location==paste0("Stock_",stock)))+
    geom_line(aes(time,Biomass/(dim(fit$B)[1]/dim(fit$B_stocks)[1])),colour= "blue",size=2)+
    scale_colour_manual(values= rgb.palette(10),guide= guide_legend(ncol= 5),name= "")+
    coord_cartesian(xlim=c(0,totiter+2))+
    scale_y_continuous(breaks= c(0,0.25*subset(B0.df,location==paste0("Stock_",stock))$B0/(dim(fit$B)[1]/dim(fit$B_stocks)[1]),subset(B0.df,location==paste0("Stock_",stock))$B0/(dim(fit$B)[1]/dim(fit$B_stocks)[1])),labels= c(0,expression(paste(0.25,B[0])),expression(B[0])))+
    ylab("scaled biomass")+
    geom_hline(aes(yintercept= B0/(dim(fit$B)[1]/dim(fit$B_stocks)[1])),data=  subset(B0.df,location==paste0("Stock_",stock)),linetype="dotted")+
    geom_hline(aes(yintercept= 0.25*B0/(dim(fit$B)[1]/dim(fit$B_stocks)[1])),colour= "red",data=  subset(B0.df,location==paste0("Stock_",stock)))+
    geom_vline(xintercept = warmup1)+
    theme_bw()+ plot.options+theme(legend.position= "top")
  frame()
  vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
  pushViewport(viewport(layout=grid.layout(100,100)))
  print(plot1, vp=vplayout(1:49,1:100), more = T)
  print(plot2, vp=vplayout(51:100,1:100), more = T)  
  grid.text("A", y=unit(1, "npc") - unit(.05, "npc"),x=unit(1, "npc") - unit(.95, "npc"),hjust=0)
  grid.text("B", y=unit(1, "npc") - unit(.55, "npc"),x=unit(1, "npc") - unit(.95, "npc"),hjust=0)
}

### low straying
plot_assessment(fit1[[1]])
dev.off()
### high straying
plot_assessment(fit1[[2]])
dev.off()
par(mfrow= c(1,2))
pacf(I(fit1[[1]]$assessment[3,,72,1]-fit1[[1]]$B_stocks[3,,1])[1:72], main = "Low Straying")
pacf(I(fit1[[2]]$assessment[3,,72,1]-fit1[[1]]$B_stocks[3,,1])[1:72], main = "High Straying")







graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ran_stray_prob")
### * ran_stray_prob

flush(stderr()); flush(stdout())

### Name: ran_stray_prob
### Title: Generate a series of random stray matrices using a mean stray
###   probability matrix
### Aliases: ran_stray_prob

### ** Examples

stray_mat <- spat_cor_mat(n_loc=10,spat_sd=2, sumto1=TRUE)
stray_mat

rand_strays1 <- ran_stray_prob(stray_mat,n_iter=3,scale=100)
rand_strays2 <- ran_stray_prob(stray_mat,n_iter=3,scale=1000)
rand_strays3 <- ran_stray_prob(stray_mat,n_iter=3,scale=10^5)

df1 <- melt(rand_strays1,varnames = c("site1","site2","replicate"),value.name="stray_probability")
df2 <- melt(rand_strays2,varnames = c("site1","site2","replicate"),value.name="stray_probability")
df3 <- melt(rand_strays3,varnames = c("site1","site2","replicate"),value.name="stray_probability")

df4 <- rbind(df1,df2,df3)
df4$precision <- factor(rep(c(1,2,3),each=nrow(df4)/3), labels= c("Precision = 10^2", "Precision =10^3", "Precision =10^5"))
df4$site1 <- factor(df4$site1)

ggplot(aes(site2,stray_probability),data =df4)+
  geom_line(aes(colour=site1))+
  labs(colour="source site")+ylab("straying probability (random Dirichlet samples")+
  scale_x_continuous(breaks= c(1:10))+
  xlab("destination site")+
  facet_grid(precision~replicate)+theme_bw()



cleanEx()
nameEx("ran_surv_prob")
### * ran_surv_prob

flush(stderr()); flush(stdout())

### Name: ran_surv_prob
### Title: Generate a series of correlated, random survival matrices
### Aliases: ran_surv_prob

### ** Examples

stray_mat <- spat_cor_mat(n_loc=10,spat_sd=2, sumto1=TRUE)
stray_mat

rand_strays1 <- ran_stray_prob(stray_mat,n_iter=3,scale=100)
rand_strays2 <- ran_stray_prob(stray_mat,n_iter=3,scale=1000)
rand_strays3 <- ran_stray_prob(stray_mat,n_iter=3,scale=10^5)

df1 <- melt(rand_strays1,varnames = c("site1","site2","replicate"),value.name="stray_probability")
df2 <- melt(rand_strays2,varnames = c("site1","site2","replicate"),value.name="stray_probability")
df3 <- melt(rand_strays3,varnames = c("site1","site2","replicate"),value.name="stray_probability")

df4 <- rbind(df1,df2,df3)
df4$precision <- factor(rep(c(1,2,3),each=nrow(df4)/3), labels= c("Precision = 10^2", "Precision =10^3", "Precision =10^5"))
df4$site1 <- factor(df4$site1)

ggplot(aes(site2,stray_probability),data =df4)+
  geom_line(aes(colour=site1))+
  labs(colour="source site")+ylab("straying probability (random Dirichlet samples")+
  scale_x_continuous(breaks= c(1:10))+
  xlab("destination site")+
  facet_grid(precision~replicate)+theme_bw()



cleanEx()
nameEx("spat_cor_mat")
### * spat_cor_mat

flush(stderr()); flush(stdout())

### Name: spat_cor_mat
### Title: Generate a spatially correlated matrix
### Aliases: spat_cor_mat

### ** Examples

#a correlation matrix with 10 sites and sd of 4
round(spat_cor_mat(n_loc=10,spat_sd=4,sumto1=FALSE),3)

#a stray th 10 sites and sd of 4
round(spat_cor_mat(n_loc=10,spat_sd=4,sumto1=TRUE),3)

# illustrate different parameter values
library(RColorBrewer)

a <- data.frame(n_loc= rep(10,3),spat_scale=c(1,4,10),spat_sd=c(1,4,10))

b <- array(cbind(with(a,mapply(spat_cor_mat,spat_scale=spat_scale,n_loc=n_loc)),
           with(a,mapply(spat_cor_mat,spat_sd=spat_sd,n_loc=n_loc))),dim= c(10,10,6))

df1 <- melt(b,varnames = c("site1","site2","parameter"),value.name="correlation")
df1$site1 <- factor(df1$site1)
df1$parameter <- factor(df1$parameter, labels= 
  c(paste(rep(c("Cauchy, scale=","Gaussian, sd="),each= 3),rep(c(1,4,10),2))))

ggplot(aes(site2,correlation),data =df1)+
  geom_line(aes(colour=site1))+
  labs(colour="site")+
  scale_x_continuous(breaks= c(1:10))+
  xlab("site")+
  facet_wrap(~parameter)+theme_bw()



cleanEx()
nameEx("spat_temp_ts")
### * spat_temp_ts

flush(stderr()); flush(stdout())

### Name: spat_temp_ts
### Title: Generate a spatially and temporally correlated time series
### Aliases: spat_temp_ts

### ** Examples


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



cleanEx()
nameEx("ssr_linear")
### * ssr_linear

flush(stderr()); flush(stdout())

### Name: ssr_linear
### Title: Discrete survival, straying, and recruitment projection.
### Aliases: ssr_linear

### ** Examples


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





cleanEx()
nameEx("ssr_linear_ode")
### * ssr_linear_ode

flush(stderr()); flush(stdout())

### Name: ssr_linear_ode
### Title: Continuous time survival, straying, and harvest with discrete
###   recruitment
### Aliases: ssr_linear_ode

### ** Examples


### set up parameters for the simulation test survival & migration function with two
### locations ### general parameters ###
n_loc <- 5  ### number of locations
n_iter <- 100  ### number of years
n_stages <- 10  ### number of stages including eggs and recruits
stage_mat <- 3  ### stage # that indicates maturity 

### 9 x L mortality matrix (can also be a single number or a single row with L
### columns)
mort <- matrix(0.334, ncol = n_loc, nrow = length((stage_mat - 1):(n_stages)))

### create identifier of matrix entries that correspond to mortality for quick
### replacement in the array
mort_mat_id <- array(matrix(1:n_stages^2, ncol = n_stages) %in% 
                       c(diag(matrix(1:n_stages^2, ncol = n_stages)[-c(1), ]), n_stages^2),
                       dim = c(n_stages, n_stages, n_loc))[-1,-1, ]

### fecundity at age 3-10
fec_at_age <- fecundity_age(stage_mat:(n_stages))

### biomass at age 3-10 in tons
tons_at_age <- LVBweight(stage_mat:(n_stages))/1000

### straying probability matrix (must be square with L x L dimensions)
spat_scale <- 0.3

### matrix of mean stray probabilities
stray_probs <- spat_cor_mat(n_loc, spat_scale = spat_scale, sumto1 = TRUE)

### dirichlet scale parameter higher = less stochastic probabiliies
### dirichlet scale parameter higher = less stochastic probabiliies
C <- 1000

### stochastic realizations of the stray matrix for each year using the dirichlet
### distribution
Crand <- ran_stray_prob(stray_mat=stray_probs,n_iter=n_iter,scale= C)

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
site_sd = 0.5  ### overall log-scale recruitment standard deviation 

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
X <- array(as.numeric(NA), dim = c(n_stages, n_loc, n_iter), dimnames = list(d1, d2, d3))
X[, , 1:2] <- matrix(rep(1e+08, each = n_loc), ncol = n_loc, byrow = T)

## set up the adult biomass array
B <- array(NA, dim = c(n_loc, n_iter+1, 4), dimnames = list(d2, 1:(n_iter+1), c("Actual", "Observed","Estimate","Forecast")))
B[, 1, ] <- tons_at_age %*% X[3:10, , 1]/1000
B[, 2, ] <- tons_at_age %*% X[3:10, , 2]/1000

## set up the estimated spawner biomass array
BF <- array(NA, dim = c(n_loc, n_iter, n_iter), dimnames = list(d2, d3, d3))

## set up the estimated forecast spawner biomass array
BF2 <- array(NA, dim = c(n_loc, n_iter+1, n_iter), dimnames = list(d2, 1:(n_iter+1),d3))

## set up the forecast abundance at age array
AF <- array(NA, dim = c(n_stages,n_loc,  n_iter+1, n_iter), dimnames = list(d1,d2,1:(n_iter+1),d3))

## set up the size frequnecy array
S_freq <- array(NA, dim = c(n_stages-2,n_loc, n_iter), dimnames = list(d1[-c(1,2)],d2,d3))

### log-scale observer errorr
obs_sd <- 0.3


### project the population with stochastic recruitment ###

  for (i in 1:2) {
    S_freq[,,i]  <- t(t(X[-c(1,2),,i])/colSums(X[-c(1,2),,i]))
  }
system.time(for (i in 3:n_iter) {
    X[, , i] <- ssr_linear_ode(alpha=alpha,beta=beta, fec_at_age = fec_at_age, 
                             n_loc = n_loc, n_stages = n_stages, stage_mat = stage_mat, 
                             eggs = X[1, , i - 2], 
                             X0 = X[(stage_mat - 1):n_stages, , i - 1], 
                             Z = mort, stray = Crand[, , i - 1],
                             inst_h= 0, 
                             errors = errors$ts[i - 1, ])$X
  
    ### age frequency ###
    S_freq[,,i]  <- t(t(X[-c(1,2),,i])/colSums(X[-c(1,2),,i]))
    
    ### actual biomass
    B[, i, 1] <- tons_at_age %*% X[3:n_stages, , i]/1000
    
    ### observed biomass with observer error
    B[, i, 2] <- B[, i, 1] * exp(rnorm(n_loc, 0, obs_sd) - 0.5 * obs_sd^2)
  }
)
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

dataL = melt(BF2, id="x",varnames= c("Location","Time","Iteration"),value.name= "Tons")

ggplot(aes(Time, Tons), data = df_ts2) + 
  geom_line(aes(Time,Tons,group= Iteration),colour= "blue",data= dataL,size= 0.3)+
  geom_line(aes(colour = Source, linetype = Source),size= 0.5) +
  geom_point(aes(colour = Source,size= Source)) +
  scale_size_manual(values= c(0,1,0,2))+
  scale_linetype_manual(values = c(1, 3, 1,0)) + 
  plot.options + facet_grid(Location ~ .) + 
  scale_colour_manual(values = c("red", "black", "blue","blue")) + 
  ylab("Tons (Thousands)") + 
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))
  

### correlation in recruitment among sites
levelplot(cor(t(X[2, , ])))

### correlation in adult biomass among sites
levelplot(cor(t(B[, 1:n_iter, 1])))








### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
