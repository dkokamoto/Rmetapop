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
warmup1 <- 71           ### Warmup prior to inducing fishing 
warmup2 <-71            ### Warmup prior to measuring statistics 
totiter <- 72           ### Total number of iterations 

### function to convert rows of a dataframe to a list for parallel processing
rows.to.list <- function( df ) {
  ll<-apply(df,1,list)
  ll<-lapply(ll,unlist)
}

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
                         n_iter=3,warmup= 2,point.estimate=FALSE,
                         h_crit=0.2, Fmort=0.2)
             
parallel.sim <- function(x){  
    x <- as.data.frame(t(x))
    
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
  B0.df <- fit[7]$B
  pops <- melt(fit$B)
  names(pops) <- c("population","time","Biomass")
  
  pops <- join(pops,data.frame(list(population=unique(pops$population),location=paste("Stock",rep(1:(dim(fit$B_stocks)[1]),each = dim(fit$B)[1]/dim(fit$B_stocks)[1]),sep= "_"))))
  
  B0.df$location <- paste("Stock",B0.df$location,sep= "_")
  pred <- melt(apply(fit$assessment[,-1,,c(2)],1,diag),
               varnames=c("time","site"),value.name= "pred")
  
  CSUL <- melt(apply(fit$assessment[,-1,,c(3)],1,diag),
               varnames=c("time","site"),value.name= "CSUL")
  CSLL <- melt(apply(fit$assessment[,-1,,c(3)],1,diag),
               varnames=c("time","site"),value.name= "CSLL")
  forecast <- cbind(pred,CSUL,CSLL)
  
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
    geom_line(aes(y= ESS,x= time,linetype= PosteriorMean,colour=PosteriorMean))+       
    geom_point(aes(y= ESS,x= time,size= PosteriorMean,colour=PosteriorMean,size = PosteriorMean))+       
    geom_point(aes(time,Biomass,colour= Source,size= Source),
               data= subset(surveys,location==paste0("Stock_",stock)))+
    geom_line(aes(time,Biomass,colour= Source,size= Source,linetype= Source),
              data=  subset(surveys,location==paste0("Stock_",stock)))+
    # facet_wrap(~location,ncol= 2)+
    theme_bw()+
    coord_cartesian(xlim= c(0,totiter+1))+
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
    coord_cartesian(xlim=c(0,totiter+1))+
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

plot_assessment(fit1[[1]])

pacf(I(fit1[[1]]$assessment[3,,72,1]-fit1[[1]]$B_stocks[3,,1])[1:72])


