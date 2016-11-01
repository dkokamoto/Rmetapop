
plist <- c("Rmetapop","deSolve","ggplot2","grid","gtools","mvtnorm","plyr",
               "reshape2","slam","scales","coda","rstan","VGAM","MASS",
           "grid","scales","parallel","gridExtra","dplyr")
sapply(plist, FUN = function(X) {
  do.call("require", list(X)) 
})


### stock assessment parameters (Table E7- 2014/15 assessment)
E0_assess <-  59463*1e6*200  ### eggs at equilibrium
R0_assess <- 500*1e6  ### modified from table E7`
h_assess <- 0.8  ### steepness
M_assess <- 0.5  ### natural mortality
stage_mat <- 3  ### stage at maturity
a_bh_assess <- E0_assess  * (1 - h_assess)/(4 * h_assess  * R0_assess )
b_bh_assess <- (5 * h_assess - 1)/(4 * h_assess  * R0_assess)

#### additional key parameters
C <- 1000000  ### little stochasticity in transitions
n_stages= 10
obs_sd = 0.4
h_crit= 0.25
### other unique parameters
n_stocks <- 2
n_loc = n_stocks*10

### set cauchy scale to achieve desired retention 
test <- rev(seq(from= 0.11,to = 0.99, length.out= 21))
s_to_ret <- function(x) (mean(diag(spat_cor_mat(n_loc*n_stocks,spat_scale=x,sumto1=T))))
stray_s_scale <- optim(stray_s_scale,function(p) sum((sapply(p,s_to_ret)-test)^2),
                       lower= rep(0.000001,length(stray_s_scale)),method="L-BFGS-B")$par                                               

### generate parameter list
stray_s_scale <-  c(0.045,0.144,0.203,0.25,0.302,0.345,0.398,0.45,0.507,0.57,0.643,
                    0.724,0.811,0.958,1.101,1.407,1.609,2.249,5.85,10.927,15.379)

### range of parameters for the simulation
rec_corr_sd <-  exp(seq(from = log(0.4),to = log(15), length.out= 21))
surv_rho  <- c(0.01)
rec_sd  <- c(0.8)
phi  <- c(0.5)
Fmort <- seq(0,0.5,by = 0.1)
dd <- c(1,0)
spat_alloc <- c(1,3)

rows.to.list <- function( df ) {
  ll<-apply(df,1,list)
  ll<-lapply(ll,unlist)
}

dep_length <- function(x) {
  a <- rle(as.numeric(x))
  ifelse(length(a$lengths[a$values==1])==0,NA,mean(a$lengths[a$values==1]))
}

param.df <- expand.grid(stray_s_scale= stray_s_scale,rec_corr_sd=rec_corr_sd,
                        surv_rho= surv_rho,phi=phi,Fmort= Fmort,
                        dd= dd,spat_alloc=spat_alloc,M_assess=M_assess)

param.df <- param.df[rep(1:nrow(param.df),length(a_bh_assess)),]
param.df$a_bh_assess <- rep(a_bh_assess,each= nrow(param.df)/length(a_bh_assess))
param.df$b_bh_assess <- rep(b_bh_assess,each= nrow(param.df)/length(a_bh_assess))


### adjust recruitment SD
coefs <- c(0.934164949,-0.306527755,-0.001890675,0.059833737,-0.020983389,0.002484043)
param.df$rec_sd <- cbind(1,poly(log(param.df$rec_corr_sd),5,raw= TRUE))%*%t(t(coefs))

### convert df of parameters to list for parallel processing
param.list <- rows.to.list(param.df)    

### compile the model 
fit.compile <- nlss_compile()

num = 950
### extract parameter of interest
  x <- param.df[num,]
  
  system.time(rep <- fishery_simulate(n_loc=n_loc,      ### total number of subpopulations
                                      n_stages=n_stages,            ### number of age/stage classes
                                      stage_mat=stage_mat,          ### stage/age at maturity
                                      a_bh=x$a_bh_assess,                    ### BH alpha parameter
                                      b_bh=x$b_bh_assess,                    ### BH beta parameter
                                      phi=x$phi,                    ### lag 1 recruitment autocorrelation 
                                      M=x$M_assess,                   ### natural mortality (a constant, vector or matrix)  
                                      assessment= TRUE,            ### whether to assess the stock
                                      spat_alloc= x$spat_alloc,
                                      collective_dd=x$dd,           ### group density dependence
                                      spat_scale=x$stray_s_scale,   ### spatial Cauchy straying scale
                                      spat_sd=x$rec_corr_sd,        ### spatial SD of recruitment synchrony
                                      site_sd=x$rec_sd,             ### degree of stochasticity in alpha
                                      C=C,                          ### inverse of degree of stochasticity in straying 
                                      obs_sd=obs_sd,                ### observation error on surveys of total biomass
                                      stock_IDs= rep(1:n_stocks,each= n_loc/n_stocks), ### ids representing to which stock the populations belong
                                      n_iter=53,                    ### total number of iterations
                                      warmup= 8,  
                                      warmup2= 13,                    ### warmup prior to harvest
                                      point.estimate=FALSE,         ### assessment point estimate vs full posterior
                                      h_crit=h_crit,                ### lower harvest cutoff
                                      Fmort=x$Fmort,                ### proportional fishing rate
                                      fit.compile=fit.compile,      ### compiled assessment model
                                      ret_ts= TRUE))
  
#system.time(fit_MM <- mclapply(param.list,parallel.sim,mc.cores= 8))
#save(fit_MM,param.df,file= "fit_4.27.16.Rdata")

Assess<- melt(rep$TS$forecast[,,,1:2],value.name= "Estimate",varnames= c("Stock","Year","Assess","type"))%>%
        filter((Year==Assess+1&type=="Prediction")|type=="Estimated Stock Size")
Harvest <- rep$TS$Harvest%>%
  melt(varnames= c("Year","Stock"), value.name= "Harvest")

Biomass<- rep$TS$B_stocks[,,]%>%
  melt(value.name= "Biomass",varnames=c("Stock","Year","Type"))
Loc_Biomass <- melt(rep$TS$B[1:10,],value.name= "Biomass",varnames= c("Loc","Year"))
Cred_Intvl <- data.frame(rep$TS$assess[1,,53,c(3,4)])%>%
  mutate(Year=1:nrow(.))

options <- theme(
    strip.text = element_text(size= 12),
    panel.grid.minor = element_line(colour = NA),
    panel.grid.major = element_line(colour = NA),
    strip.background = element_rect(fill= NA,colour= NA),
    legend.key = element_blank(),
    legend.key.height = unit(0.6, "cm"),
    legend.background = element_rect(fill = NA),
    legend.position= "top",
    legend.direction= "horizontal",
    legend.title= element_blank(),
    axis.title.x= element_blank(),
    legend.box='horizontal',
    legend.text = element_text(size = 12,face= "plain"),
    legend.key.width = unit(1,"lines"),
    axis.ticks.length = unit(.25,"lines"),
    axis.text.y= element_text(size=12,colour= "black",angle=90,hjust=0.5),
    axis.text.x= element_text(size=12,colour= "black"),
    axis.title = element_blank(),
    axis.ticks.margin = unit(.25,"lines"),
    panel.background=element_rect(fill= NA,colour=NA),
    plot.background=element_rect(fill= NA,colour=NA),
    panel.border = element_rect(colour="black",fill=NA),
    panel.margin= unit(1,"lines"),
    plot.margin = unit(rep(.2,4), "inches"))
  
p1 <- ggplot(data= subset(Assess,type== "Estimated Stock Size"&Stock=="Stock_1"))+
  geom_ribbon(aes(ymin= ESS_UP,ymax= ESS_LOW,x= Year),data= Cred_Intvl,fill= "grey80")+
  geom_line(aes(Year,Estimate,group= factor(Assess):type,linetype= type), alpha= 0.5)+
  geom_point(aes(Year,Biomass, fill= Type),data= subset(Biomass,Type!= "actual"&Stock=="Stock_1"))+
  geom_line(aes(Year,Biomass,size= Type),colour= "blue",
            data= subset(Biomass,Type== "actual"&Stock=="Stock_1"))+
  scale_size_manual(values= 2, name= "", labels= c("Actual Biomass"))+
  scale_linetype_manual(values= 1, name= "", labels= c("Assessment Estimate"))+
  scale_fill_manual(values= "grey", name= "", labels= c("Surveyed Biomass"))+
  ylab("scaled biomass")+
  geom_vline(xintercept=8,linetype= "dotted")+
  geom_hline(yintercept= rep$TS$B0[1])+
  geom_hline(yintercept= rep$TS$B0[1]*0.25,col= "red")+
  theme_bw()+
  options+
  coord_cartesian(ylim= c(0,1.5*rep$TS$B0[1]),expand= c(0,0))+
  scale_y_continuous(breaks= c(0,rep$TS$B0[1]*0.25,rep$TS$B0[1]),
                     labels= c(0,expression(paste(0.25,B[0])),expression(paste(B[0]))))

p2 <- ggplot(aes(Year,Biomass,colour= Loc),data= Loc_Biomass)+
  geom_hline(yintercept= rep$TS$B0[1])+
  geom_path(aes(Year,Biomass/10,size= Type),colour= "blue",
            data= subset(Biomass,Type== "actual"&Stock=="Stock_1"))+
  geom_path(size=0.5,alpha=0.8)+geom_hline(yintercept= rep$TS$B0[1]*0.25,col= "red")+
  ylab("")+
  theme_bw()+
  options+
  geom_vline(xintercept=8,linetype= "dotted")+
  geom_hline(yintercept= rep$TS$B0[1]/10)+
  scale_colour_manual(labels=paste("Stocklet",1:10),
                      values= tim.colors(10))+
  geom_hline(yintercept= rep$TS$B0[1]*0.25/10,col= "red")+
  coord_cartesian(ylim= c(0,1.5*rep$TS$B0[1]/10),expand= c(0,0))+
  scale_y_continuous(breaks= c(0,rep$TS$B0[1]*0.25/10,rep$TS$B0[1]/10),
                     labels= c(0,expression(paste(0.25,B[0])),expression(paste(B[0]))))+
  scale_size_manual(values= 2, name= "", labels= c("Stock Mean"))


pdf(width =8, height = 7, file = paste0("/home/okamoto/SIMS_EX_",num,".pdf"), family = "serif",pointsize = 9)
#quartz(width =8, height = 7)
grid.arrange(p1,p2,left = textGrob("Adult Biomass",rot= 90, gp=gpar(fontsize=12,font=12)))
grid.text("A) Stock Scale", y=unit(1, "npc") - unit(.035, "npc"),x=unit(1, "npc") - unit(.92, "npc"),hjust=0, gp=gpar(fontsize=12))
grid.text("B) Stocklet Scale", y=unit(1, "npc") - unit(.52, "npc"),x=unit(1, "npc") - unit(.92, "npc"),hjust=0, gp=gpar(fontsize=12))
dev.off()

