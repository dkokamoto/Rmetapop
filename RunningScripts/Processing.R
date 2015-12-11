load("~/Data/combined_results_2015-05-28_1.RData")

stray_s_scale <- c(0.01,0.05,0.1,0.5,1,2,5,10)
rec_corr_sd <- c(0.01,0.05,0.1,0.5,1,2,5,10)
surv_rho  <- c(0.05)
rec_sd  <- c(0.6)
phi  <- c(0,0.5)
Fmort <- c(0,0.2,0.4)

param.df <- expand.grid(stray_s_scale= stray_s_scale,
                        rec_corr_sd=rec_corr_sd,surv_rho= surv_rho, 
                        rec_sd=rec_sd,phi=phi,Fmort= Fmort)

closures <- function(x) mean(colMeans(x$Harvest[32:72,]==0))
pop_biomass_cv <- function(x) mean(apply(x$B[,32:72],2,sd)/apply(x$B[,32:72],2,mean))
stock_biomass_cv <- function(x) mean(apply(x$B_stocks[,32:72,1],2,sd)/apply(x$B_stocks[,32:72,1],2,mean))

stock_biomass_mu <- function(x) mean(apply(x$B_stocks[,32:72,1],2,mean))
pop_biomass_mu <- function(x) mean(apply(x$B[,32:72],2,mean))

stock_overfished <- function(x) mean(apply((x$B_stocks[,32:72,1])>(x$B0[,3]),2,mean))
pop_overfished <- function(x) mean(apply(x$B[,32:72]>(rep(x$B0[,3]/10,each= 10)),2,mean))

param.df$p_cv <- sapply(fit1,pop_biomass_cv)
param.df$s_cv <- sapply(fit1,stock_biomass_cv)
param.df$s_mu <- sapply(fit1,stock_biomass_mu)
param.df$p_mu <- sapply(fit1,pop_biomass_mu)
param.df$p_of <- sapply(fit1,pop_overfished)
param.df$s_of <- sapply(fit1,stock_overfished)

ggplot(aes(x=stray_s_scale,y=p_cv),data= param.df)+
  geom_line(aes(y= s_cv),colour= "red")+
  geom_line()+
  facet_grid(rec_corr_sd~Fmort)+
  theme_bw()

ggplot(aes(x=stray_s_scale,y=p_of),data= param.df)+
  geom_line(aes(y= s_of),colour= "red")+
  geom_line()+
  facet_grid(rec_corr_sd~Fmort,scales= "free_y")+
  theme_bw()

