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
