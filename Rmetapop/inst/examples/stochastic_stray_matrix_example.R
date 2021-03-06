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
