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

