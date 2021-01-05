library(tidyverse)
library(lattice)
library(viridisLite)
library(Rfast)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(scales)

#-----------------------------------------------------------
#1. Heat Maps for 3 Guild Variables 
load("guildinfodata_0522.RData")

guildinfodata<-rbind(guildinfodata_0522[1:49,],
                     c(50,"2011-01-29",1,50,rep(0,5)),
                     guildinfodata_0522[50:1549,])

K<-50
guild_size<-matrix(0,K,31)
guild_spend<-matrix(0,K,31)
guild_meminteract<-matrix(0,K,31)
guild_noofpurch<-matrix(0,K,31)
guild_avggame<-matrix(0,K,31)

for(i in 1:K){
  
  temp<-data.frame(guildinfodata[guildinfodata$guildindex==i,
                                 5:9],rownames.force = NA)
  
  temp[is.na(temp)]<-0
  
  guild_size[i,]<-as.numeric(temp[,5])
  guild_spend[i,]<-as.numeric(temp[,3])
  guild_meminteract[i,]<-as.numeric(temp[,1])
  guild_noofpurch[i,]<-as.numeric(temp[,4])
  guild_avggame[i,]<-as.numeric(temp[,2])
}
rm('guildinfodata_0522','temp')

x.scale <- list(at=seq(0,30,5))

coul <- heat.colors(100)
p1<-levelplot(t(guild_meminteract[1:(K-1),1:30]),at=seq(0, 50, length.out=30),
              col.regions=coul,
              scales=list(x=x.scale),xlab="Days",ylab="Guilds",main="No. of guild members who played in team")

p2<-levelplot(t(guild_avggame[1:(K-1),1:30]),at=seq(0, 8, length.out=30),
              col.regions=coul,
              scales=list(x=x.scale),xlab="Days",ylab="Guilds",main="No. of daily average game sessions")

p3<-levelplot(t(guild_noofpurch[1:(K-1),1:30]),at=seq(0, 10, length.out=30),
          col.regions=coul,
          scales=list(x=x.scale),xlab="Days",ylab="Guilds",main="No. of purchases by guild members")

grid.arrange(p1,p2,p3,ncol=3)

#--------------------------------------------------------------

id.remove<-c('412707','626494','632811','981252',
             '987762','8710791','9569599','9600937',
             '10182266','10360113','10502933','10632317',
             '10638416','10786900','10921298','10987066',
             '11509405','18892676','19201223','19222498')

#2.Now with matrix x of player specific variables
load("data_0522.RData")
iv<-data_0522
for(ii in 1:length(id.remove)){
  
  idx<-which(iv$characterid==id.remove[ii])
  iv<-iv[-idx,]
}

#3. Get login ind. and login time
library(readr)
dv_0522 <- read_csv("dv_0522.txt")
dv<-dv_0522
for(ii in 1:length(id.remove)){
  
  idx<-which(dv$characterid==id.remove[ii])
  dv<-dv[-idx,]
}

#4. Do some clean up

#first: logindummy=1 and logintime = 0. Replace logintime with 1
dv$logintime[dv$logindummy==1 & dv$logintime==0]<- 1

#second: logindummy = 0 and pvpplaytime or pvpkillpoint >0. Replace logindummy = 1
# and logintime = max(1,pvpplaytime+pveplaytime)
data.full<-cbind(iv,dv)
data.full[is.na(data.full)]<-0
data.full$logindummy[(data.full$pvpplaytime>0 | data.full$pvpkillpoint>0) & data.full$logindummy==0]<- 1
data.full$logintime[data.full$logindummy==1 & data.full$logintime==0]<- data.full$pvpplaytime[data.full$logindummy==1 & data.full$logintime==0]+
  data.full$pvetime[data.full$logindummy==1 & data.full$logintime==0]
data.full$logintime[data.full$logindummy==1 & data.full$logintime==0]<- 1

#third: level should be 1 for the players above
data.full$level[data.full$logindummy==1 & data.full$level==0]<- 1
data.full$friendmeanlevel[data.full$logindummy==1 & data.full$friendmeanlevel==0]<- 1
rm('data_0522','dv_0522')

#4. Add new variables (time_since, week, holiday)
mydate<- unique(as.Date(data.full$date,format="%Y-%m-%d"))
week<- matrix(0,length(mydate),1)
week<- 1*(weekdays(mydate)=="Saturday" | weekdays(mydate)=="Sunday")
holiday<- matrix(0,length(mydate),1)
holiday[mydate=="2011-02-03" | mydate=="2011-02-14"]<- 1

# Get the Guildmembership matrix
load("data_0522_guildmemship.RData")
guild_mem<- data_0522_guildmemship
rm('data_0522_guildmemship')
for(ii in 1:length(id.remove)){
  
  idx<-which(guild_mem$characterid==id.remove[ii])
  guild_mem<-guild_mem[-idx,]
}
K<-50
gm.all<-guild_mem
rm('guild_mem')
idx<-(gm.all$guilddummy>0)
gm.all<-gm.all[idx,]
rm('idx')

new_dv<-cbind(data.full,gm.all$characterid,gm.all$guildindex)
rm('gm.all')
rm('dv')

m<-31
plotmat_login<-matrix(0,K,m)
plotmat_purch<-matrix(0,K,m)
plotmat_act<-matrix(0,K,m)

for(k in 1:K){
  for(j in 1:m){
  
  temp<-new_dv[new_dv$dayid==j & new_dv$`gm.all$guildindex`==k,]
  plotmat_login[k,j] = 100*sum(temp$logindummy)/(dim(temp)[1])
  plotmat_purch[k,j] = 100*sum(temp$purchasedummy)/sum(temp$logindummy)
  plotmat_act[k,j] = sum(temp$logintime)/(60*sum(temp$logindummy))
  }
}

plotmat_login<- plotmat_login[1:(K-1),1:(m-1)]
rownames(plotmat_login) = paste("guild", seq(K-1), sep="")
colnames(plotmat_login) = paste("day", seq(m-1), sep="")
dat_login = as.data.frame(plotmat_login)
dat_login$guild = rownames(dat_login)
mdat_login = melt(dat_login, id.vars="guild")
mdat_login$day = as.numeric(gsub("day", "", mdat_login$variable))

plotmat_purch[plotmat_purch==0]<- 0.01
plotmat_purch<- plotmat_purch[1:(K-1),1:(m-1)]
rownames(plotmat_purch) = paste("guild", seq(K-1), sep="")
colnames(plotmat_purch) = paste("day", seq(m-1), sep="")
dat_purch = as.data.frame(plotmat_purch)
dat_purch$guild = rownames(dat_purch)
mdat_purch = melt(dat_purch, id.vars="guild")
mdat_purch$day = as.numeric(gsub("day", "", mdat_purch$variable))

plotmat_act<- plotmat_act[1:(K-1),1:(m-1)]
rownames(plotmat_act) = paste("guild", seq(K-1), sep="")
colnames(plotmat_act) = paste("day", seq(m-1), sep="")
dat_act = as.data.frame(plotmat_act)
dat_act$guild = rownames(dat_act)
mdat_act = melt(dat_act, id.vars="guild")
mdat_act$day = as.numeric(gsub("day", "", mdat_act$variable))

g1 = ggplot()+theme_bw()+ylab("% of Players Logged In")+
  xlab("Days")+
  geom_smooth(data=mdat_login, aes(x=day, y=value, group=guild),
              color='violet',size=0.5,span=0.2,se=FALSE)+
  geom_smooth(data=mdat_login, aes(x=day, y=value),
              color='red',size=2, alpha=0.5,span=0.01)+
  scale_x_continuous(breaks = seq(0,30,2))+
  geom_vline(xintercept = which(week==1), 
             color = "black", size=0.1)+
  geom_vline(xintercept = which(holiday==1),linetype="dotted", 
             color = "blue", size=0.1)
  

g2 = ggplot() + theme_bw() + 
  ylab("% of Players with Purchase >0")+xlab("Days")+
  geom_smooth(data=mdat_purch, aes(x=day, y=value, group=guild),
              color='orange',size=0.5,alpha = 0.2,span=0.2,se=FALSE)+
  geom_smooth(data=mdat_purch, aes(x=day, y=value),
              color='red',size=2, alpha=0.5,span=0.01)+
  scale_x_continuous(breaks = seq(0,30,2))+
  scale_y_continuous(limit=c(0,NA),oob=squish)+
  geom_vline(xintercept = which(week==1), 
             color = "black", size=0.1)+
  geom_vline(xintercept = which(holiday==1),linetype="dotted", 
             color = "blue", size=0.1)

g3 = ggplot() + theme_bw() + scale_x_continuous(breaks = seq(0,30,2))+
  ylab("Duration (minutes)")+xlab("Days")+
  geom_smooth(data=mdat_act, aes(x=day, y=value, group=guild),
            color='green',size=0.5,alpha = 0.2,span=0.2,se=FALSE)+
  geom_smooth(data=mdat_act, aes(x=day, y=value),
              color='red',size=2, alpha=0.5,span=0.01)+
  geom_vline(xintercept = which(week==1), 
             color = "black", size=0.1)+
  geom_vline(xintercept = which(holiday==1),linetype="dotted", 
             color = "blue", size=0.1)

grid.arrange(g1,g3,g2,ncol=3)


#------------------------- P matrix ---------------------------------

Pmat<-function(TT,t.val){
  s<-matrix(0,TT,TT)
  for(i in 1:TT){
    s[i,]<-1*(abs(i-c(1:TT))<=t.val)
  }
  e1<-cbind(s,s,s)
  e<-rbind(e1,e1,e1)
  P <- matrix(1, 3*TT, 3*TT)-e
  return(P)
}
TT<- 5
P1<- Pmat(TT,0)
P2<- Pmat(TT,1)
P3<- Pmat(TT,2)

x.scale <- list(at=seq(1,3*TT,1))
out1<-levelplot(P1,col.regions = c(gray(0.75),gray(0.25)),
                scales=list(x=x.scale,y=x.scale),
                xlab="t' = 0",ylab="",colorkey=FALSE)
out2<-levelplot(P2,col.regions = c(gray(0.75),gray(0.25)),
                scales=list(x=x.scale,y=x.scale),
                xlab="t' = 1",ylab="",colorkey=FALSE)
out3<-levelplot(P3,col.regions = c(gray(0.75),gray(0.25)),
                scales=list(x=x.scale,y=x.scale),
                xlab="t' = 2",ylab="",colorkey=FALSE)
grid.arrange(out1,out2,out3,ncol=3)
