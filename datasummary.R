

# These ids are being removed because they have NA values, as well as other anomalies.
id.remove<-c('412707','626494','632811','981252',
             '987762','8710791','9569599','9600937',
             '10182266','10360113','10502933','10632317',
             '10638416','10786900','10921298','10987066',
             '11509405','18892676','19201223','19222498')

m<- 31

#1. Get matrix g of guild variables
load("guildinfodata_0522.RData")

guildinfodata<-rbind(guildinfodata_0522[1:49,],
                     c(50,"2011-01-29",1,50,rep(0,5)),
                     guildinfodata_0522[50:1549,])
K<-50
g.full<-matrix(0,K*m,5)

for(i in 1:K){
  
  temp<-data.frame(guildinfodata[guildinfodata$guildindex==i,
                                  5:9],rownames.force = NA)
  s<-1+m*(i-1)
  e<-m*i
  g.full[s:e,]<-cbind(as.numeric(temp[,1]),as.numeric(temp[,2]),as.numeric(temp[,3]),
                  as.numeric(temp[,4]),as.numeric(temp[,5]))
}
g.full[is.na(g.full)]<-0
rm('guildinfodata_0522','temp')

summary_guild<- matrix(0,5,7)

for(i in 1:5){
  
  summary_guild[i,]<- c(100*sum(1*(g.full[,i]==0))/1550,mean(g.full[,i]),
                        quantile(g.full[,i],c(0.25,0.5,0.75,0.95)),
                        sd(g.full[,i]))
}
row.names(summary_guild)<- names(guildinfodata)[c(5:9)]
write.csv(summary_guild,file='summary_guild.csv')


#2. Now with matrix x of player specific variables
load("data_0522.RData")
iv<-data_0522
for(ii in 1:length(id.remove)){
  
  idx<-which(iv$characterid==id.remove[ii])
  iv<-iv[-idx,]
}
unique.id.x<-unique(iv$characterid)
nx<-length(unique.id.x)

#3. Get login ind. and login time
library(readr)
dv_0522 <- read_csv("dv_0522.txt")
dv<-dv_0522
for(ii in 1:length(id.remove)){
  
  idx<-which(dv$characterid==id.remove[ii])
  dv<-dv[-idx,]
}
unique.id.y<-unique(dv$characterid)
ny<-length(unique.id.y)

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

#fourth. Add new variables (time_since, week, holiday)
unique.id<- unique(data.full$characterid)
n<-length(unique.id)
mydate<- as.Date(data.full$date,format="%Y-%m-%d")
week<- weekdays(mydate)
data.full$week<- 1*(week=="Saturday" | week=="Sunday")
data.full$holiday<- 0
data.full$holiday[mydate=="2011-02-03" | mydate=="2011-02-14"]<- 1
rm('mydate','week')

data.full$timesince<- 0

for(i in 1:n){
  
  ii<-unique.id[i]
  temp1<- data.full[data.full$characterid==ii,c('logindummy','dayid')]
  temp2<- matrix(0,m)
  for(j in 2:(m)){
    
    if(temp1$logindummy[j-1]>0){
      temp2[j]<- 0
    }
    if(temp1$logindummy[j-1]==0){
      temp2[j]<- 1+temp2[j-1]
    }
  }
  data.full$timesince[data.full$characterid==ii]<- temp2
}
rm('temp1','temp2')

x.full<-matrix(0,nx*m,21)

for(i in 1:nx){
  ii<-unique.id.x[i]
  temp<-data.matrix(data.full[data.full$characterid==ii,
                       c(4,5,7,9:12,15:21,23:25,32,38:40)],rownames.force = NA)
  s<-1+m*(i-1)
  e<-m*i
  x.full[s:e,]<-temp
}
x.full[is.na(x.full)]<-0

rm('data_0522','temp','dv_0522')

#5. Now get the summary stats for active players

idx<- which(x.full[,18]>0)
x.full.idx<-x.full[idx,]
x.full.idx<-x.full.idx[,-18]
nn<-dim(x.full.idx)[1]

summary_x<- matrix(0,20,7)

for(i in 1:20){
  
  summary_x[i,]<- c(100*sum(1*(x.full.idx[,i]==0))/nn,mean(x.full.idx[,i]),
                        quantile(x.full.idx[,i],c(0.25,0.5,0.75,0.95)),
                        sd(x.full.idx[,i]))
}

summary_x[20,]<- c(100*sum(1*(data.full$timesince==0))/dim(data.full)[1],
                   mean(data.full$timesince),
                   quantile(data.full$timesince,c(0.25,0.5,0.75,0.95)),
                   sd(data.full$timesince))
row.names(summary_x)<- names(data.full)[c(4,5,7,9:12,15:21,23:25,38:40)]
write.csv(summary_x,file='summary_x.csv')
