

m1<-16
m2<-16

id.remove<-c('412707','626494','632811','981252',
             '987762','8710791','9569599','9600937',
             '10182266','10360113','10502933','10632317',
             '10638416','10786900','10921298','10987066',
             '11509405','18892676','19201223','19222498')

#1.Get the matrix x of player specific variables
load("data_0522.RData")
iv<-data_0522
for(ii in 1:length(id.remove)){
  
  idx<-which(iv$characterid==id.remove[ii])
  iv<-iv[-idx,]
}

#2. Get login ind., login time and purchase ind
library(readr)
dv_0522 <- read_csv("dv_0522.txt")
dv<-dv_0522
for(ii in 1:length(id.remove)){
  
  idx<-which(dv$characterid==id.remove[ii])
  dv<-dv[-idx,]
}

#3. Do some clean up

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
  temp2<- matrix(0,m1+m2-1)
  for(j in 2:(m1+m2-1)){
    
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

#5. Start with y (login ind. and login time) and e (puchase ind) first

y.train<-matrix(0,n*(m1-1),1)
e.train<-y.train
y.pred<- matrix(0,n*(m2-1),1)
e.pred<-y.pred
week.train<- y.train
week.pred<-y.pred
holiday.train<-y.train
holiday.pred<-y.pred

for(i in 1:n){
  
  s1<- 2+31*(i-1)
  e1<- 16+31*(i-1)
  s2<- 17+31*(i-1)
  e2<- 31+31*(i-1)
  s<-1+15*(i-1)
  e<-15*i
  y.train[s:e]<- data.full$logintime[s1:e1]
  e.train[s:e]<- data.full$purchasedummy[s1:e1]
  week.train[s:e]<- data.full$week[s1:e1]
  holiday.train[s:e]<- data.full$holiday[s1:e1]
  
  y.pred[s:e]<- data.full$logintime[s2:e2]
  e.pred[s:e]<- data.full$purchasedummy[s2:e2]
  week.pred[s:e]<- data.full$week[s2:e2]
  holiday.pred[s:e]<- data.full$holiday[s2:e2]
}

#6. Now with matrix g of guild variables
load("guildinfodata_0522.RData")

guildinfodata<-rbind(guildinfodata_0522[1:49,],
                     c(50,"2011-01-29",1,50,rep(0,5)),
                     guildinfodata_0522[50:1549,])
K<-50
g1<-matrix(0,K*(m1-1),5)
g2<-matrix(0,K*(m2-1),5)

for(i in 1:K){
  
  temp<-data.frame(guildinfodata[guildinfodata$guildindex==i,
                                 5:9],rownames.force = NA)
  s<-1+15*(i-1)
  e<-15*i
  g1[s:e,]<-cbind(as.numeric(temp[1:15,1]),as.numeric(temp[1:15,2]),
                  as.numeric(temp[1:15,3]),
                  as.numeric(temp[1:15,4]),as.numeric(temp[1:15,5]))
  g2[s:e,]<-cbind(as.numeric(temp[16:30,1]),as.numeric(temp[16:30,2]),
                  as.numeric(temp[16:30,3]),
                  as.numeric(temp[16:30,4]),as.numeric(temp[16:30,5]))
}
g1[is.na(g1)]<-0
g2[is.na(g2)]<-0
rm('guildinfodata_0522','temp','idx')

#6. Now with matrix x and z of player specific variables
unique.id<-unique(data.full$characterid)
x.train<-matrix(0,n*(m1-1),20)
x.pred<-x.train
for(i in 1:n){
  ii<-unique.id[i]
  temp<-data.matrix(data.full[data.full$characterid==ii,
                              c(4,5,7,9:12,15:21,23:25, 38:40)],rownames.force = NA)
  s<-1+15*(i-1)
  e<-15*i
  x.train[s:e,]<-temp[1:15,]
  x.pred[s:e,]<-temp[16:30,]
}
x.train[is.na(x.train)]<-0
x.pred[is.na(x.pred)]<-0
x.train[,18:19]<- cbind(week.train,holiday.train)
x.pred[,18:19]<- cbind(week.pred,holiday.pred)
rm('temp')

#7. Get the Guildmembership matrix (train and pred)
load("data_0522_guildmemship.RData")
guild_mem<- data_0522_guildmemship
for(ii in 1:length(id.remove)){
  
  idx<-which(guild_mem$characterid==id.remove[ii])
  guild_mem<-guild_mem[-idx,]
}

gm.train<-matrix(0,n*(m1-1),K)
gm.pred<-matrix(0,n*(m2-1),K)
unique.id<-unique(guild_mem$characterid)
for(i in 1:n){
  
  ii<-unique.id[i]
  temp<-data.matrix(guild_mem[guild_mem$characterid==ii,
                              7],rownames.force = NA)
  temp<-matrix(temp,ncol=K,byrow = T)
  s<-1+15*(i-1)
  e<-15*i
  gm.train[s:e,]<-temp[1:15,]
  gm.pred[s:e,]<-temp[16:30,]
}

rm('temp','data_0522_guildmemship','s1','e1','s','e',
   's2','e2','idx','ii','i')

## Transformations on y
y.train<-y.train/60
y.pred<-y.train/60

# Transformations on x, z
x.train[,c(1,3:5,7:13,15:17,20)]<-asinh(x.train[,c(1,3:5,7:13,15:17,20)])
x.train[,c(2,6)]<-asinh(x.train[,c(2,6)]/60)
x.pred[,c(1,3:5,7:13,15:17,20)]<-asinh(x.pred[,c(1,3:5,7:13,15:17,20)])
x.pred[,c(2,6)]<-asinh(x.pred[,c(2,6)]/60)
z.train<- x.train[,c(1,2,3:13,15:17,20)]
z.pred<- x.pred[,c(1,2,3:13,15:17,20)]

#Transformations on guild variables
g1<-asinh(g1)
g2<-asinh(g2)

out<-list('g.train'=g1,'g.pred'=g2,
          'z.train'=z.train,'z.pred'=z.pred,
          'x.train'=x.train,'x.pred'=x.pred,
          'y.train'=y.train,'y.pred'=y.pred,
          'e.train'=e.train,'e.pred'=e.pred,
          'gm.train'=gm.train,'gm.pred'=gm.pred,
          'data.full'=data.full)

save(out,file='out.RData')

