#Script for analysing trends in average weight

library(lme4) #mixed models
library(lmerTest)  #for anova for glmm model
#library(AICcmodavg) #for prediction SE for glmm model
library(epiR) #Concordanance between predicted and observed
library(PBSmapping)
data(worldLLhigh)
library(arm)  #form pseudo R2
library(dplyr)


#Standard or First run (paper)
Run="Standard"

#---DATA SECTION-----
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive("Analyses/Catch and effort"))

#1. DAILY LOGBOOKS
Logbook=read.csv("Logbook.data.mean.weight.csv")

Wei.range=read.csv(handl_OneDrive("Data/Length_Weights/Data.Ranges.csv"),stringsAsFactors = F)
Wei.range.names=read.csv(handl_OneDrive("Data/Length_Weights/Species.names.csv"),stringsAsFactors = F)
Wei.range=merge(Wei.range,Wei.range.names,by="Sname",all.x=T)


#2. OBSERVERS DATA
if(Run=="First") Survey=read.csv("Survey.weight.csv")



#----PARAMETERS SECTIONS----

#All parameters sourced from "DoF LW relationships.xlx" in DATA

#FL to TL pars             
b.w=13.171   
b.g=4.6424
b.d=2.9835
b.s=0.2133

a.w=1.0044
a.g=1.0837
a.d=1.1849
a.s=1.2185


#TL to TW pars 
bwt.g=0.0000004623  
bwt.w=0.0000027500     
bwt.d=0.0000034694
bwt.s=0.0000021684
awt.g=3.47701
awt.w=3.08059
awt.d=3.10038
awt.s=3.20688


#MAx FL observed by onboard observers
Max.d=275
Max.s=166
Max.w=146
Max.g=158

Max.FL.w=Max.FL.g=180
Max.FL.d=300
Max.FL.s=240

#Max TWT
Max.w.d=ceiling(bwt.d*(b.d+a.d*Max.d)^awt.d)
Max.w.s=ceiling(bwt.s*(b.s+a.s*Max.s)^awt.s)
Max.w.w=ceiling(bwt.w*(b.w+a.w*Max.w)^awt.w)
Max.w.g=ceiling(bwt.g*(b.g+a.g*Max.g)^awt.g)

#MIN TWT (from size at birth)
Min.w.d=bwt.d*70^awt.d
Min.w.s=bwt.s*55^awt.s
Min.w.w=bwt.w*25^awt.w
Min.w.g=bwt.g*30^awt.g


#Other species
Max.w.wob=subset(Wei.range,Sname=='Wobbegong (general)')$TW.max
Min.w.wob=subset(Wei.range,Sname=='Wobbegong (general)')$TW.min

Max.w.sH=subset(Wei.range,Sname=='Smooth Hammerhead')$TW.max
Min.w.sH=subset(Wei.range,Sname=='Smooth Hammerhead')$TW.min

Max.w.Sp=subset(Wei.range,Sname=='Blacktip, Longnose Grey, Spinner Shark')$TW.max
Min.w.Sp=subset(Wei.range,Sname=='Blacktip, Longnose Grey, Spinner Shark')$TW.min

Max.w.Tig=subset(Wei.range,Sname=='Tiger Shark')$TW.max
Min.w.Tig=subset(Wei.range,Sname=='Tiger Shark')$TW.min


#Species ranges (as per CPUE standardisation)
Dusky.range=c(-28,128)
Sandbar.range=c(-27,117)
Whiskery.range=c(-29,128)
Gummy.range=c(115,129)




#---PROCEDURE SECTIONS----


#1. Create logbook datasets
Logbook=Logbook[-match("X",names(Logbook))]
names(Logbook)[match(c("LatDeg","LongDeg"),names(Logbook))]=c("LAT","LONG")
Logbook=subset(Logbook,!is.na(nfish))
Logbook$Source="Logbook"
Logbook$SHEET_NO=with(Logbook,paste(date,vessel,SNo))

  #Drop vars
DROP=c("TYPE.DATA","Bioregion","conditn","factor","year.c","sname1","Same.return","hooks",
       "LatMin","LongMin","landwt","Same.return.SNo","ID","SNo","DSNo","TSNo","day","block10")
Logbook=Logbook[,-match(DROP,names(Logbook))]

  #Put blocks in same format
Logbook$blockx=substr(Logbook$blockx,1,4)
Logbook$finyear=as.character(Logbook$finyear)
Logbook$vessel=as.character(Logbook$vessel)
Logbook$method=as.character(Logbook$method)
Logbook$zone=as.character(Logbook$zone)
Logbook$FL=NA
Logbook=Logbook[,-match(c("Estuary","Bioregion.old"),names(Logbook))]

  #Core distribution
Mean.w.whiskery=subset(Logbook,LAT <= Whiskery.range[1] & LONG <= Whiskery.range[2])
Mean.w.gummy=subset(Logbook,LONG >= Gummy.range[1]  & LONG <= Gummy.range[2]  & LAT <=(-30))
Mean.w.dusky=subset(Logbook,LAT <= Dusky.range[1]  & LONG <= Dusky.range[2] )
Mean.w.sandbar=subset(Logbook,LAT <= Sandbar.range[1] & LONG <= Sandbar.range[2] )

  #Core dist for other species with data
Other.species=c(19000,18023,18022,18001)
Mean.w.smooth.HH=subset(Logbook,species== 19000 & LONG>=115 & LONG<=128 & LAT<=(-32) & LAT>=(-35))
Mean.w.Spinner=subset(Logbook,species== 18023 & LONG>=113 & LONG<=119 & LAT<=(-27) & LAT>=(-35))
Mean.w.Tiger=subset(Logbook,species== 18022 & LONG>=113 & LONG<=116 & LAT<=(-27) & LAT>=(-34))
Mean.w.Copper=subset(Logbook,species== 18001 & LONG>=115 & LONG<=129 & LAT<=(-31) & LAT>=(-35))



#2. Create survey datasets
if(Run=="First")
{
  #Add weight 
  #Survey$FL=with(Survey,ifelse(is.na(FL)| FL==0,CALCULATED.FL,FL))
  
  Survey$FL=with(Survey,ifelse(is.na(FL)& !is.na(TL) & SPECIES=="WH",(TL-b.w)/a.w,
                               ifelse(is.na(FL)& !is.na(TL) & SPECIES=="GM",(TL-b.g)/a.g,
                                      ifelse(is.na(FL)& !is.na(TL) & SPECIES=="BW",(TL-b.d)/a.d,
                                             ifelse(is.na(FL)& !is.na(TL) & SPECIES=="TK",(TL-b.s)/a.s,FL)))))
  
  Survey=subset(Survey,FL>0 & !is.na(FL))
  Survey$TL=with(Survey,ifelse(SPECIES=="WH",b.w+a.w*FL,
                               ifelse(SPECIES=="GM",b.g+a.g*FL,
                                      ifelse(SPECIES=="BW",b.d+a.d*FL,
                                             ifelse(SPECIES=="TK",b.s+a.s*FL,NA)))))
  Survey$livewt=with(Survey,ifelse(SPECIES=="WH",bwt.w*TL^awt.w,
                                   ifelse(SPECIES=="GM",bwt.g*TL^awt.g,
                                          ifelse(SPECIES=="BW",bwt.d*TL^awt.d,
                                                 ifelse(SPECIES=="TK",bwt.s*TL^awt.s,0)))))
  #keep mesh size 6.5-7
  Survey=subset(Survey,MESH_SIZE%in%c("6.5","7"))
  
  Survey=subset(Survey,!BLOCK==0)
  Survey$nfish=1
  
  #add finyear
  Survey$finyear=with(Survey,ifelse(Month<=6,year-1,year))
  Survey=subset(Survey,!finyear==2013)    #drop 2013/14 due to very few observations
  Dummy=as.character(Survey$finyear+1)
  Dummy=substr(Dummy,3,4)
  Survey$finyear=paste(Survey$finyear,"-",Dummy,sep="")
  Survey$species=with(Survey,ifelse(SPECIES=="WH",17003,ifelse(SPECIES=="GM",17001,
                                                               ifelse(SPECIES=="BW",18003,ifelse(SPECIES=="TK",18007,99999999)))))
  
  Survey$shots=1
  Survey$hours=Survey$SOAK_TIME
  Survey$finyear=as.character(Survey$finyear)
  Survey$BOAT=as.character(Survey$BOAT)
  Survey$Method=as.character(Survey$Method)
  Survey$ZONE=as.character(Survey$zone)
  Survey$ZONE=with(Survey,ifelse(ZONE=="1","Zone1",ifelse(ZONE=="2","Zone2",ifelse(ZONE=="WC","West",NA))))
  Survey$ZONE=as.character(with(Survey,ifelse(is.na(ZONE) & LONG>=116.5 & LAT<=(-26),"Zone2",
                                              ifelse(is.na(ZONE) & LONG<116.5 & LAT<=(-33),"Zone1",
                                                     ifelse(is.na(ZONE) & LAT>(-33) & LAT<=(-26) & LONG<116.5,"West",
                                                            ifelse(is.na(ZONE) & LAT>(-26) & LONG<114,"Closed",
                                                                   ifelse(is.na(ZONE) & LAT>(-26) & LONG>=114 & LONG<123.75,"North",
                                                                          ifelse(is.na(ZONE) & LAT>(-26) & LONG>=123.75,"Joint",ZONE))))))))
  Survey$Source="Observer"
  Survey$bdays=Survey$fdays=Survey$nlines=NA
  Survey$BOTDEPTH=Survey$depthMax
  
  #Drop some vars  
  Keep=match(c("finyear","year","month","date","ZONE","BOTDEPTH","BOAT","fdays","Method","BLOCK","bdays",
               "hours","shots","nlines","NET_LENGTH","species","nfish","livewt","LAT","LONG","Source",
               "SHEET_NO","FL","Mid.Lat","Mid.Long"),names(Survey))
  Survey=Survey[,Keep]
  names(Survey)[1:23]=names(Mean.w.sandbar)
  
  #Core distribution
  Survey.whi=subset(Survey,LAT <= Whiskery.range[1] & LONG <= Whiskery.range[2])
  Survey.whi=subset(Survey.whi,FL>30 & FL <300)
  Survey.gum=subset(Survey,LONG >= Gummy.range[1]  & LONG <= Gummy.range[2]  & LAT <=(-30))
  Survey.gum=subset(Survey.gum,FL>30 & FL <300)
  Survey.dus=subset(Survey,LAT <= Dusky.range[1]  & LONG <= Dusky.range[2] )
  Survey.dus=subset(Survey.dus,FL>30 & FL <300)
  Survey.san=subset(Survey,LAT <= Sandbar.range[1] & LONG <= Sandbar.range[2] )
  Survey.san=subset(Survey.san,FL>30 & FL <300)
  
  
  
  #3. Explore simple mean average weight by year for survey and logbook
  fn=function(dat1,dat2,SPEC,MAX.W,MIN.W,DO.LEG,LAB)
  {
    dat1=subset(dat1,species==SPEC)
    dat2=subset(dat2,species==SPEC)
    dat1$Mean.W=dat1$livewt/dat1$nfish
    dat2$Mean.W=dat2$livewt/dat2$nfish
    dat1=subset(dat1,Mean.W<=MAX.W & !(is.na(Mean.W)|Mean.W<MIN.W))
    dat2=subset(dat2,Mean.W<=MAX.W & !(is.na(Mean.W)|Mean.W<MIN.W))
    ag1=aggregate(Mean.W~year,dat1,mean,na.rm=T)
    ag1.sd=aggregate(Mean.W~year,dat1,sd,na.rm=T)
    
    ag2=aggregate(Mean.W~year,dat2,mean,na.rm=T)
    ag2.sd=aggregate(Mean.W~year,dat2,sd,na.rm=T)
    MIN=min(c(min(ag2$Mean.W-ag2.sd$Mean.W,na.rm=T),min(ag1$Mean.W-ag1.sd$Mean.W,na.rm=T)))
    MAX=max(c(max(ag2$Mean.W+ag2.sd$Mean.W,na.rm=T),max(ag1$Mean.W+ag1.sd$Mean.W,na.rm=T)))  
    YR=sort(unique(c(ag2$year,ag1$year)))
    plot(ag1$year,ag1$Mean.W,pch=19,cex=2,ylim=c(MIN,MAX),xlim=c(min(YR),max(YR)),
         ylab="",xlab="Year",col=4,cex.axis=1.25)
    segments(ag1$year,ag1$Mean.W-ag1.sd$Mean.W,ag1$year,ag1$Mean.W+ag1.sd$Mean.W,col=4)
    points(ag2$year,ag2$Mean.W,pch=19,cex=2,ylim=c(MIN,MAX),col=2)
    segments(ag2$year,ag2$Mean.W-ag2.sd$Mean.W,ag2$year,ag2$Mean.W+ag2.sd$Mean.W,col=2)
    if(DO.LEG=="YES")
    {
      mtext(LAB,3,cex=1.5) 
      legend("topright",c("logbook","survey"),bty='n',pch=19,col=c(2,4),cex=1.75)
    }
    legend("bottomleft",unique(dat1$zone),bty='n',cex=1.75)
    axis(1,YR,F)
    
  }
  
  handle=handl_OneDrive("Analyses/Mean weight/")
  tiff(file=paste(handle,"Mean.data.whiskery.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")  
  par(mfcol=c(3,1),mar=c(2,3.5,1.5,1),oma=c(1,1,1,1),las=1,mgp=c(2,.7,0))
  fn(subset(Mean.w.whiskery,zone=="West"),subset(Survey.whi,zone=="West"),17003,Max.w.w,Min.w.w,"YES","Whiskery shark")
  fn(subset(Mean.w.whiskery,zone=="Zone1"),subset(Survey.whi,zone=="Zone1"),17003,Max.w.w,Min.w.w,"NO",NA)
  fn(subset(Mean.w.whiskery,zone=="Zone2"),subset(Survey.whi,zone=="Zone2"),17003,Max.w.w,Min.w.w,"NO",NA)
  mtext("Year",1,outer=T,line=0,cex=1.5)
  mtext("Live weight +/-SD (kg)",2,outer=T,cex=1.5,las=3,lin=-1)
  dev.off()
  
  
  tiff(file=paste(handle,"Mean.data.dusky.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")  
  par(mfcol=c(3,1),mar=c(2,3.5,1.5,1),oma=c(1,1,1,1),las=1,mgp=c(2,.7,0))
  fn(subset(Mean.w.dusky,zone=="West"),subset(Survey.dus,zone=="West"),18003,Max.w.d,Min.w.d,"YES","Dusky shark")
  fn(subset(Mean.w.dusky,zone=="Zone1"),subset(Survey.dus,zone=="Zone1"),18003,Max.w.d,Min.w.d,"NO",NA)
  fn(subset(Mean.w.dusky,zone=="Zone2"),subset(Survey.dus,zone=="Zone2"),18003,Max.w.d,Min.w.d,"NO",NA)
  mtext("Year",1,outer=T,line=0,cex=1.5)
  mtext("Live weight +/-SD (kg)",2,outer=T,cex=1.5,las=3,lin=-1)
  dev.off()
  
  
  
  tiff(file=paste(handle,"Mean.data.sandbar.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")  
  par(mfcol=c(3,1),mar=c(2,3.5,1.5,1),oma=c(1,1,1,1),las=1,mgp=c(2,.7,0))
  fn(subset(Mean.w.sandbar,zone=="West"),subset(Survey.san,zone=="West"),18007,Max.w.s,Min.w.s,"YES","Sandbar shark")
  fn(subset(Mean.w.sandbar,zone=="Zone1"),subset(Survey.san,zone=="Zone1"),18007,Max.w.s,Min.w.s,"NO",NA)
  fn(subset(Mean.w.sandbar,zone=="Zone2"),subset(Survey.san,zone=="Zone2"),18007,Max.w.s,Min.w.s,"NO",NA)
  mtext("Year",1,outer=T,line=0,cex=1.5)
  mtext("Live weight +/-SD (kg)",2,outer=T,cex=1.5,las=3,lin=-1)
  dev.off()
  
  tiff(file=paste(handle,"Mean.data.gummy.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")  
  par(mfcol=c(1,1),mar=c(2,3.5,1.5,1),oma=c(1,1,1,1),las=1,mgp=c(2,.7,0))
  fn(subset(Mean.w.gummy,zone=="Zone2"),subset(Survey.gum,zone=="Zone2"),17001,Max.w.g,Min.w.g,"YES","Gummy shark")
  mtext("Year",1,outer=T,line=0,cex=1.5)
  mtext("Live weight +/-SD (kg)",2,outer=T,cex=1.5,las=3,lin=-1)
  dev.off()
  
  
  
  
  #4. Mean weight by year and weight distribution by year from commercial logbooks
  fn.Wght.Freq=function(LaT,SPEC,Gear,MAX.w,BRK,XMAX,SPEC.name)
  {
    dat=subset(Logbook,LAT<=(LaT) & species==SPEC & method==Gear)
    dat$Aver.w=dat$livewt/dat$nfish
    yrs=unique(dat$finyear)
    dat=subset(dat,!is.na(Aver.w))
    dat=subset(dat,Aver.w<=MAX.w)
    n.yrs=length(yrs)
    storE=vector('list',n.yrs)
    MeaN=SD=rep(NA,n.yrs)
    yrs1=unique(dat$finyear)
    for (i in 1:n.yrs) 
    {
      dat1=subset(dat,finyear==yrs1[i])
      storE[[i]]=hist(dat1$Aver.w,plot=F,breaks=BRK)
      MeaN[i]=mean(dat1$Aver.w)
      SD[i]=sd(dat1$Aver.w)
    }
    tiff(file=paste(handle,"/Weight.Freq.",SPEC.name,".tiff",sep=""),
         width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")  
    par(mfcol=c(2,1),las=1,mai=c(1,1.25,.1,.1),oma=c(.3,.5,.1,.1),mgp=c(3.25,.7,0))
    COLs=1:n.yrs
    
    #plot histograms
    plot(storE[[1]]$mids,storE[[1]]$counts/sum(storE[[1]]$counts),type='l',xlim=c(0,XMAX),ylim=c(0,1),xaxt='n',
         ylab="Relative frequency",xlab="Average live weight (kg)",cex.lab=1.4,lwd=2)
    for (i in 2:n.yrs) lines(storE[[i]]$mids,storE[[i]]$counts/sum(storE[[i]]$counts),col=COLs[i],lwd=2)
    axis(1,at=BRK,labels=BRK)
    legend("topright",as.character(yrs),lty=1,col=COLs,bty='n',lwd=2)
    
    #plot mean and sd
    plot(1:n.yrs,MeaN,xaxt='n',ylim=c(min(MeaN-SD),max(MeaN+SD)),ylab="Average live weight (+/-SD)",
         xlab="Financial year",cex.lab=1.5,pch=19,cex=2)
    arrows(1:n.yrs, MeaN, 1:n.yrs, MeaN+SD, angle=90, length=0.1)
    arrows(1:n.yrs, MeaN, 1:n.yrs, MeaN-SD, angle=90, length=0.1)
    #add trendline
    dat$finyear1=factor(dat$finyear)
    levels(dat$finyear1)=1:length(levels(dat$finyear1))
    dat$finyear=as.numeric(as.character(dat$finyear1))
    fit <- glm(dat$Aver.w~dat$finyear)
    abline(fit, col="blue",lty=2,  lwd=2)
    axis(1,at=1:n.yrs,labels=yrs)
    text(mean(1:n.yrs),mean(c(min(MeaN),max(MeaN+SD))),paste("slope=",round(coef(fit)[2],2)),
         cex=1.5,col="blue",font=2)
    dev.off()
  }
  handle=handl_OneDrive("Analyses/Mean weight")
  fn.Wght.Freq(-26.5,17003,"GN",MAX.w=Max.w.w,BRK=seq(0,20,by=5),XMAX=Max.w.w,"Whiskery")
  fn.Wght.Freq(-26.5,17001,"GN",MAX.w=Max.w.g,BRK=seq(0,Max.w.g,by=5),XMAX=25,"Gummy")
  fn.Wght.Freq(-26.5,18003,"GN",MAX.w=Max.w.d,BRK=seq(0,Max.w.d,by=5),XMAX=50,"Dusky")
  fn.Wght.Freq(-26.5,18007,"GN",MAX.w=Max.w.s,BRK=seq(0,Max.w.s,by=5),XMAX=25,"Sandbar")
  
  
  #Merge survey and logbook data
  ID=match("FL",names(Survey.whi))
  
  #Merge.survey="YES"   #merge onboard observing on gillnet boats with TDGDLF logbook records
  Merge.survey="NO"
  
  if(Merge.survey=="YES")
  {
    Mean.w.whiskery=rbind(Mean.w.whiskery,Survey.whi)
    Mean.w.gummy=rbind(Mean.w.gummy,Survey.gum)
    Mean.w.dusky=rbind(Mean.w.dusky,Survey.dus)
    Mean.w.sandbar=rbind(Mean.w.sandbar,Survey.san)  
  }
  
}


#Table of numbers by species
if(Run=="First")
{
  TABLE=function(dat1,dat2,SPEC,MAXFL,MAX.W,MIN.W)
  {
    dat=rbind(dat1,dat2)
    dat=subset(dat,species==SPEC)
    dat.logbook=subset(dat,Source=="Logbook")
    dat.observer=subset(dat,Source=="Observer")
    Ns=aggregate(nfish~Source,dat,sum)
    Total.N.survey=Ns$nfish[2]
    Total.N.logbook=Ns$nfish[1]
    Total.N=Total.N.survey+Total.N.logbook
    a=subset(dat,!is.na(FL) & FL<=MAXFL)
    Mean.FL=round(mean(a$FL,na.rm=T))
    SD.FL=round(sd(a$FL,na.rm=T))
    N.measured=sum(a$nfish)
    
    Average.w=dat$livewt/dat$nfish
    Average.w.all=subset(Average.w,!is.na(Average.w))
    Average.w=subset(Average.w.all,Average.w<=MAX.W & Average.w>MIN.W)
    Mean.ave.w=mean(Average.w,na.rm=T)
    SD.ave.w=sd(Average.w,na.rm=T)
    YR.range=range(dat.logbook$finyear)
    YR.range.logbook=paste(YR.range[1],".to.",YR.range[2],sep="")
    
    YR.range=range(dat.observer$finyear)
    YR.range.observer=paste(YR.range[1],".to.",YR.range[2],sep="")
    
    return(list(tab=data.frame(species=SPEC,Total.N=Total.N,Total.N.logbook=Total.N.logbook,
                               Total.N.survey=Total.N.survey,N.measured=N.measured,
                               Mean.FL=Mean.FL,SD.FL=SD.FL,Mean.ave.w=Mean.ave.w,SD.ave.w=SD.ave.w,
                               years.logbook=YR.range.logbook,years.observer=YR.range.observer),
                percent.used=100*length(Average.w)/length(Average.w.all)))
  }
  
  Table1.whi=TABLE(Mean.w.whiskery,Survey.whi[,-(24:25)],17003,Max.FL.w,Max.w.w,Min.w.w)
  Table1.gum=TABLE(Mean.w.gummy,Survey.gum[,-(24:25)],17001,Max.FL.g,Max.w.g,Min.w.g)
  Table1.dus=TABLE(Mean.w.dusky,Survey.dus[,-(24:25)],18003,Max.FL.d,Max.w.d,Min.w.d)
  Table1.san=TABLE(Mean.w.sandbar,Survey.san[,-(24:25)],18007,Max.FL.s,Max.w.s,Min.w.s)
  Table1=rbind(Table1.whi$tab,Table1.gum$tab,Table1.dus$tab,Table1.san$tab)
  
  Used=c(Table1.whi$percent.used,Table1.gum$percent.used,Table1.dus$percent.used,Table1.san$percent.used)
  
  setwd(handle)
  write.csv(Table1,"paper/Table1.csv",row.names=F)
  
}


#Plot overlap of data sets
if(Run=="First")
{
  fn=function(dat,dat1)
  {
    plot(1:10,1:10,col="transparent",ylim=c(-36,-26),xlim=c(112,119))
    points(dat$LONG,dat$LAT,pch=19,col=2)
    points(dat1$Mid.Long,dat1$Mid.Lat,pch=19,col=3)
    legend('topright',c("logbook","survey"),pch=19,col=2:3,bty='n')
  }
  fn(Mean.w.whiskery,Survey.whi)
  fn(Mean.w.gummy,Survey.gum)
  fn(Mean.w.dusky,Survey.dus)
  fn(Mean.w.sandbar,Survey.san)
}


#Create data sets for analysis
data.fn=function(dat,MAX.w,MIN.w)
{
  dat$Mean.wght=dat$livewt/dat$nfish
  dat=subset(dat,Mean.wght<=MAX.w & Mean.wght>MIN.w)
  return(dat)
}
Agg.w.whiskery=data.fn(subset(Mean.w.whiskery,species==17003),Max.w.w,Min.w.w)
Agg.w.gummy=data.fn(subset(Mean.w.gummy,species==17001),Max.w.g,Min.w.g)
Agg.w.dusky=data.fn(subset(Mean.w.dusky,species==18003),Max.w.d,Min.w.d)
Agg.w.sandbar=data.fn(subset(Mean.w.sandbar,species==18007),Max.w.s,Min.w.s)

Agg.w.smooth.HH=data.fn(Mean.w.smooth.HH,Max.w.sH,Min.w.sH)
Agg.w.Spinner=data.fn(Mean.w.Spinner,Max.w.Sp,Min.w.Sp)
Agg.w.Tiger=data.fn(Mean.w.Tiger,Max.w.Tig,Min.w.Tig)
Agg.w.Copper=data.fn(Mean.w.Copper,Max.w.d,Min.w.d)


if(Run=="First")
{
  Agg.w.whiskery.survey=data.fn(subset(Survey.whi,species==17003),Max.w.w,Min.w.w)
  Agg.w.gummy.survey=data.fn(subset(Survey.gum,species==17001),Max.w.g,Min.w.g)
  Agg.w.dusky.survey=data.fn(subset(Survey.dus,species==18003),Max.w.d,Min.w.d)
  Agg.w.sandbar.survey=data.fn(subset(Survey.san,species==18007),Max.w.s,Min.w.s)
}


# #aggregate by shot
# Agg.w.whiskery=aggregate(cbind(nfish,livewt)~finyear+year+month+date+zone+depthMax+vessel
#               +SHEET_NO+blockx+species+Source,subset(Mean.w.whiskery,species==17003),sum)
# 
# Agg.w.gummy=aggregate(cbind(nfish,livewt)~finyear+year+month+date+zone+depthMax+vessel
#                          +SHEET_NO+blockx+species+Source,subset(Mean.w.gummy,species==17001),sum)
# 
# Agg.w.dusky=aggregate(cbind(nfish,livewt)~finyear+year+month+date+zone+depthMax+vessel
#                          +SHEET_NO+blockx+species+Source,subset(Mean.w.dusky,species==18003),sum)
# 
# Agg.w.sandbar=aggregate(cbind(nfish,livewt)~finyear+year+month+date+zone+depthMax+vessel
#                          +SHEET_NO+blockx+species+Source,subset(Mean.w.sandbar,species==18007),sum)
# 
# 
# #Calculate mean weight
# Agg.w.whiskery$Mean.wght=Agg.w.whiskery$livewt/Agg.w.whiskery$nfish
# Agg.w.gummy$Mean.wght=Agg.w.gummy$livewt/Agg.w.gummy$nfish
# Agg.w.dusky$Mean.wght=Agg.w.dusky$livewt/Agg.w.dusky$nfish
# Agg.w.sandbar$Mean.wght=Agg.w.sandbar$livewt/Agg.w.sandbar$nfish
# 
# Agg.w.whiskery=subset(Agg.w.whiskery,Mean.wght<=Max.w.w & Mean.wght>Min.w.w)
# Agg.w.gummy=subset(Agg.w.gummy,Mean.wght<=Max.w.g & Mean.wght>Min.w.g)
# Agg.w.dusky=subset(Agg.w.dusky,Mean.wght<=Max.w.d & Mean.wght>Min.w.d)
# Agg.w.sandbar=subset(Agg.w.sandbar,Mean.wght<=Max.w.s & Mean.wght>Min.w.s)




#check mean weights
fn.check=function(dat,SPEC,MAX.W,MIN.W)
{
  par(mfcol=c(2,1))
  dat=subset(dat,species==SPEC)
  NN=nrow(dat)
  dat.over=subset(dat,Mean.wght>MAX.W | Mean.wght<MIN.W)
  dat.over=subset(dat.over,!is.na(Mean.wght))
  dat=subset(dat,Mean.wght<=MAX.W & !(is.na(Mean.wght)|Mean.wght<MIN.W))
  
  
  boxplot(Mean.wght~as.factor(species),dat)
  hist(dat$Mean.wght)
  a=round(range(dat$Mean.wght,na.rm=T),3)
  legend("right",paste("range=",a[1],"-",a[2],sep=""),bty="n")
  
  dat.over=subset(dat.over,species==SPEC & livewt>MAX.W )
  dat.below=subset(dat.over,species==SPEC & Mean.wght<MIN.W)
  Tab.over=table(dat.over$Source)
  Tab.below=table(dat.below$Source)
  return(list(Total.records=NN,Good.records=nrow(dat),Tab.over=Tab.over,
              Tab.below=Tab.below))
}
Chk.whis=fn.check(Agg.w.whiskery,17003,Max.w.w,Min.w.w)
Chk.gum=fn.check(Agg.w.gummy,17001,Max.w.g,Min.w.g)
Chk.dus=fn.check(Agg.w.dusky,18003,Max.w.d,Min.w.d)
Chk.san=fn.check(Agg.w.sandbar,18007,Max.w.s,Min.w.s)


FIN.YRS=sort(unique(Mean.w.whiskery$finyear))

#Explore
if(Run=="First")
{
  explr.fun=function(dat,SPEC,LAB,MAX.W,MIN.W)
  {
    dat=subset(dat,species==SPEC & Mean.wght<=MAX.W & !(is.na(Mean.wght)|Mean.wght<MIN.W))
    
    
    Mean=aggregate(Mean.wght~year,dat,mean)
    plot(dat$year,dat$Mean.wght,col=1,pch=19,ylab="",xlab="")
    points(Mean$year,Mean$Mean.wght,pch=19,cex=2.5,col=2)
    mtext(LAB,3,cex=1.75)
    
    #tables
    Mn.Yr=table(dat$month,dat$finyear)
    Blk.Yr=table(dat$blockx,dat$finyear)
    return(list(Blk.Yr=Blk.Yr,Mn.Yr=Mn.Yr))
  }
  
  tiff(file="Mean weight by year.logbook.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1,3,2,1),oma=c(1,.5,1,.1),las=1,mgp=c(.1,.7,0))
  EXPL.Whi=explr.fun(Agg.w.whiskery,17003,"Whiskery shark",Max.w.w,Min.w.w)
  EXPL.Gum=explr.fun(Agg.w.gummy,17001,"Gummy shark",Max.w.g,Min.w.g)
  EXPL.Dus=explr.fun(Agg.w.dusky,18003,"Dusky shark",Max.w.d,Min.w.d)
  EXPL.San=explr.fun(Agg.w.sandbar,18007,"Sandbar shark",Max.w.s,Min.w.s)
  mtext("Mean weigth",side=2,line=-1,las=0,cex=2,outer=T)
  dev.off()
  
  
  tiff(file="Mean weight by year.survey.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1,3,2,1),oma=c(1,.5,1,.1),las=1,mgp=c(.1,.7,0))
  EXPL.Whi=explr.fun(Agg.w.whiskery.survey,17003,"Whiskery shark",Max.w.w,Min.w.w)
  EXPL.Gum=explr.fun(Agg.w.gummy.survey,17001,"Gummy shark",Max.w.g,Min.w.g)
  EXPL.Dus=explr.fun(Agg.w.dusky.survey,18003,"Dusky shark",Max.w.d,Min.w.d)
  EXPL.San=explr.fun(Agg.w.sandbar.survey,18007,"Sandbar shark",Max.w.s,Min.w.s)
  mtext("Mean weigth",side=2,line=-1,las=0,cex=2,outer=T)
  dev.off()
  
  
  # #Check tables
  # EXPL.Whi$Blk.Yr
  # EXPL.Whi$Mn.Yr
  # 
  # EXPL.Gum$Blk.Yr
  # EXPL.Gum$Mn.Yr
  # 
  # EXPL.Dus$Blk.Yr
  # EXPL.Dus$Mn.Yr
  # 
  # EXPL.San$Blk.Yr
  # EXPL.San$Mn.Yr
  
  
  #Mean weigth by block for logbooks
  fn.plot.weight.spatially=function(SPEC,MAXW)
  {
    a=subset(Logbook,species==SPEC)
    test=aggregate(cbind(nfish,livewt)~finyear+blockx,a,sum)
    test$Average.W=test$livewt/test$nfish
    test=subset(test,Average.W<=MAXW)
    test$LAT=-(as.numeric(substr(test$blockx,1,2)))
    test$LONG=100+(as.numeric(substr(test$blockx,3,4)))
    
    
    N=sort(unique(a$finyear))
    
    par(mfcol=c(3,3),mai=c(.6,.5,.1,.1))
    for(i in 1:length(N))
    {
      b=subset(test,finyear==N[i])
      plot(b$LONG,b$LAT,cex=b$Average.W/mean(test$Average.W),main=N[i],ylab="",xlab="",pch=19,col=2)
    }
    plot(1:10,1:10,col="transparent",ylab="",xlab="",xaxt='n',yaxt='n')
    legend("center",c("5 kg","10 kg","15 kg"),pch=19,cex=1.5,pt.cex=c(5,10,15)/mean(test$Average.W),bty='n',
           col=2)
    legend("top",paste(SPEC),bty='n',cex=2)
  }
  
  tiff(file="Mean.weight.bk.logbk.whis.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  fn.plot.weight.spatially(17003,Max.w.w)
  dev.off()
  
  tiff(file="Mean.weight.bk.logbk.gum.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  fn.plot.weight.spatially(17001,Max.w.g)
  dev.off()
  
  tiff(file="Mean.weight.bk.logbk.dus.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  fn.plot.weight.spatially(18001,Max.w.d)
  dev.off()
  
  tiff(file="Mean.weight.bk.logbk.san.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  fn.plot.weight.spatially(18007,Max.w.s)
  dev.off()
  
  
  #Coplots
  # coplot(log.weight ~ finyear|blockx, type="b", data=dat) # Lines
  # coplot(log.weight ~ finyear|month, type="b", data=dat) # Points and lines
  
}


#Remove blocks with few observations
Threshold=99  #upper 99 of observations
fn.block.cum=function(dat,SPEC,Threshold,MAX.W,MIN.W)
{
  dat=subset(dat,species==SPEC & Mean.wght<=MAX.W & !(is.na(Mean.wght)|Mean.wght<MIN.W))
  
  dat$dummy=1
  TABLE13=aggregate(dummy~blockx,data=dat,sum)
  TABLE13=TABLE13[order(-TABLE13$dummy),]
  TABLE13$CumCatch=cumsum(TABLE13$dummy)
  TABLE13$PerCumCatch=round(TABLE13$CumCatch*100/sum(TABLE13$dummy),2)
  This.blk=subset(TABLE13,PerCumCatch<=Threshold)
  return(This.blk$blockx)
}
Whis.blks=fn.block.cum(Agg.w.whiskery,17003,Threshold,Max.w.w,Min.w.w)
Gum.blks=fn.block.cum(Agg.w.gummy,17001,Threshold,Max.w.g,Min.w.g)
Dus.blks=fn.block.cum(Agg.w.dusky,18003,Threshold,Max.w.d,Min.w.d)
San.blks=fn.block.cum(Agg.w.sandbar,18007,Threshold,Max.w.s,Min.w.s)

SH.blks=fn.block.cum(Agg.w.smooth.HH,19000,Threshold,Max.w.sH,Min.w.sH)
Spin.blks=fn.block.cum(Agg.w.Spinner,18023,Threshold,Max.w.Sp,Min.w.Sp)
Tig.blks=fn.block.cum(Agg.w.Tiger,18022,Threshold,Max.w.Tig,Min.w.Tig)
Cop.blks=fn.block.cum(Agg.w.Copper,18001,Threshold,Max.w.d,Min.w.d)

 

if(Run=="First")
{
  Whis.blks.survey=c(3415)
  Gum.blks.survey=c(3418,3419)
  Dus.blks.survey=c(3415)
  San.blks.survey=c(3415)
}


#Remove years with spatially uneven sampling (i.e. biased to a particular geographic area)
if(Run=="First")
{
  check.yr=function(dat,SPEC,MAX.W,MIN.W,BLKS)
  {
    dat=subset(dat,species==SPEC & Mean.wght<=MAX.W & blockx%in%BLKS & !(is.na(Mean.wght)|Mean.wght<MIN.W))
    return(list(Yr_Blk=table(dat$finyear,dat$blockx),Yr_Month=table(dat$finyear,dat$month)))
  }
  
  check.yr(Agg.w.whiskery.survey,17003,Max.w.w,Min.w.w,Whis.blks.survey)
  check.yr(Agg.w.gummy.survey,17001,Max.w.g,Min.w.g,Gum.blks.survey)
  check.yr(Agg.w.dusky.survey,18003,Max.w.d,Min.w.d,Dus.blks.survey)
  check.yr(Agg.w.sandbar.survey,18007,Max.w.s,Min.w.s,San.blks.survey)
  
  #Drop these suvery years due to no sampling or very spatially skewed sampling
  no.yrs.whis=c("1999-00","2000-01","2001-02","2003-04","2011-12")
  no.yrs.gum=c("1999-00","1993-94","2003-04","2011-12","2012-13")
  no.yrs.dus=c("1999-00","2003-04")
  no.yrs.san=c("1993-94","1999-00","2000-01","2003-04","2004-05","2011-12")
}



#Hausman test of choosing a fixed or random term
# phtest_glmer <- function (glmerMod, glmMod, ...)
#   {  
#   coef.wi <- coef(glmMod)
#   coef.re <- fixef(glmerMod)  ## changed coef() to fixef() for glmer
#   vcov.wi <- vcov(glmMod)
#   vcov.re <- vcov(glmerMod)
#   names.wi <- names(coef.wi)
#   names.re <- names(coef.re)
#   coef.h <- names.re[names.re %in% names.wi]
#   dbeta <- coef.wi[coef.h] - coef.re[coef.h]
#   df <- length(dbeta)
#   dvcov <- vcov.re[coef.h, coef.h] - vcov.wi[coef.h, coef.h]
#   stat <- abs(t(dbeta) %*% as.matrix(solve(dvcov)) %*% dbeta)  ## added as.matrix()
#   pval <- pchisq(stat, df = df, lower.tail = FALSE)
#   names(stat) <- "chisq"
#   parameter <- df
#   names(parameter) <- "df"
#   alternative <- "one model is inconsistent"
#   res <- list(statistic = stat, p.value = pval, parameter = parameter, 
#               method = "Hausman Test",  alternative = alternative,
#               data.name=deparse(glmerMod@call$data))  ## changed
#   class(res) <- "htest"
#   return(res)
# }
# 
# fn.fit=function(dat)
# {
#   dat=subset(dat,!is.na(depthMax))
#   dat$finyear=as.factor(dat$finyear)
#   dat$month=as.factor(dat$month)
#   dat$blockx=as.factor(dat$blockx)
#   dat$vessel=as.factor(dat$vessel)
#   dat$log.weight=log(dat$Mean.wght)
#   
#   model.ran=lmer(log.weight~finyear+blockx+month+depthMax+(1 |vessel), data = dat)
#   model.fix=glm(log.weight~finyear+blockx+month+depthMax+vessel,family=gaussian, data = dat)
#   
#   #Hausman test
#   Hausman.test=phtest_glmer(model.ran,model.fix)
#   
#   return(Hausman.test)
# }
# 
# Best.model.whi=fn.fit(subset(Agg.w.whiskery,blockx%in%Whis.blks))
# Best.model.gum=fn.fit(subset(Agg.w.gummy,blockx%in%Gum.blks))
# Best.model.dus=fn.fit(subset(Agg.w.dusky,blockx%in%Dus.blks))
# Best.model.san=fn.fit(subset(Agg.w.sandbar,blockx%in%San.blks))  


#Ho is that the random effects model is OK
#Ha is that the results from random effects are different that the fixed effects
#therefore omitting those fixed effects is biasing the estimates of the other fixed
#effects coefficients. 
#Hence, if test is significant use the fixed effects model

#Model selection
fn.best.mod=function(dat,Model.type)
{
  dat=subset(dat,!is.na(depthMax))
  dat$finyear=as.factor(dat$finyear)
  dat$month=as.factor(dat$month)
  dat$blockx=as.factor(dat$blockx)
  dat$vessel=as.factor(dat$vessel)
  dat$log.weight=log(dat$Mean.wght)
  
  
  if(Model.type=="fixed")
  {
    mod1=glm(log.weight~finyear+blockx+month+depthMax+vessel,family=gaussian, data = dat)
    mod2=glm(log.weight~finyear+blockx+month+depthMax,family=gaussian, data = dat)
    mod3=glm(log.weight~finyear+blockx+month,family=gaussian, data = dat)
    mod4=glm(log.weight~finyear+blockx,family=gaussian, data = dat)
    mod5=glm(log.weight~finyear,family=gaussian, data = dat)
    MODS=list(mod1,mod2,mod3,mod4,mod5)
    AICs=rep(NA,length(MODS))
    for(a in 1:length(AICs)) AICs[a]=AIC(MODS[[a]])
    id=which(AICs==min(AICs))
    Best=MODS[[id]]
  }
  
  if(Model.type=="Mixed")
  {
    mod1=lmer(log.weight~finyear+blockx+month+depthMax+(1 |vessel), data = dat)
    mod2=lmer(log.weight~finyear+blockx+month+(1 |vessel), data = dat)
    mod3=lmer(log.weight~finyear+blockx+(1 |vessel), data = dat)
    mod4=lmer(log.weight~finyear+(1 |vessel), data = dat)
    MODS=list(mod1,mod2,mod3,mod4)
    AICs=rep(NA,length(MODS))
    for(a in 1:length(AICs)) AICs[a]=AIC(MODS[[a]])
    id=which(AICs==min(AICs))
    Best=MODS[[id]]
    
  }
  return(Best@call)
}

Whis.fit=fn.best.mod(subset(Agg.w.whiskery,blockx%in%Whis.blks),"Mixed")
Gum.fit=fn.best.mod(subset(Agg.w.gummy,blockx%in%Gum.blks),"Mixed")
Dus.fit=fn.best.mod(subset(Agg.w.dusky,blockx%in%Dus.blks),"Mixed")
San.fit=fn.best.mod(subset(Agg.w.sandbar,blockx%in%San.blks),"Mixed")  

SmH.fit=fn.best.mod(subset(Agg.w.smooth.HH,blockx%in%SH.blks),"Mixed") 
Spi.fit=fn.best.mod(subset(Agg.w.Spinner,blockx%in%Spin.blks),"Mixed") 
Tig.fit=fn.best.mod(subset(Agg.w.Tiger,blockx%in%Tig.blks),"Mixed")
Cop.fit=fn.best.mod(subset(Agg.w.Copper,blockx%in%Cop.blks),"Mixed")



if(Run=="First")
{
  fn.best.mod=function(dat)
  {
    dat=subset(dat,!is.na(depthMax))
    dat$finyear=as.factor(dat$finyear)
    dat$month=as.factor(dat$month)
    dat$vessel=as.factor(dat$vessel)
    dat$log.weight=log(dat$Mean.wght)
    
    mod1=lmer(log.weight~finyear+month+depthMax+(1 |vessel), data = dat)
    mod2=lmer(log.weight~finyear+month+(1 |vessel), data = dat)
    mod3=lmer(log.weight~finyear+(1 |vessel), data = dat)
    MODS=list(mod1,mod2,mod3)
    AICs=rep(NA,length(MODS))
    for(a in 1:length(AICs)) AICs[a]=AIC(MODS[[a]])
    id=which(AICs==min(AICs,na.rm=T))
    Best=MODS[[id]]
    
    return(Best@call)
  }
  Whis.fit.survey=fn.best.mod(subset(Agg.w.whiskery.survey,blockx%in%Whis.blks.survey & !finyear%in%no.yrs.whis))
  Dus.fit.survey=fn.best.mod(subset(Agg.w.dusky.survey,blockx%in%Dus.blks.survey & !finyear%in%no.yrs.dus))
  San.fit.san.survey=fn.best.mod(subset(Agg.w.sandbar.survey,blockx%in%San.blks.survey & !finyear%in%no.yrs.san))
  
  fn.best.mod=function(dat)
  {
    dat=subset(dat,!is.na(depthMax))
    dat$finyear=as.factor(dat$finyear)
    dat$month=as.factor(dat$month)
    dat$vessel=as.factor(dat$vessel)
    dat$log.weight=log(dat$Mean.wght)
    
    mod1=lmer(log.weight~finyear+depthMax+(1 |vessel), data = dat)
    mod2=lmer(log.weight~finyear+(1 |vessel), data = dat)
    MODS=list(mod1,mod2)
    AICs=rep(NA,length(MODS))
    for(a in 1:length(AICs)) AICs[a]=AIC(MODS[[a]])
    id=which(AICs==min(AICs,na.rm=T))
    Best=MODS[[id]]
    
    return(Best@call)
  }
  Gum.fit.survey=fn.best.mod(subset(Agg.w.gummy.survey,blockx%in%Gum.blks.survey & !finyear%in%no.yrs.gum))
  
}



#Fit best model
#Model type
Model.type="Mixed"

  #Logbook
fn.fit=function(dat,BLKS,SPEC,Formula)
{
  dat=subset(dat,blockx%in%BLKS)
  
  dat=subset(dat,!is.na(depthMax))
  
  dat$finyear=as.factor(dat$finyear)
  dat$month=as.factor(dat$month)
  dat$blockx=as.factor(dat$blockx)
  dat$vessel=as.factor(dat$vessel)
  dat$zone=as.factor(dat$zone)
  
  dat$log.weight=log(dat$Mean.wght)
  
#   if(Model.type=="fixed")
#   {
#     model=glm(log.weight~finyear+blockx+month+depthMax+vessel,family=gaussian, data = dat)
#     null.model=glm(log.weight~1,family=gaussian, data = dat)
#   }
  
  Term.devs.mixed=NA
  null.model=lmer(log.weight~1+(1 |vessel), data = dat)
  model=eval(Formula)
  
  if(Run=="First")
  {
    if(SPEC%in%c(17001,17003))
    {
      model=lmer(log.weight~finyear+blockx+month+(1 |vessel), data = dat)
      
      #Terms deviance
      if(Run=="First")
      {
        Finyear.dev=deviance(lmer(log.weight~finyear+(1 |vessel), data = dat))
        Block.dev=deviance(lmer(log.weight~finyear+blockx+(1 |vessel), data = dat))
        Month.dev=deviance(model) 
        Term.devs.mixed=data.frame(Term=c("Null model","finyear","blockx","month"),
                                   deviance=c(deviance(null.model),Finyear.dev,Block.dev,Month.dev))
        
      }
    }
    
    if(SPEC%in%c(18003,18007))
    {
      model=lmer(log.weight~finyear+blockx+month+depthMax+(1 |vessel), data = dat)
      
      #Terms deviance
      if(Run=="First")
      {
        Finyear.dev=deviance(lmer(log.weight~finyear+(1 |vessel), data = dat))
        Block.dev=deviance(lmer(log.weight~finyear+blockx+(1 |vessel), data = dat))
        Month.dev=deviance(lmer(log.weight~finyear+blockx+month+(1 |vessel), data = dat))
        Depth.dev=deviance(model) 
        Term.devs.mixed=data.frame(Term=c("Null model","finyear","blockx","month","depthMax"),
                                   deviance=c(deviance(null.model),Finyear.dev,Block.dev,Month.dev,Depth.dev))
        
      }
    }
  }

  
  #Overall deviance
  if(Run=="First")
  {
    Dev.exp= 100*abs((deviance(null.model) - deviance(model)) / deviance(null.model))
  }
  
  if(Run=="First") return(list(model=model,Dev.exp=Dev.exp,Term.devs.mixed=Term.devs.mixed,dat=dat)) else
    list(model=model,dat=dat)
}

Whis.fit=fn.fit(Agg.w.whiskery,Whis.blks,17003,Whis.fit)  
Gum.fit=fn.fit(Agg.w.gummy,Gum.blks,17001,Gum.fit)
Dus.fit=fn.fit(Agg.w.dusky,Dus.blks,18003,Dus.fit)
San.fit=fn.fit(Agg.w.sandbar,San.blks,18007,San.fit)  

SmH.fit=fn.fit(Agg.w.smooth.HH,SH.blks,19000,SmH.fit) 
Spi.fit=fn.fit(Agg.w.Spinner,Spin.blks,18023,Spi.fit) 
Tig.fit=fn.fit(Agg.w.Tiger,Tig.blks,18022,Tig.fit)
Cop.fit=fn.fit(Agg.w.Copper,Cop.blks,18001,Cop.fit)


    #Dev. explained by model
if(Run=="First")
{
  Whis.fit$Dev.exp
  Gum.fit$Dev.exp
  Dus.fit$Dev.exp
  San.fit$Dev.exp
  
  #Anova
  fn.anova=function(mod,Model.type)
  {
    if(Model.type=="Mixed")
    {
      ANO.fix=anova(mod)  #p values for fixed terms
      ANO.ran=rand(mod)  #p values for random term
      return(list(ANO.fix=ANO.fix,ANO.ran=ANO.ran))
    }
    
    if(Model.type=="fixed") return(anova(mod,test="Chisq"))
  }
  Anova.Whis=fn.anova(Whis.fit$model,"Mixed")
  Anova.Gum=fn.anova(Gum.fit$model,"Mixed")
  Anova.Dus=fn.anova(Dus.fit$model,"Mixed")
  Anova.San=fn.anova(San.fit$model,"Mixed")
  
  
  #Deviance explained by each term
  
  if(Model.type=="Mixed")
  {
    fn.get.dev.exp.term=function(mod)
    {
      
      Null=mod$deviance[1]
      Dev.exp.increase=rep(NA,length(2:nrow(mod)))
      for(q in 2:nrow(mod)) Dev.exp.increase[q-1]=abs(100*(Null-mod$deviance[q])/Null)
      Term.dev.exp=data.frame(Term=mod$Term[2:nrow(mod)])
      Term.dev.exp$Percent.dev.exp=NA
      Term.dev.exp$Percent.dev.exp[1]=Dev.exp.increase[1]
      for(q in 2:nrow(Term.dev.exp)) Term.dev.exp$Percent.dev.exp[q]=Dev.exp.increase[q]-Dev.exp.increase[q-1]
      return(Term.dev.exp)
    }
    
    print("Whiskery")
    print(fn.get.dev.exp.term(Whis.fit$Term.devs.mixed))
    print("--------------------------")
    print("Gummy")
    print(fn.get.dev.exp.term(Gum.fit$Term.devs.mixed))
    print("--------------------------")
    print("Dusky")
    print(fn.get.dev.exp.term(Dus.fit$Term.devs.mixed))
    print("--------------------------")
    print("Sandbar")
    print(fn.get.dev.exp.term(San.fit$Term.devs.mixed))
    
  }
  
  
}

  #Survey
if(Run=="First")
{
  fn.fit.survey=function(dat,SPEC)
  {
    dat=subset(dat,!is.na(depthMax))
    dat$finyear=as.factor(dat$finyear)
    dat$month=as.factor(dat$month)
    dat$vessel=as.factor(dat$vessel)
    dat$zone=as.factor(dat$zone)
    dat$log.weight=log(dat$Mean.wght)
    
    Term.devs.mixed=NA
    null.model=lmer(log.weight~1+(1 |vessel), data = dat)
    if(SPEC%in%c(18007,17003))
    {
      model=lmer(log.weight~finyear+month+(1 |vessel), data = dat)
      
      #Terms deviance
      Finyear.dev=deviance(lmer(log.weight~finyear+(1 |vessel), data = dat))
      Month.dev=deviance(model) 
      Term.devs.mixed=data.frame(Term=c("Null model","finyear","month"),
                                 deviance=c(deviance(null.model),Finyear.dev,Month.dev))
    }
    
    if(SPEC%in%c(18003))
    {
      model=lmer(log.weight~finyear+month+depthMax+(1 |vessel), data = dat)
      
      #Terms deviance
      Finyear.dev=deviance(lmer(log.weight~finyear+(1 |vessel), data = dat))
      Month.dev=deviance(lmer(log.weight~finyear+month+(1 |vessel), data = dat))
      Depth.dev=deviance(model) 
      Term.devs.mixed=data.frame(Term=c("Null model","finyear","month","depthMax"),
                                 deviance=c(deviance(null.model),Finyear.dev,Month.dev,Depth.dev))
    }
    
    if(SPEC%in%c(17001))
    {
      model=lmer(log.weight~finyear+depthMax+(1 |vessel), data = dat)
      
      #Terms deviance
      Finyear.dev=deviance(lmer(log.weight~finyear+(1 |vessel), data = dat))
      Depth.dev=deviance(model) 
      Term.devs.mixed=data.frame(Term=c("Null model","finyear","depthMax"),
                                 deviance=c(deviance(null.model),Finyear.dev,Depth.dev))
    }
    
    #Overall deviance
    Dev.exp= 100*abs((deviance(null.model) - deviance(model)) / deviance(null.model))
    
    return(list(model=model,Dev.exp=Dev.exp,Term.devs.mixed=Term.devs.mixed,dat=dat))
  }
  
  Whis.fit.survey=fn.fit.survey(subset(Agg.w.whiskery.survey,blockx%in%Whis.blks.survey & !finyear%in%no.yrs.whis),17003)
  Gum.fit.survey=fn.fit.survey(subset(Agg.w.gummy.survey,blockx%in%Gum.blks.survey & !finyear%in%no.yrs.gum),17001)
  Dus.fit.survey=fn.fit.survey(subset(Agg.w.dusky.survey,blockx%in%Dus.blks.survey & !finyear%in%no.yrs.dus),18003)
  San.fit.survey=fn.fit.survey(subset(Agg.w.sandbar.survey,blockx%in%San.blks.survey & !finyear%in%no.yrs.san),18007)  
  
  
  #Dev. explained by model
  Whis.fit.survey$Dev.exp
  Gum.fit.survey$Dev.exp
  Dus.fit.survey$Dev.exp
  San.fit.survey$Dev.exp
  
  #Anova
  Anova.Whis=fn.anova(Whis.fit.survey$model,"Mixed")
  Anova.Gum=fn.anova(Gum.fit.survey$model,"Mixed")
  Anova.Dus=fn.anova(Dus.fit.survey$model,"Mixed")
  Anova.San=fn.anova(San.fit.survey$model,"Mixed")
  
  #Deviance explained by each term
  if(Model.type=="Mixed")
  {
    
    print("Whiskery")
    print(fn.get.dev.exp.term(Whis.fit.survey$Term.devs.mixed))
    print("--------------------------")
    print("Gummy")
    print(fn.get.dev.exp.term(Gum.fit.survey$Term.devs.mixed))
    print("--------------------------")
    print("Dusky")
    print(fn.get.dev.exp.term(Dus.fit.survey$Term.devs.mixed))
    print("--------------------------")
    print("Sandbar")
    print(fn.get.dev.exp.term(San.fit.survey$Term.devs.mixed))
    
  }
  
  # fn.dev.exp.term=function(dat,mod,Model.type)
  # {
  #   if(Model.type=="fixed")return(100*(dat$Deviance[2:length(dat$Deviance)]/mod$null.deviance))
  #   
  #   if(Model.type=="Mixed")
  #   {
  #     a=mod$deviance[2:length(mod$deviance)]
  #     names(a)=mod$Term[2:length(mod$deviance)]
  #     for (i in 1:length(a))a[i]=max(abs(mod$deviance[1]-a[i]),0)
  #     
  #     return(100*(a/abs(mod$deviance[i])))
  #   }
  # }
  
  # if(Model.type=="fixed")
  # {
  #   NN=2:length(rownames(Anova.Whis))
  #   DaT=data.frame(Term=rownames(Anova.Whis)[NN])
  #   DaT.whi=DaT.gum=DaT.san=DaT.dus=DaT
  #   
  #   DaT.whi$Term.Dev.exp=fn.dev.exp.term(Anova.Whis,Whis.fit$model,"fixed")
  #   DaT.whi$p=Anova.Whis$"Pr(>Chi)"[NN]
  #   DaT.whi$Tot.dev.exp=Whis.fit$Dev.exp
  #   write.csv(DaT.whi,"paper/Anova.whi.csv",row.names=F)
  #   
  #   DaT.gum$Term.Dev.exp=fn.dev.exp.term(Anova.Gum,Gum.fit$model,"fixed")
  #   DaT.gum$p=Anova.Gum$"Pr(>Chi)"[NN]
  #   DaT.gum$Tot.dev.exp=Gum.fit$Dev.exp
  #   write.csv(DaT.gum,"paper/Anova.gum.csv",row.names=F)
  #   
  #   DaT.dus$Term.Dev.exp=fn.dev.exp.term(Anova.Dus,Dus.fit$model,"fixed")
  #   DaT.dus$p=Anova.Dus$"Pr(>Chi)"[NN]
  #   DaT.dus$Tot.dev.exp=Dus.fit$Dev.exp
  #   write.csv(DaT.dus,"paper/Anova.dus.csv",row.names=F)
  #   
  #   DaT.san$Term.Dev.exp=fn.dev.exp.term(Anova.San,San.fit$model,"fixed")
  #   DaT.san$p=Anova.San$"Pr(>Chi)"[NN]
  #   DaT.san$Tot.dev.exp=San.fit$Dev.exp
  #   write.csv(DaT.san,"paper/Anova.san.csv",row.names=F)
  # }
}

#Model fit diagnostic plots
if(Run=="First")
{
  fn.plot.diag=function(MODEL,SPECIES,Model.type)
  {
    if(Model.type=="fixed")RES=MODEL$residuals   #residuals
    if(Model.type=="Mixed")RES=resid(MODEL)
    Std.RES=RES/sd(RES)   #standardised residuals (res/SD(res))
    PREDS=predict(MODEL)
    
    par(mfcol=c(2,2),las=1,mar=c(3,3,2,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(2,.5,0),cex.axis=.8,cex.lab=1.1)
    qqnorm(RES,main="",ylim=c(-5,5),xlim=c(-5,5),ylab="Residuals",xlab="Quantiles of standard normal distribution")
    qqline(RES, col = 'grey40',lwd=1.5,lty=2)
    
    hist(Std.RES,xlim=c(-5,5),ylab="Frequency",xlab="Stan. residuals",main="",col="grey",breaks=50)
    box()
    
    plot(PREDS,Std.RES,ylim=c(-4,12),ylab="Stan. residuals",xlab="Expected values")
    abline(0,0,lwd=1.5,lty=2,col='grey40')
    
    plot(PREDS,sqrt(abs(Std.RES)),ylim=c(0,2.6),ylab="Square root of stan. residuals",xlab="Expected values")
    mtext(paste(SPECIES),3,outer=T,lin=-1)
  }
  #Logbook
  fn.plot.diag(Whis.fit$model,"Whiskery shark","Mixed")
  fn.plot.diag(Gum.fit$model,"gummy shark","Mixed")
  fn.plot.diag(Dus.fit$model,"dusky shark","Mixed")
  fn.plot.diag(San.fit$model,"sandbar shark","Mixed")
  
  
  #survey
  fn.plot.diag(Whis.fit.survey$model,"Whiskery shark","Mixed")
  fn.plot.diag(Gum.fit.survey$model,"gummy shark","Mixed")
  fn.plot.diag(Dus.fit.survey$model,"dusky shark","Mixed")
  fn.plot.diag(San.fit.survey$model,"sandbar shark","Mixed")
  
}

#Plot observations of mean annual weight
if(Run=="First")
{
  fn.plot.observed=function(dat)
  {
    MeaN=aggregate(Mean.wght~finyear,dat,mean)
    SD=aggregate(Mean.wght~finyear,dat,sd)
    yrs=MeaN$finyear
    n.yrs=length(yrs)  
    MeaN=MeaN$Mean.wght
    SD=SD$Mean.wght
    
    plot(1:n.yrs,MeaN,xaxt='n',ylim=c(min(MeaN-SD),max(MeaN+SD)),ylab="Average live weight (+/-SD)",
         xlab="Financial year",cex.lab=1.5,pch=19,cex=2)
    arrows(1:n.yrs, MeaN, 1:n.yrs, MeaN+SD, angle=90, length=0.1)
    arrows(1:n.yrs, MeaN, 1:n.yrs, MeaN-SD, angle=90, length=0.1)
    axis(1,at=1:n.yrs,labels=yrs)
    
  }
  #Logbook
  par(mfcol=c(2,2),mar=c(1,4,2,.1),oma=c(3,1,.1,.1),las=1,mgp=c(.1,.9,0))
  fn.plot.observed(Whis.fit$dat)
  fn.plot.observed(Gum.fit$dat)
  fn.plot.observed(Dus.fit$dat)
  fn.plot.observed(San.fit$dat)
  
  #Survey
  par(mfcol=c(2,2),mar=c(1,4,2,.1),oma=c(3,1,.1,.1),las=1,mgp=c(.1,.9,0))
  fn.plot.observed(Whis.fit.survey$dat)
  fn.plot.observed(Gum.fit.survey$dat)
  fn.plot.observed(Dus.fit.survey$dat)
  fn.plot.observed(San.fit.survey$dat)
  
}


# Predicted annual mean weights

#Blocks by zone
Blk.zn=Logbook%>%select(zone,blockx)%>%
                 distinct(blockx,.keep_all =T)%>%
                filter(zone%in%c('West','Zone1','Zone2'))

#create species data
  #Logbooks
Pred.dat.w=subset(Agg.w.whiskery,blockx%in%Whis.blks)
Pred.dat.g=subset(Agg.w.gummy,blockx%in%Gum.blks)
Pred.dat.d=subset(Agg.w.dusky,blockx%in%Dus.blks)
Pred.dat.s=subset(Agg.w.sandbar,blockx%in%San.blks)

N.w=table(Pred.dat.w$finyear)
N.g=table(Pred.dat.g$finyear)
N.d=table(Pred.dat.d$finyear)
N.s=table(Pred.dat.s$finyear)

Pred.dat.smh=subset(Agg.w.smooth.HH,blockx%in%SH.blks)
Pred.dat.spi=subset(Agg.w.Spinner,blockx%in%Spin.blks)
Pred.dat.tig=subset(Agg.w.Tiger,blockx%in%Tig.blks)
Pred.dat.cop=subset(Agg.w.Copper,blockx%in%Cop.blks)

N.smh=table(Pred.dat.smh$finyear)
N.spi=table(Pred.dat.spi$finyear)
N.tig=table(Pred.dat.tig$finyear)
N.cop=table(Pred.dat.cop$finyear)


  #Survey
if(Run=="First")
{
  Pred.dat.w.survey=subset(Agg.w.whiskery.survey,blockx%in%Whis.blks.survey & !finyear%in%no.yrs.whis)
  Pred.dat.g.survey=subset(Agg.w.gummy.survey,blockx%in%Gum.blks.survey & !finyear%in%no.yrs.gum)
  Pred.dat.d.survey=subset(Agg.w.dusky.survey,blockx%in%Dus.blks.survey & !finyear%in%no.yrs.dus)
  Pred.dat.s.survey=subset(Agg.w.sandbar.survey,blockx%in%San.blks.survey & !finyear%in%no.yrs.san)
  N.w.survey=table(Pred.dat.w.survey$finyear)
  N.g.survey=table(Pred.dat.g.survey$finyear)
  N.d.survey=table(Pred.dat.d.survey$finyear)
  N.s.survey=table(Pred.dat.s.survey$finyear)
  
}

  #observations VS predictions
if(Run=="First")
{
  fn.obs.pred=function(dat,MODEL,SPECIES)
  {
    
    PREDS=predict(MODEL)
    OBS=log(dat$Mean.wght)
    OBS=OBS[1:length(PREDS)]
    
    plot(OBS,PREDS,ylab="Expected value",xlab="Observed value")
    lines(OBS,OBS,col=2,lwd=2)
    mtext(paste(SPECIES),3,line=1,cex=2)
    
    #Concordanance between predicted and observed
    Concordance=c(rho=epi.ccc(OBS, PREDS)$rho.c,C.b=epi.ccc(OBS, PREDS)$C.b) 
    return(Concordance)
  }
  fn.obs.pred(Pred.dat.w,Whis.fit$model,"Whiskery shark")
  fn.obs.pred(Pred.dat.g,Gum.fit$model,"gummy shark")
  fn.obs.pred(Pred.dat.d,Dus.fit$model,"dusky shark")
  fn.obs.pred(Pred.dat.s,San.fit$model,"sandbar shark")
  
  #Concordance
  #No deviation from the 45 degree line occurs when C.b = 1
  
}


  #function for creating new data set for predicting
fn.new.dat=function(dat,vars)
{
  dat=subset(dat,!is.na(depthMax))
  
  dat$finyear=as.factor(dat$finyear)
  dat$month=as.factor(dat$month)
  dat$blockx=as.factor(dat$blockx)
  dat$vessel=as.factor(dat$vessel)

  finyear=levels(dat$finyear)
  
  Mode=table(dat$vessel)  
  vessel=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$vessel))
  
  Mode=table(dat$month)  
  month=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$month))
  
  Mode=table(dat$blockx)  
  blockx=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$blockx))
  
  depthMax=mean(dat$depthMax)
  NEW=data.frame(finyear=finyear,month=month,vessel=vessel,blockx=blockx,depthMax=depthMax)
  NEW=NEW[,match(vars,names(NEW))]
  return(NEW)
  
  
}

  #function for plotting confidence bounds
CI.fun=function(YR,LOW1,UP1,Colr,Colr2)
{
  Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
  Biom.Vec <- c(LOW1, tail(UP1, 1), rev(UP1), LOW1[1])
  polygon(Year.Vec, Biom.Vec, col = Colr, border = Colr2)
}

  #function for predicting
fn.pred.mixed=function(model,NEWDATA,N)
{
  #fitted model
  fm1=model
  
  #model matrix
  mm = model.matrix(terms(fm1),NEWDATA)
  
  #mean prediction
  Log.pred=mm %*% fixef(fm1)
  
  
  #variance
  pvar1 <- diag(mm %*% tcrossprod(vcov(fm1),mm))  #fixed effect only
  pvar1=pvar1+VarCorr(fm1)$vessel[1]     #adding random effect variance
  SD=(pvar1)^0.5
  SE=SD/(c(N))^0.5
  return(list(Pred.mean=Log.pred,Pred.SE=SE))
}

Predict.fn=function(model,NEWdata,Model.type,N)
{
  if(Model.type=="fixed")
  {
    PREDS=predict(model,newdata=NEWdata,type='response',se.fit=T)
    Pred.mean=exp(PREDS$fit+(PREDS$se.fit^2)/2)    #bias correction for log transformation
    Pred.SE=exp(PREDS$se.fit)
   }
  
  if(Model.type=="Mixed")
  {
    #new data
    NEWdata$log.weight=0
    PRED.MIX=fn.pred.mixed(model,NEWdata,N)
    Pred.mean=exp(PRED.MIX$Pred.mean+(PRED.MIX$Pred.SE^2)/2)
    Pred.SE=exp(PRED.MIX$Pred.SE)
  }
    
  PRED=data.frame(Finyear=NEWdata$finyear,Pred.mean=Pred.mean,Pred.SE=Pred.SE,CV=Pred.SE/Pred.mean)
  return(PRED)
}

  #predict new data

    #Logbook
#w.vars=g.vars=c("finyear","blockx","month","vessel")
#d.vars=s.vars=c("finyear","blockx","month","vessel","depthMax")

le.vars=c("finyear","blockx","month","vessel","depthMax")
w.vars=g.vars=d.vars=s.vars=le.vars
New.w=fn.new.dat(Pred.dat.w,w.vars)
New.g=fn.new.dat(Pred.dat.g,g.vars)
New.d=fn.new.dat(Pred.dat.d,d.vars)
New.s=fn.new.dat(Pred.dat.s,s.vars)

New.smh=fn.new.dat(Pred.dat.smh,le.vars)
New.spi=fn.new.dat(Pred.dat.spi,le.vars)
New.tig=fn.new.dat(Pred.dat.tig,le.vars)
New.cop=fn.new.dat(Pred.dat.cop,le.vars)

Pred.whis=Predict.fn(Whis.fit$model,New.w,"Mixed",N.w)
Pred.gum=Predict.fn(Gum.fit$model,New.g,"Mixed",N.g)
Pred.dus=Predict.fn(Dus.fit$model,New.d,"Mixed",N.d)
Pred.san=Predict.fn(San.fit$model,New.s,"Mixed",N.s)

Pred.smh=Predict.fn(SmH.fit$model,New.smh,"Mixed",N.smh)
Pred.spi=Predict.fn(Spi.fit$model,New.spi,"Mixed",N.spi)
Pred.tig=Predict.fn(Tig.fit$model,New.tig,"Mixed",N.tig)
Pred.cop=Predict.fn(Cop.fit$model,New.cop,"Mixed",N.cop)

 
 #add sample size
Pred.smh$n=SmH.fit$dat%>%group_by(finyear)%>%summarise(n=n())%>%pull(n)
Pred.spi$n=Spi.fit$dat%>%group_by(finyear)%>%summarise(n=n())%>%pull(n)
Pred.tig$n=Tig.fit$dat%>%group_by(finyear)%>%summarise(n=n())%>%pull(n)
Pred.cop$n=Cop.fit$dat%>%group_by(finyear)%>%summarise(n=n())%>%pull(n)


  #Predict by zone 
  #Whiskery
dummy=fn.fit(subset(Agg.w.whiskery,zone=="West"),Whis.blks,17003,Whis.fit$model@call)
Pred.whis_west=Predict.fn(dummy$model,fn.new.dat(subset(Pred.dat.w,zone=="West"),w.vars),
                          "Mixed",with(subset(Pred.dat.w,zone=="West"),table(finyear)))

dummy=fn.fit(subset(Agg.w.whiskery,zone=="Zone1"),Whis.blks,17003,Whis.fit$model@call)
Pred.whis_zn1=Predict.fn(dummy$model,fn.new.dat(subset(Pred.dat.w,zone=="Zone1"),w.vars),
                         "Mixed",with(subset(Pred.dat.w,zone=="Zone1"),table(finyear)))

dummy=fn.fit(subset(Agg.w.whiskery,zone=="Zone2"),Whis.blks,17003,Whis.fit$model@call)
Pred.whis_zn2=Predict.fn(dummy$model,fn.new.dat(subset(Pred.dat.w,zone=="Zone2"),w.vars),
                         "Mixed",with(subset(Pred.dat.w,zone=="Zone2"),table(finyear)))

  #Gummy
Pred.gum_zn2=Pred.gum  #only zone 2 data for gummy

  #Dusky
dummy=fn.fit(subset(Agg.w.dusky,zone=="West"),Dus.blks,18003,Dus.fit$model@call)
Pred.dus_west=Predict.fn(dummy$model,fn.new.dat(subset(Pred.dat.d,zone=="West"),d.vars),
                          "Mixed",with(subset(Pred.dat.d,zone=="West"),table(finyear)))

dummy=fn.fit(subset(Agg.w.dusky,zone=="Zone1"),Dus.blks,18003,Dus.fit$model@call)
Pred.dus_zn1=Predict.fn(dummy$model,fn.new.dat(subset(Pred.dat.d,zone=="Zone1"),d.vars),
                         "Mixed",with(subset(Pred.dat.d,zone=="Zone1"),table(finyear)))

dummy=fn.fit(subset(Agg.w.dusky,zone=="Zone2"),Dus.blks,18003,Dus.fit$model@call)
Pred.dus_zn2=Predict.fn(dummy$model,fn.new.dat(subset(Pred.dat.d,zone=="Zone2"),d.vars),
                        "Mixed",with(subset(Pred.dat.d,zone=="Zone2"),table(finyear)))



  #Sandbar
dummy=fn.fit(subset(Agg.w.sandbar,zone=="West"),San.blks,18007,San.fit$model@call)
Pred.san_west=Predict.fn(dummy$model,fn.new.dat(subset(Pred.dat.s,zone=="West"),s.vars),
                         "Mixed",with(subset(Pred.dat.s,zone=="West"),table(finyear)))

dummy=fn.fit(subset(Agg.w.sandbar,zone=="Zone1"),San.blks,18007,San.fit$model@call)
Pred.san_zn1=Predict.fn(dummy$model,fn.new.dat(subset(Pred.dat.s,zone=="Zone1"),s.vars),
                         "Mixed",with(subset(Pred.dat.s,zone=="Zone1"),table(finyear)))

rm(dummy)

    #Survey
if(Run=="First")
{
  w.vars=s.vars=c("finyear","month","vessel")
  d.vars=c("finyear","month","vessel","depthMax")
  g.vars=c("finyear","vessel","depthMax")
  New.w=fn.new.dat(Pred.dat.w.survey,w.vars)
  New.g=fn.new.dat(Pred.dat.g.survey,g.vars)
  New.d=fn.new.dat(Pred.dat.d.survey,d.vars)
  New.s=fn.new.dat(Pred.dat.s.survey,s.vars)
  
  Pred.whis.survey=Predict.fn(Whis.fit.survey$model,New.w,"Mixed",N.w.survey)
  Pred.gum.survey=Predict.fn(Gum.fit.survey$model,New.g,"Mixed",N.g.survey)
  Pred.dus.survey=Predict.fn(Dus.fit.survey$model,New.d,"Mixed",N.d.survey)
  Pred.san.survey=Predict.fn(San.fit.survey$model,New.s,"Mixed",N.s.survey)
  
}

  #Plot predictions
YRS=1:length(FIN.YRS)

    #year as continuous variable
plot.pred.wgt=function(dat,Yr.plt,normalise)
{
  dat$Finyear=as.character(dat$Finyear)
  add.yrs=FIN.YRS[which(!FIN.YRS%in%unique(dat$Finyear))]
  if(length(add.yrs)>0)dat=rbind(dat,data.frame(Finyear=add.yrs,Pred.mean=NA,Pred.SE=NA))
  dat=dat[order(dat$Finyear),]
  
  dat$Pred.mean=exp(dat$Pred.mean)
  dat$Pred.SE=exp(dat$Pred.SE)
  
  if(normalise=="NO")
  {
    std.mean=dat$Pred.mean
    UP1=(dat$Pred.mean)+(dat$Pred.SE)*2
    LOW1=(dat$Pred.mean)-(dat$Pred.SE)*2
  }
  if(normalise=="YES")
  {
    MEAN=mean(dat$Pred.mean,na.rm=T)
    std.mean=dat$Pred.mean/MEAN
    UP1=(dat$Pred.mean/MEAN)+(dat$Pred.SE/MEAN)*2
    LOW1=(dat$Pred.mean/MEAN)-(dat$Pred.SE/MEAN)*2
  }
  
  
  YLIMs=c(min(LOW1,na.rm=T),max(UP1,na.rm=T))
  #YLIMs=c(0,2)

  plot(YRS,YRS,pch=19,cex=1.5,col="transparent",ylim=YLIMs,ylab="",xlab="",cex.axis=1.5,xaxt='n')
  for(s in 1:length(Yr.plt))
  {
    n=Yr.plt[[s]]
    CI.fun(YRS[n],UP1[n],LOW1[n],"grey75","transparent")
    lines(YRS[n],std.mean[n],lwd=2.5,lty=1)
    if(length(n)==1)
    {
      lines(c(YRS[n],YRS[n]),c(LOW1[n],UP1[n]),lwd=8,col="grey75")
      points(YRS[n],std.mean[n],pch="-",cex=1.5)
    }
  }
  axis(1,YRS,F)
  axis(1,seq(1,length(YRS),4),F,tck=-0.06)
}

    #year as dotpoints
plot.pred.wgt_dot=function(dat,Yr.plt,normalise)
{
  dat$Finyear=as.character(dat$Finyear)
  add.yrs=FIN.YRS[which(!FIN.YRS%in%unique(dat$Finyear))]
  if(length(add.yrs)>0)dat=rbind(dat,data.frame(Finyear=add.yrs,Pred.mean=NA,Pred.SE=NA))
  dat=dat[order(dat$Finyear),]
  
  dat$Pred.mean=exp(dat$Pred.mean)
  dat$Pred.SE=exp(dat$Pred.SE)
  
  if(normalise=="NO")
  {
    std.mean=dat$Pred.mean
    UP1=(dat$Pred.mean)+(dat$Pred.SE)*2
    LOW1=(dat$Pred.mean)-(dat$Pred.SE)*2
  }
  if(normalise=="YES")
  {
    MEAN=mean(dat$Pred.mean,na.rm=T)
    std.mean=dat$Pred.mean/MEAN
    UP1=(dat$Pred.mean/MEAN)+(dat$Pred.SE/MEAN)*2
    LOW1=(dat$Pred.mean/MEAN)-(dat$Pred.SE/MEAN)*2
  }
  YLIMs=c(min(LOW1,na.rm=T),max(UP1,na.rm=T))
  plot(YRS,std.mean,pch=19,ylim=YLIMs,ylab="",xlab="",cex.axis=1.5,xaxt='n',cex=2)
  arrows(YRS, std.mean, YRS, UP1, angle=90, length=0.1)
  arrows(YRS, std.mean, YRS, LOW1, angle=90, length=0.1)
  axis(1,YRS,F)
  axis(1,seq(1,length(YRS),4),F,tck=-0.06)
}

AXIS1=function()axis(1,seq(1,length(YRS),2),FIN.YRS[seq(1,length(YRS),2)],tck=-0.05,cex.axis=1.25)
Yr.plt.g=list(YRS)
Yr.plt.d=list(YRS)
Yr.plt.w=list(YRS)
Yr.plt.s=list(YRS)

#Logbook
  #year as continuous var
if(Run=="First")
{
  tiff(file="paper/pred Mean weight by finyr.logbook.tiff",width = 2400, height = 2000,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1,4,2,.1),oma=c(3,1,.1,.1),las=1,mgp=c(.1,.9,0))
  plot.pred.wgt(Pred.gum,Yr.plt.g,"YES")
  mtext("Gummy shark",3,cex=1.25)
  plot.pred.wgt(Pred.dus,Yr.plt.d,"YES")
  mtext("Dusky shark",3,cex=1.25)
  AXIS1()
  plot.pred.wgt(Pred.whis,Yr.plt.w,"YES")
  mtext("Whiskery shark",3,cex=1.25)
  plot.pred.wgt(Pred.san,Yr.plt.s,"YES")
  mtext("Sandbar shark",3,cex=1.25)
  mtext("Financial year",1,line=1.5,cex=1.5,outer=T)
  mtext("Normalised mean weight",2,line=-1,cex=1.5,outer=T,las=3)
  AXIS1()
  dev.off()
}


  #year as factor (i.e. dotpoints)
if(Run=="First")
{
  tiff(file="paper/pred Mean weight by finyr.logbook_dots.tiff",width = 2400, height = 2000,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1,4,2,.1),oma=c(3,1,.1,.1),las=1,mgp=c(.1,.9,0))
  plot.pred.wgt_dot(Pred.gum,Yr.plt.g,"YES")
  mtext("Gummy shark",3,cex=1.25)
  plot.pred.wgt_dot(Pred.dus,Yr.plt.d,"YES")
  mtext("Dusky shark",3,cex=1.25)
  AXIS1()
  plot.pred.wgt_dot(Pred.whis,Yr.plt.w,"YES")
  mtext("Whiskery shark",3,cex=1.25)
  plot.pred.wgt_dot(Pred.san,Yr.plt.s,"YES")
  mtext("Sandbar shark",3,cex=1.25)
  mtext("Financial year",1,line=1.5,cex=1.5,outer=T)
  mtext("Normalised mean weight",2,line=-1,cex=1.5,outer=T,las=3)
  AXIS1()
  dev.off()
}

#Survey
if(Run=="First")
{
  FIN.YRS=sort(c(unique(Survey$finyear),"2007-08","2008-09","2009-10","2010-11"))
  YRS=1:length(FIN.YRS)
  Yr.plt.g=list(2:6)
  Yr.plt.d=list(1:6,8:10,12:14,19:20)
  Yr.plt.w=list(1:6,10,12:14,20)
  Yr.plt.s=list(2:6,9:10,13:14,20)
  AXIS1=function()axis(1,seq(1,length(YRS),4),FIN.YRS[seq(1,length(YRS),4)],tck=-0.05,cex.axis=1.25)
  
  #year as continuous var
  tiff(file="paper/pred Mean weight by finyr.survey.tiff",width = 2400, height = 2000,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1,4,2,.1),oma=c(3,1,.1,.1),las=1,mgp=c(.1,.9,0))
  plot.pred.wgt(Pred.gum.survey,Yr.plt.g,"YES")
  mtext("Gummy shark",3,cex=1.25)
  plot.pred.wgt(Pred.dus.survey,Yr.plt.d,"YES")
  mtext("Dusky shark",3,cex=1.25)
  AXIS1()
  plot.pred.wgt(Pred.whis.survey,Yr.plt.w,"YES")
  mtext("Whiskery shark",3,cex=1.25)
  plot.pred.wgt(Pred.san.survey,Yr.plt.s,"YES")
  mtext("Sandbar shark",3,cex=1.25)
  mtext("Financial year",1,line=1.5,cex=1.5,outer=T)
  mtext("Normalised mean weight",2,line=-1,cex=1.5,outer=T,las=3)
  AXIS1()
  dev.off()
  
  
  #year as factor (i.e. dotpoints)
  tiff(file="paper/pred Mean weight by finyr.survey_dots.tiff",width = 2400, height = 2000,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(1,4,2,.1),oma=c(3,1,.1,.1),las=1,mgp=c(.1,.9,0))
  plot.pred.wgt_dot(Pred.gum.survey,Yr.plt.g,"YES")
  mtext("Gummy shark",3,cex=1.25)
  plot.pred.wgt_dot(Pred.dus.survey,Yr.plt.d,"YES")
  mtext("Dusky shark",3,cex=1.25)
  AXIS1()
  plot.pred.wgt_dot(Pred.whis.survey,Yr.plt.w,"YES")
  mtext("Whiskery shark",3,cex=1.25)
  plot.pred.wgt_dot(Pred.san.survey,Yr.plt.s,"YES")
  mtext("Sandbar shark",3,cex=1.25)
  mtext("Financial year",1,line=1.5,cex=1.5,outer=T)
  mtext("Normalised mean weight",2,line=-1,cex=1.5,outer=T,las=3)
  AXIS1()
  dev.off()
  
}


#Effect of blocks
if(Run=="First")
{
  fn.new.dat=function(dat)
  {
    dat=subset(dat,!is.na(depthMax))
    
    dat$finyear=as.factor(dat$finyear)
    dat$month=as.factor(dat$month)
    dat$blockx=as.factor(dat$blockx)
    dat$vessel=as.factor(dat$vessel)
    
    Mode=table(dat$finyear)  
    finyear=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$finyear))
    
    
    Mode=table(dat$vessel)  
    vessel=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$vessel))
    
    Mode=table(dat$month)  
    month=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$month))
    
    blockx=levels(dat$blockx)
    
    depthMax=mean(dat$depthMax)
    
    return(data.frame(finyear=finyear,month=month,vessel=vessel,blockx=blockx,depthMax=depthMax))
    
    
  }
  
  #function for predicting
  Predict.fn.blk=function(model,NEWdata,Model.type)
  {
    if(Model.type=="fixed")
    {
      biasCorr <- model$deviance/model$df.residual/2      
      PREDS=predict(model,newdata=NEWdata,type='response',se.fit=T)
      
      Pred.mean=exp(PREDS$fit+biasCorr)
      Pred.SE=exp(PREDS$se.fit+biasCorr)      
    }
    
    if(Model.type=="Mixed") Pred.mean=exp(predict(model,newdata=NEWdata,type='response'))
    
    Log.pred=data.frame(Block=NEWdata$blockx,Pred.mean=Pred.mean)
    return(Log.pred)
  }
  
  
  #do predictions
  Pred.whis=Predict.fn.blk(Whis.fit$model,fn.new.dat(Pred.dat.w),"Mixed")
  Pred.gum=Predict.fn.blk(Gum.fit$model,fn.new.dat(Pred.dat.g),"Mixed")
  Pred.dus=Predict.fn.blk(Dus.fit$model,fn.new.dat(Pred.dat.d),"Mixed")
  Pred.san=Predict.fn.blk(San.fit$model,fn.new.dat(Pred.dat.s),"Mixed")
  
  
  
  #Spatial display
  fun.plot.map=function(what,LEG)
  {
    
    LAT=-(as.numeric(substr(what$Block,1,2)))
    LONG=100+(as.numeric(substr(what$Block,3,4)))
    
    plotMap(worldLLhigh, xlim=XAXIS,ylim=YAXIS,
            plt = NULL,col="grey80",tck = 0.025, tckMinor = 0.0125,
            xlab="",ylab="",axes=F)
    
    points(LONG+0.5,LAT-0.5,cex=(what$Pred.mean/mean(what$Pred.mean))*2.5,pch=19,col="grey20")
    box()
    axis(2,seq(YAXIS[1],YAXIS[2],1),F,tck=-0.02)
    axis(2,seq(YAXIS[1],YAXIS[2],4),F,tck=-0.04)
    axis(1,seq(XAXIS[1],XAXIS[2],1),F,tck=-0.02)
    axis(1,seq(XAXIS[1],XAXIS[2],4),F,tck=-0.04)
    
    if(LEG=="YES")legend("topright",c("5 kg","10 kg","15 kg"),cex=1.25,pch=19,
                         pt.cex=(c(5,10,15)/mean(what$Pred.mean))*2.5,bty="n",col="grey20")
    
  }
  AXIS1=function()axis(2,seq(YAXIS[1],YAXIS[2],4),-seq(YAXIS[1],YAXIS[2],4),cex.axis=1.25,las=2,tck=-0.04)
  AXIS2=function()axis(1,seq(XAXIS[1],XAXIS[2],4),seq(XAXIS[1],XAXIS[2],4),cex.axis=1.25,tck=-0.04)
  
  
  XAXIS=c(112,122)
  YAXIS=c(-36,-26)
  
  tiff(file="paper/pred Mean weight by block.tiff",width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfcol=c(2,2),mar=c(2,3,2,.1),oma=c(2,1,.1,.1),las=1,mgp=c(.1,.8,0))
  
  fun.plot.map(Pred.gum,"YES")
  mtext("Gummy shark",3,.5,cex=1.5)
  AXIS1()
  AXIS2()
  
  fun.plot.map(Pred.dus,"NO")
  mtext("Dusky shark",3,.5,cex=1.5)
  AXIS1()
  AXIS2()
  
  fun.plot.map(Pred.whis,"NO")
  mtext("Whiskery shark",3,.5,cex=1.5)
  AXIS1()
  AXIS2()
  
  fun.plot.map(Pred.san,"NO")
  mtext("Sandbar shark",3,.5,cex=1.5)
  AXIS1()
  AXIS2()
  
  #text(119,-29,"Western",cex=1.75)
  #text(119,-30,"Australia",cex=1.75)
  
  mtext("Latitude (?S)",side=2,line=-0.8,font=1,las=0,cex=1.65,outer=T)
  mtext("Longitude (?E)",side=1,line=0.4,font=1,las=0,cex=1.65,outer=T)
  
  par(fig=c(.5,.92,.5,.92), new = T,mgp=c(.1,.4,0))
  plotMap(worldLLhigh, xlim=c(110,155),ylim=c(-45,-10),plt = c(.6, 1, .6, 1),
          col="grey80",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  box()
  polygon(x=c(XAXIS,rev(XAXIS)),y=c(YAXIS[2],YAXIS[2],YAXIS[1],YAXIS[1]),lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
  text(134,-23.5,("Australia"),col="black", cex=1.35)
  
  dev.off()
  
  
  
  #Plot predictions
  fn.FL.plot=function(what)
  {
    N=length(what$Block)
    
    MEAN=mean(what$Pred.mean,na.rm=T)
    std.mean=what$Pred.mean/MEAN
    UP1=(what$Pred.mean/MEAN)+(what$Pred.SE/MEAN)*2
    LOW1=(what$Pred.mean/MEAN)-(what$Pred.SE/MEAN)*2
    
    #MAXY=max(UP1)
    #MINY=min(LOW1)
    MAXY=2.5
    MINY=0
    plot(1:N,std.mean,ylab="",xlab="",xaxt='n',cex.axis=1.25,
         ylim=c(MINY,MAXY),pch=19,cex=1.25)
    segments(1:N,LOW1,1:N,UP1,lwd=1.5)
    axis(1,1:N,F)
    axis(1,seq(1,N,2),what$Block[seq(1,N,2)],tck=-0.03,cex=1.25)
  }
  
  
  par(mfcol=c(2,2),mar=c(1,4,2,.1),oma=c(3,1,.1,.1),las=1,mgp=c(.1,.7,0))
  
  fn.FL.plot(Pred.gum)
  mtext("Gummy shark",3,cex=1.25)
  
  fn.FL.plot(Pred.dus)
  mtext("Dusky shark",3,cex=1.25)
  
  fn.FL.plot(Pred.whis)
  mtext("Whiskery shark",3,cex=1.25)
  
  fn.FL.plot(Pred.san)
  mtext("Sandbar shark",3,cex=1.25)
  
  mtext("Fishing block",1,line=1.5,cex=1.5,outer=T)
  mtext("Normalised mean weight",2,line=-1,cex=1.5,outer=T,las=3)
  
  
  # #LSMEAS approach   (package lmerTest)
  # fn.lsmeans=function(model,var) LSMEAS=lsmeans(model, test.effs=var)
  # 
  # Lsmeans.whis=fn.lsmeans(Whis.fit$model,"finyear")
  # Lsmeans.gum=fn.lsmeans(Gum.fit$model,"finyear")
  # Lsmeans.dus=fn.lsmeans(Dus.fit$model,"finyear")
  # Lsmeans.san=fn.lsmeans(San.fit$model,"finyear")
  # 
  # 
  # plot(Lsmeans.whis)
  # plot(Lsmeans.gum)
  # plot(Lsmeans.dus)
  # plot(Lsmeans.san)
  
  
}


#FL CHANGES BY SPACE AND TIME
if(Run=="First")
{
  depth.range=seq(0,120,20)
  depthMax=seq(10,150,10)
  fn.FL=function(dat,Model.type)
  {
    SPEC=unique(dat$species)
    dat$finyear=as.factor(dat$finyear)
    dat$month=as.factor(dat$month)
    dat$blockx=as.factor(dat$blockx)
    dat$vessel=as.factor(dat$vessel)
    dat$logFL=log(dat$FL)
    
    
    #   dat$Depth.bin=cut(dat$depthMax,depth.range)
    #   depths=levels(dat$Depth.bin)
    #   par(mfcol=c(3,2))
    #   for( i in 1:length(depths))
    #   {
    #     a=subset(dat,Depth.bin==depths[i])
    #     hist(a$FL,breaks=50,xlim=c(min(dat$FL),max(dat$FL)),main=paste("depth=",depths[i]))
    #   }
    #   mtext(SPEC,3,outer=T,line=-2)
    
    
    if(Model.type=="fixed")
    {
      model=glm(logFL~finyear+blockx+month+depthMax+vessel,family=gaussian, data = dat)
      null.model=glm(logFL~1,family=gaussian, data = dat)
    }
    
    if(Model.type=="Mixed")
    {
      model=lmer(logFL~finyear+blockx+month+depthMax+(1 |vessel), data = dat)
      null.model=lmer(logFL~1+(1 |vessel), data = dat)
    }
    
    
    #Anova
    ANOVA=fn.anova(model,Model.type)
    
    #Overall deviance
    Dev.exp= 100*abs((deviance(null.model) - deviance(model)) / deviance(null.model))
    
    #Terms deviance
    terms.dev.exp=NA
    if(Model.type=="fixed")
    {
      terms.dev.exp=fn.dev.exp.term(ANOVA,model)
      names(terms.dev.exp)=row.names(ANOVA)[2:nrow(ANOVA)]
    }
    
    
    #Predictions
    #1. create data to predict
    Mode=table(dat$finyear)  
    finyear=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$finyear))
    
    Mode=table(dat$vessel)  
    vessel=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$vessel))
    
    Mode=table(dat$month)  
    month=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$month))
    
    Mode=table(dat$blockx)  
    blockx=factor(names(Mode[which(Mode==max(Mode))]),levels(dat$blockx))
    
    
    #2. predict data
    #2.2.1 depth  
    Depth.round=round(dat$depthMax/10)*10
    N.depth=c(unlist(table(Depth.round)))
    N.depth=N.depth[2:length(N.depth)]
    Depth.range=as.numeric(names(N.depth))
    NEWdata=data.frame(finyear=finyear,month=month,vessel=vessel,
                       blockx=blockx,depthMax=Depth.range)
    
    
    if(Model.type=="fixed")
    {
      PREDS=predict(model,newdata=NEWdata,type='response',se.fit=T)
      biasCorr <- model$deviance/model$df.residual/2   
      Pred.mean=PREDS$fit+biasCorr
      Pred.SE=PREDS$se.fit+biasCorr
    }
    
    if(Model.type=="Mixed") 
    {
      NEWdata$logFL=0
      PRED.MIX=fn.pred.mixed(model,NEWdata,N.depth)
      Pred.mean=PRED.MIX$Pred.mean
      Pred.SE=PRED.MIX$Pred.SE
    }
    
    
    #2.2.2 block
    NEWdata=data.frame(finyear=finyear,month=month,vessel=vessel,
                       blockx=levels(dat$blockx),depthMax=mean(dat$depthMax))
    N.blk=c(unlist(table(dat$blockx)))
    
    if(Model.type=="fixed")
    {
      PREDS=predict(model,newdata=NEWdata,type='response',se.fit=T)
      Pred.mean.BLK=PREDS$fit+biasCorr
      Pred.SE.BLK=PREDS$se.fit+biasCorr
    }
    
    if(Model.type=="Mixed") 
    {
      NEWdata$logFL=0
      PRED.MIX=fn.pred.mixed(model,NEWdata,N.blk)
      Pred.mean.BLK=PRED.MIX$Pred.mean
      Pred.SE.BLK=PRED.MIX$Pred.SE
    }
    
    return(list(Dev.exp=Dev.exp,Anova=ANOVA,terms.dev.exp=terms.dev.exp,Depth.range=Depth.range, 
                pred.depth=Pred.mean,pred.depth.SE=Pred.SE,Pred.mean.BLK=Pred.mean.BLK,
                Pred.SE.BLK=Pred.SE.BLK,BLKS=levels(dat$blockx)))
  }
  
  FL.whi=fn.FL(subset(Survey.whi,species==17003 & !is.na(depthMax) & FL<Max.FL.w &!finyear%in%c("1999-00","2003-04")),"Mixed")
  FL.gum=fn.FL(subset(Survey.gum,species==17001 & !is.na(depthMax) & FL<Max.FL.g &!finyear%in%c("1999-00","2003-04")),"Mixed")
  FL.dus=fn.FL(subset(Survey.dus,species==18003 & !is.na(depthMax) & FL<Max.FL.d &!finyear%in%c("1999-00","2003-04")),"Mixed")
  FL.san=fn.FL(subset(Survey.san,species==18007 & !is.na(depthMax) & FL<Max.FL.s &!finyear%in%c("1999-00","2003-04")),"Mixed")
  
  
  
  #Plot example of difference in size despite using same gear
  
  #Spatial display
  fun.plot.map=function(what,LEG)
  {
    
    LAT=-(as.numeric(substr(what$BLKS,1,2)))
    LONG=100+(as.numeric(substr(what$BLKS,3,4)))
    
    plotMap(worldLLhigh, xlim=XAXIS,ylim=YAXIS,
            plt = NULL,col="grey70",tck = 0.025, tckMinor = 0.0125,
            xlab="",ylab="",axes=F)
    
    points(LONG+0.5,LAT-0.5,cex=(exp(what$Pred.mean.BLK)/mean(exp(what$Pred.mean.BLK)))*3,pch=19,col="grey20")
    box()
    axis(2,seq(YAXIS[1],YAXIS[2],1),F,tck=-0.02)
    axis(2,seq(YAXIS[1],YAXIS[2],4),F,tck=-0.04)
    axis(1,seq(XAXIS[1],XAXIS[2],1),F,tck=-0.02)
    axis(1,seq(XAXIS[1],XAXIS[2],4),F,tck=-0.04)
    
    if(LEG=="YES")legend("right",c("80 cm","100 cm","120 cm"),cex=1.25,pch=19,pt.cex=(c(80,100,120)/mean(exp(what$Pred.mean.BLK)))*3,bty="n",col="grey20")
    
  }
  
  fn.FL.plot.depth=function(what)
  {
    MAXY=1.2
    MINY=0.8
    
    what$pred.depth=exp(what$pred.depth)
    what$pred.depth.SE=exp(what$pred.depth.SE)
    MEAN=mean(what$pred.depth)
    UP1=(what$pred.depth/MEAN)+(what$pred.depth.SE/MEAN)*2
    LOW1=(what$pred.depth/MEAN)-(what$pred.depth.SE/MEAN)*2
    
    plot(what$Depth.range,what$pred.depth/MEAN,ylab="",xlab="",xaxt='n',cex.axis=1.25,
         ylim=c(MINY,MAXY),pch=19,cex=1.25,col="transparent")
    CI.fun(what$Depth.range,UP1,LOW1,"grey75","transparent")
    lines(what$Depth.range,what$pred.depth/MEAN,lwd=2)
    
  }
  
  XAXIS=c(112,129)
  YAXIS=c(-36,-26)
  
  tiff(file="paper/pred.FL.block_depth.tiff",width = 2400, height = 2000,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(2,2),mar=c(2,3,2,.1),oma=c(2,1,.1,.1),las=1,mgp=c(.1,.8,0))
  
  #Blocks
  fun.plot.map(FL.dus,"YES")
  mtext("Dusky shark",3,.5,cex=1.5)
  axis(2,seq(YAXIS[1],YAXIS[2],4),-seq(YAXIS[1],YAXIS[2],4),cex.axis=1.25,las=2,tck=-0.04)
  axis(1,seq(XAXIS[1],XAXIS[2],4),seq(XAXIS[1],XAXIS[2],4),cex.axis=1.25,tck=-0.04)
  mtext("Latitude (?S)",side=2,line=2.25,font=1,las=0,cex=1.65)
  
  fun.plot.map(FL.san,"NO")
  mtext("Sandbar shark",3,.5,cex=1.5)
  axis(1,seq(XAXIS[1],XAXIS[2],4),seq(XAXIS[1],XAXIS[2],4),cex.axis=1.25,tck=-0.04)
  axis(2,seq(YAXIS[1],YAXIS[2],4),-seq(YAXIS[1],YAXIS[2],4),cex.axis=1.25,las=2,tck=-0.04)
  mtext("Longitude (?E)                                             ",side=1,
        line=2.25,font=1,las=0,cex=1.65)
  
  #Depth
  fn.FL.plot.depth(FL.dus)
  mtext("Normalised mean FL",side=2,line=2.25,font=1,las=0,cex=1.65)
  axis(1,seq(20,140,20),seq(20,140,20),cex.axis=1.25,tck=-0.04)
  axis(1,depthMax,F,tck=-0.02)
  
  fn.FL.plot.depth(FL.san)
  axis(1,seq(20,140,20),seq(20,140,20),cex.axis=1.25,tck=-0.04)
  axis(1,depthMax,F,tck=-0.02)
  mtext(" Depth (m)",1,line=0.25,cex=1.65,outer=T)
  
  #Add inset
  par(fig=c(.5,.92,.5,.92), new = T,mgp=c(.1,.4,0))
  plotMap(worldLLhigh, xlim=c(110,155),ylim=c(-45,-10),plt = c(.6, .99, .6, .99),
          col="grey90",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  box()
  polygon(x=c(XAXIS,rev(XAXIS)),y=c(YAXIS[2],YAXIS[2],YAXIS[1],YAXIS[1]),lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
  text(134,-23.5,("Australia"),col="black", cex=1.35)
  dev.off()
  
  
  
  
  fn.FL.plot=function(what)
  {
    N=length(what$BLKS)
    MAXY=max(what$Pred.mean.BLK+2*(what$Pred.SE.BLK))
    MINY=min(what$Pred.mean.BLK-2*(what$Pred.SE.BLK))
    plot(1:N,what$Pred.mean.BLK,ylab="",xlab="",xaxt='n',cex.axis=1.25,
         ylim=c(MINY,MAXY),pch=19,cex=1.25)
    segments(1:N,what$Pred.mean.BLK-2*(what$Pred.SE.BLK),
             1:N,what$Pred.mean.BLK+2*(what$Pred.SE.BLK),lwd=1.5)
    axis(1,1:N,F)
    axis(1,seq(1,N,2),what$BLKS[seq(1,N,2)],tck=-0.03,cex=1.25)
  }
  
  
  par(mfcol=c(2,2),mar=c(1,4,2,.1),oma=c(3,1,.1,.1),las=1,mgp=c(.1,.7,0))
  
  fn.FL.plot(FL.gum)
  mtext("Gummy shark",3,cex=1.25)
  
  fn.FL.plot(FL.dus)
  mtext("Dusky shark",3,cex=1.25)
  
  fn.FL.plot(FL.whi)
  mtext("Whiskery shark",3,cex=1.25)
  
  fn.FL.plot(FL.san)
  mtext("Sandbar shark",3,cex=1.25)
  
  mtext("Fishing block",1,line=1.5,cex=1.5,outer=T)
  mtext("Mean fork length (cm)",2,line=-1,cex=1.5,outer=T,las=3)
  
}



# EXPORT QUANTITIES OF INTEREST FOR POP DYN MODEL -------------------------
setwd(handl_OneDrive('Analyses\\Data_outs'))

#normalised
write.csv(Pred.gum%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%select(Finyear,mean,CV),
          "Gummy shark/Gummy shark.annual.mean.size_relative.csv",row.names = F)
write.csv(Pred.dus%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%select(Finyear,mean,CV),
          "Dusky shark/Dusky shark.annual.mean.size_relative.csv",row.names = F)
write.csv(Pred.whis%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%select(Finyear,mean,CV),
          "Whiskery shark/Whiskery shark.annual.mean.size_relative.csv",row.names = F)
write.csv(Pred.san%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%select(Finyear,mean,CV),
          "Sandbar shark/Sandbar shark.annual.mean.size_relative.csv",row.names = F)

  
write.csv(Pred.smh%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%select(Finyear,mean,CV),
          "Smooth hammerhead/Smooth hammerhead.annual.mean.size_relative.csv",row.names = F)
write.csv(Pred.spi%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%select(Finyear,mean,CV),
          "Spinner shark/Spinner shark.annual.mean.size_relative.csv",row.names = F)
write.csv(Pred.tig%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%select(Finyear,mean,CV),
          "Tiger shark/Tiger shark.annual.mean.size_relative.csv",row.names = F)
write.csv(Pred.cop%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%select(Finyear,mean,CV),
          "Copper shark/Copper shark.annual.mean.size_relative.csv",row.names = F)


#by zone
fn.out=function(LisT,NM)
{
  for(l in 1:length(LisT))
  {
   write.csv(LisT[[l]]%>%mutate(mean=Pred.mean/mean(Pred.mean))%>%
                select(Finyear,mean,CV),
              paste(getwd(),'/',NM,'/',NM,".annual.mean.size_relative_",names(LisT)[l],".csv",sep=""),row.names = F)
  }
}

fn.out(list(west=Pred.whis_west,zone1=Pred.whis_zn1,zone2=Pred.whis_zn2),
       NM="Whiskery shark")
fn.out(list(zone2=Pred.gum_zn2),NM="Gummy shark")
fn.out(list(west=Pred.dus_west,zone1=Pred.dus_zn1,zone2=Pred.dus_zn2),
       NM="Dusky shark")
fn.out(list(west=Pred.san_west,zone1=Pred.san_zn1),NM="Sandbar shark")




