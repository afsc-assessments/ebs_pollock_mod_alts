rm(list=ls())
R
ls()
source("R/prelims.R")
#-------------------------------------------------------------------------------
# Visual compare runs
#-------------------------------------------------------------------------------
source("R/compareRuns.r")

# Read in the output of the assessment
mod14.1 <- readList("mod14.1/For_R.rep")
mod0.0 <- readList("mod0.0/For_R.rep")
names(mod14.1)
# Read in model results
i=0
for (i in 1:1) {
  rn=paste0("mod14.",i,"/For_R.rep"); mn=paste0("mod14.",i)
  assign(mn,readList(rn))
  print(rn)
}
for (i in 0:1) {
  rn=paste0("mod16.",i,"/For_R.rep"); mn=paste0("mod16.",i)
  assign(mn,readList(rn))
  print(rn)
}
lstOuts   <- list( Model_14.1= mod14.1, Model_16.0= mod16.0)
lstOuts   <- list( "2015_Model"= mod0.0, Model_16.0= mod16.0)
#lstOuts   <- list( Model_14.1= mod14.1, Model_16.0= mod16.0, Model_16.1= mod16.1, Model_16.2= mod16.2)
#lstOuts   <- list( Model_14.1= mod14.1, Model_16.0= mod16.0, Model_16.1= mod16.1)
#lstOuts   <- list( Model_0.0=mod0.0, Model_0.1=mod0.1 , Model_0.2= mod0.2, Model_0.3= mod3)

tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
tab
write.csv(tab,file=file.path(inputPath,"LikelihoodTable2015.csv"),row.names=F)
system(paste("open ",inputPath,"/LikelihoodTable2015.csv",sep=""))
#----parameter table-----------------
rt <- get_params(lstOuts)
rt

#-------------------------------------------------------------------------------
# Fit to survey data
IndexFit(mod14.1,yf=1990,yl=2016,f=1,main="Model 14.1",ylab="Survey biomass (t)")
IndexFit(mod16.0,yf=1990,yl=2016,f=1,main="Model 16.0",ylab="Survey biomass (t)")
AgeFitsSrv(mod16.0,rec_age=1,case_label="Model 16.0")

#-------------------------------------------------------------------------------
# Fit to fishery data
AgeFits(mod16.0,rec_age=1,case_label="Model 16.0 ")
mod16.0$Obs_Survey_1[,1:5] 
dt <- data.table(mod14.1$Obs_Survey_1[,1:4],mod16.0$Obs_Survey_1[,3],mod16.1$Obs_Survey_1[,3])
#dt <- data.frame(Model_0.1=mod0.1$Obs_Survey_1[,3], Model_0.1=mod0.2$Obs_Survey_1[,3], Model_0.3=mod0.3$Obs_Survey_1[,3])
names(dt) <- c("Year","Observed","Model_14.1","sd","Model_16.0","Model_16.1")
dt
dt <- data.table(dt)[Year>1986,.(Year, Observed, Model_14.1,Model_16.0, Model_16.1, lb=Observed-1.96*sd,ub=Observed+1.96*sd),]
ggplot(dt,aes(x=Year,y=Model_14.1) ) + geom_point(data=dt[variable=="Observed"],size=3) + geom_line(data=dt[variable!="Observed"]) + ylim(c(0,1800000)) + labs(x="Year",y="Survey biomass") + mytheme #+
dtp <-melt(dt,id="Year")
dt
,aes(x=Age,y=value,colour=variable)) 
dt[variable!="ub"&variable!="lb"] %>% 
ggplot(aes(x=Year,y=value,col=variable)) + geom_point(data=dt[variable=="Observed"],size=3) + geom_line(data=dt[variable!="Observed"]) + ylim(c(0,1800000)) + labs(x="Year",y="Survey biomass") + mytheme #+
            geom_line(data=)
#-------------------------------------------------------------------------------
#Selectivity
df <- data.frame(mod14.1$sel_fsh_1[,2:13] ); names(df) <- c("yr",1:11)
df <- data.frame(mod16.0$sel_fsh_1[,2:13] ); names(df) <- c("yr",1:11)
sdf <- gather(df,age,sel,2:12) %>% filter(yr>1990) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
p1 <- ggplot(sdf, aes(x=age,y=sel/2+yr,group=yr)) +                  geom_ribbon(aes(ymin=yr,ymax=sel/1.9+yr),fill="tan",col="grey60",alpha=.3)  + ylab("")                        + xlab("Age") + guides(fill=FALSE,alpha=FALSE,col=FALSE) + mytheme
p1 <- p1 + scale_y_continuous(breaks=seq(1990,2015,2)) + scale_x_continuous(breaks = c(1,3,5,7,9,11))
p1

df <- data.frame(mod3$sel_fsh_1[,2:13] ); names(df) <- c("yr",1:11)
sdf <- gather(df,age,sel,2:12) %>% filter(yr>1990) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
p2 <- ggplot(sdf, aes(x=age,y=sel/2+yr,group=yr)) +                  geom_ribbon(aes(ymin=yr,ymax=sel/1.9+yr),fill="tan",col="grey60",alpha=.3)  + ylab("")                        + xlab("Age") + guides(fill=FALSE,alpha=FALSE,col=FALSE) + mytheme
p2 <- p2 + scale_y_continuous(breaks=seq(1990,2015,2)) + scale_x_continuous(breaks = c(1,3,5,7,9,11))
grid.arrange(p1,p2,ncol=2)

df
p1
#--Projections---------
pfn <- "5yr"
pdt <- data.table(read.table(paste0("proj/",pfn,"_out/bigfile.out"),header=TRUE)) 
pdt1 <- data.table(read.table(paste0("proj/",pfn,"_out/catch_fixed_2017_18.out"),header=TRUE)) 
pdt1[Yr>2016&Yr<2019,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(SSB,.2) ,ub=quantile(SSB,.8) ),.(Yr,Alternative)]
pfn <- "10yr"
pdt2 <- data.table(read.table(paste0("proj/",pfn,"_out/bigfile.out"),header=TRUE)) 
pdt$Alternative <- as.factor(pdt$Alternative)
pdt2$Alternative <- as.factor(pdt2$Alternative)
pdt$config <-"5yr"
pdt2$config <-"10yr"

pt <- pdt[Yr>2016,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(SSB,.2) ,ub=quantile(SSB,.8) ),.(Yr,Alternative)]
pt
ggplot(pt,aes(x=Yr,y=SSB,fill=Alternative)) + geom_line() + mytheme + ylim(c(0,400)) + geom_ribbon(aes(ymin=lb,ymax=ub,fill=Alternative),alpha=0.25) + labs(y="Spawning biomass (kt)",x="Year") + scale_x_continuous(breaks=seq(2015,2032,2)) 
c1 <- ggplot(pt,aes(x=Yr,y=Catch,color=Alternative,size=1.)) + geom_line(size=1.5) + mytheme + labs(y="Catch (kt)",x="Year") + scale_x_continuous(breaks=seq(2015,2032,2)) 
c1 <- c1 + geom_line(aes(x=Yr,y=ABC),size=1)
#c1 <- c1 + geom_line(data=pt[as.numeric(Alternative)==2,.(Yr,ABC)],aes(x=Yr,y=ABC))
c1
pt[as.numeric(Alternative)==2,.(Yr,ABC)]
pt <- pdt[Yr>2016,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL)),.(Yr,Alternative)] 
pt
ggplot(pt,aes(x=Yr,y=OFL,color=Alternative)) + geom_line() + mytheme
pdt
pdx <-rbind(pdt2,pdt)
setkey(pdx,Yr,Alternative)
pt <- pdx[.(Yr>2016,(Alternative)==1),.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB)  ),.(Yr,config)]
pt <- pdx[Yr>2016&Alternative==1,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(Catch,.1) ,ub=quantile(Catch,.9)  ),.(Yr,config)]
pt

ggplot(pt,aes(x=Yr,y=Catch,color=config,shape=config)) + geom_line() + geom_point() + mytheme + geom_ribbon(aes(ymin=lb,ymax=ub,fill=config),alpha=0.25) + labs(x="Year")+ scale_x_continuous(breaks=seq(2015,2032,2))  + ylim(c(0,130))

AgeFits(mod4,rec_age=1,case_label="Model 4 assessment")
IndexFit(mod2,yf=1980,yl=2015,f=1,main="Model 2",ylab="Survey biomass (t)")
AgeFits(mod1,rec_age=1,case_label="Model 1 assessment")
AgeFits(mod2,rec_age=1,case_label="Model 2 assessment")
AgeFitsSrv(mod1,rec_age=1,case_label="2016 assessment")
AgeFitsSrv(mod16.0,rec_age=1,case_label="2016 assessment")

#--Selectivity==============================================
df <- data.frame(Model="Model 16.0",mod16.0$sel_fsh_1[,2:13] ); 
for (i in 2:7) df <- rbind(df, data.frame(Model=paste0("Model ",i),lstOuts[[i]]$sel_fsh_1[,2:13] ))
names(df) <- c("Model","yr",1:11)
sdf <- gather(df,age,sel,3:13)
tbl_df(sdf)
sdf$age <- as.numeric(sdf$age)
sdf <- data.table(sdf)
sdf <- sdf[yr>2010, (Selectivity=mean(sel)),.(Model,age)]
ggplot(sdf ,aes(x=age,y=V1,colour=Model),size=1.2) + geom_line(size=2) +ylab("Fishery selectivity") +xlab("Age")+ mytheme #+ facet_grid(yr~.)
#-------------------------------------------------------------------------------u

p1 <- dplyr::filter(sdf,yr>1979) %>% arrange(yr,age) %>%ggplot(aes(x=age,y=sel/2+yr,group=yr)) + geom_ribbon(aes(ymin=yr,ymax=sel/1.9+yr),fill="tan",col="grey60",alpha=.3)  + ylab("Selectivity by year") + xlab("Age") + guides(fill=FALSE,alpha=FALSE,col=FALSE) + mytheme
p1
dev.off()
compareTime(lstOuts,"SSB",SD=T,Sum=NULL,legendPos="top",startYear=1980)
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,startYear=2000,ylim=c(0,3.4e5))
compareTime(lstOuts,"R",SD=F,ylim=c(0,1500),startYear=2000)
compareTime(lstOuts,"R",SD=T,cex.lab=2)
compareTime(lstOuts,"TotBiom",SD=T)
compareMatrix(lstOuts,"TotF",SD=TRUE,Sum=NULL,YrInd=mod1$Yr,Apply=mean,legendPos="right")
#-------------------------------------------------------------------------------

p1 <- dplyr::filter(sdf,yr>1979) %>% arrange(yr,age) %>%ggplot(aes(x=age,y=sel/2+yr,group=yr)) + geom_ribbon(aes(ymin=yr,ymax=sel/1.9+yr),fill="tan",col="grey60",alpha=.3)  + ylab("Selectivity by year") + xlab("Age") + guides(fill=FALSE,alpha=FALSE,col=FALSE) + mytheme
p1
p2 <- dplyr::filter(sdf,yr>2002) %>% ggplot(aes(x=age,y=sel/2+yr,group=yr)) +                  geom_ribbon(aes(ymin=yr,ymax=sel/1.9+yr),fill="tan",col="grey60",alpha=.3)  + ylab("")                        + xlab("Age") + guides(fill=FALSE,alpha=FALSE,col=FALSE) + mytheme
grid.arrange(p1,p2,ncol=2)

p2 <- dplyr::filter(sdf,yr>1977) %>% ggplot(aes(x=age,y=sel/2+yr,group=yr)) + geom_ribbon(aes(ymin=yr,ymax=sel/1.6+yr,fill="salmon",col="grey",alpha=.2))  + ylab("")                        + xlab("Age") + guides(fill=FALSE,alpha=FALSE,col=FALSE) + mytheme

#Survey selectivity==============================================
srv_sel1 <- mod14.1$sel_ind_1[1,3:13]
srv_sel1 <- data.frame(Age=1:11,Selectivity=srv_sel1/max(srv_sel1))
srv_sel1$Model <-"Model 14.1"

srv_sel2 <- mod16.0$sel_ind_1[1,3:13]
srv_sel2<- data.frame(Age=1:11,Selectivity=srv_sel2/max(srv_sel2))
srv_sel2$Model <-"Model 16.0"
srv_sel <- rbind(srv_sel1,srv_sel2)

mytheme = mytheme + theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1))
ggplot(srv_sel) + geom_line(aes(x=Age,y=Selectivity,color=Model),size=2) +mytheme


#SSB ============================================================
df  <- data.table(Model = "16.0", mod16.0$SSB )
df <- rbind(df, data.table(Model="14.1",mod14.1$SSB))
for (i in 2:4) df <- rbind(df, data.table(Model=paste0("Model ",i),lstOuts[[i]]$SSB))
names(df) <- c("Model","yr","SSB","SE","lb","ub")
bdf <- filter(df,yr>1980,yr<=2016) %>% arrange(yr)
bdf

# try over same time range...
p1 <- ggplot(data=bdf,aes(x=yr,y=SSB,alpha=0.3,color=Model)) + scale_y_continuous(limits=c(0,600000)) + ylab("Spawning biomass") + xlab("Year") +  mytheme + geom_line(data=bdf,aes(x=yr,y=SSB,type=Model)) +
           geom_ribbon(data=bdf ,aes(x=yr,y=SSB,ymin=lb,ymax=ub,fill=Model),alpha=.3)  + guides(alpha=FALSE,col=FALSE) 
p1
#Recruits==========================================================
mytheme = mytheme + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
#Recruits==========================================================
mytheme = mytheme + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0))
rdf <- data.table(mod16.0$R)
rdf <- cbind(rdf,"Model 16.0")
#rdf <- data.table(mod14.1$R)
rdf0 <- data.table(mod0.0$R)
rdf0<- cbind(rdf0,"2015 Model")
names(rdf) <- c("yr","R","se","lb","ub","case")
names(rdf0) <- c("yr","R","se","lb","ub","case")
rdf <- rbind(rdf,rdf0)

#rdf  <- rdf[yr>1990&yr<2017,.(yr,R=R/mean(R),lb=lb/mean(R),ub=ub/mean(R),case)]
rdf  <- rdf[yr>1990&yr<2017,.(yr,R=R,lb=lb,ub=ub,case)]
mnR <- mean(rdf$R)
mnR
dodge <- position_dodge(width=0.8)
ggplot(rdf,aes(x=yr-1,y=R,fill=case)) + xlab("Year class") + ylab("Age 1 recruits ") + #ylim(c(0,18000)) +
       geom_bar(width=0.75,position="dodge",stat="identity") + #xlim(c(1997.5,2014.5)) +
       geom_errorbar(aes(ymin=lb,ymax=ub),width=.3,colour="blue",position=dodge) + mytheme + geom_hline(aes(yintercept=mnR))


rdf <- data.frame(mod1$R)
rdf <- cbind(rdf,"2015")
names(rdf) <- c("yr","R","se","lb","ub","case")
rdf  <- filter(rdf , yr>1970,yr<2016)
mnR <- mean(rdf$R)
mnR
rdf0 <- data.frame(mod6$R)
rdf0 <- cbind(rdf0,"mod6")
names(rdf0) <- c("yr","R","se","lb","ub","case")
rdf0 <- filter(rdf0, yr>1970,yr<2016)
rdf <- rbind(rdf0,rdf)
tbl_df(rdf)
dodfact=0.90
ggplot(rdf,aes(x=as.factor(yr-1),y=R)) + xlab("Year class") + ylab("Age 1 recruits (10e6)") + 
       geom_bar(width=0.85,position=position_dodge(width = dodfact),stat="identity",aes(fill=case)) + 
       geom_errorbar(aes(ymin=lb,ymax=ub),width=.4,col="blue",position=position_dodge(width=dodfact)) + mytheme + geom_hline(aes(yintercept=mnR),linetype="dashed")

  #scale_y_continuous(limits=c(0,1800))+
     ggtitle("2014") + mytheme
     ggtitle("2015") + mytheme

#=====================
# Retro
# get all retrospectives
#=====================
# Read in retro results
for (i in 0:15) {
  rn=paste0("mod16.0/retro/retro/r_",i,".rep")
  #rn=paste0("mod14.1/retro/retro/r_",i,".rep")
  mn=paste0("retro",i)
  assign(mn,readList(rn))
  print(rn)
}
retouts <- list()
retouts <- list( R0=retro0, R1= retro1, R2= retro2, R3= retro3, R4= retro4 , R5= retro5 , R6= retro6 , R7= retro7 , R8= retro8 , R9= retro9 , R10= retro10 )  
names(retouts$R0)
tab = list(R0=retro0$SSB)
tab
getwd()
for (i in 1:10) { tab=cbind(tab,retouts[[i]]$SSB)}

#pdf(paste(Figdir,"Retro_Mods.pdf",sep=""),width=9, height=6)
names(mod1$SSB)
rn = "mod1"
dim(tdf)
library(grid)
i=0
mytheme = mytheme + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0))
p1 <- ggplot() + scale_y_continuous(limits=c(0,510000)) + ylab("Spawning biomass") + xlab("Year") +  mytheme + geom_line(data=bdf,aes(x=yr,y=SSB),size=4) +
           geom_ribbon(data=bdf ,aes(x=yr,y=SSB,ymin=lb,ymax=ub),fill="tan",col="grey",alpha=.6)  + guides(fill=FALSE,alpha=FALSE,col=FALSE) + scale_x_continuous(breaks=seq(1978,2016,2)) 
for (i in 1:10) {
  rn=paste("retro",i,sep="");
  tdf <- data.frame(get(rn)$SSB); names(tdf) <- c("yr","SSB","SE","lb","ub"); tdf <- filter(tdf,yr>1977)
  p1 <- p1 + geom_line(data=tdf,aes(x=yr,y=SSB),col=i,linetype=i,size=1.25)
  #p1 <- p1 + geom_segment(data=tdf,aes(x=yr,xend=yr,yend=SSB,y=SSB),arrow=arrow(angle=90,length=unit(.2,"cm")),size=2,col=i)
  tdf <- tdf[dim(tdf)[1],]
  p1 <- p1 + geom_point(data=tdf,aes(x=yr,y=SSB),size=4,col=i)
  #p1 <- p1 + geom_point(get(rn)$SSB[lr,1],get(rn)$SSB[lr,2],pch=19,col=i)
}
p1
# If you want to post lines when new survey data arrive...
# p1 <- p1 +  geom_vline(xintercept=c(2014,2012,2010,2008),linetype="dashed",size=1,col="grey")
# Color blind palette
# cb_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# p1 <-    p1 + scale_fill_manual(values=cb_palette)
df <- data.frame(mod16.0$SSB )
df <- data.frame(mod14.1$SSB )
names(df) <- c("yr","SSB","SE","lb","ub")
bdf <- filter(df,yr>1977,yr<=2016)
bdf2<- filter(df,yr>1977)
bdf
bdft <- bdf
for (i in 1:10) bdft <- cbind(bdft,rep(NA,39))
names(bdft)[6:15] <- paste0("SSB_",2015:2006)
for (i in 1:10) bdft[1:(39-i),i+5] <- get(paste0("retro",i))$SSB[14:(52-i),2]
bdft
for (i in 1:10) bdft[,i+5] <- bdft[,i+5]/bdf$SSB
seq(1978,2016,3)
p2 <- ggplot(bdft,aes(x=yr,y=SSB)) + scale_x_continuous(breaks=seq(1978,2016,2)) + scale_y_continuous(limits=c(.25,1.75)) + ylab("Relative spawning biomass") + xlab("Year") +  mytheme  
for (i in 1:10) {
  tdf <- data.frame(cbind(bdft[1],SSB=bdft[5+i]))[1:(40-i),]
  names(tdf) <- c("yr", "SSB")
  p2 <- p2 + geom_line(data=tdf, aes(x=yr,y=SSB),col=i,size=1.5)
  tdf <- tdf[dim(tdf)[1],]
  p2 <- p2 + geom_point(data=tdf,aes(x=yr,y=SSB),size=4,col=i)
}    
#p2 <-p2 +geom_hline(aes(yintercept=1),size=3,linetype="dotted")
p2 <-p2 +geom_hline(aes(yintercept=1),size=1,col="grey")
p2
grid.arrange(p1,p2,nrow=2)
# Mohn's rho
rc = retro0$SSB[,2]
ntmp=0
rho=0
for (i in 1:15) {
  dtmp=get(paste0("retro",i))$SSB
  lr=length(dtmp[,1])
  ntmp= ntmp+(dtmp[lr,2] -rc[lr])/rc[lr]
  #rho = rho + (-(ALL[i,2]-ALL[*tsyrs-i,2+i]))/ALL[(j)*tsyrs-i,2]
  rho = rho + (-(dtmp[i,2]-rc[i]))/rc[lr]
  print(paste(i,ntmp/i,rho))
}    

AgeFits(mod1,rec_age=1,case_label="Model 1, 2016 assessment")
(retro3$Q_Survey_1)
(retro0$Q_Survey_1)
(retro0$Sigmar)
(retro3$Sigmar)
(retro3$m_sigmar)
(retro0$m_sigmar)
(mod1$Q_Survey_1)
IndexFit(retro3,yf=1990,yl=2016,f=1,main="Model 1",ylab="Survey biomass (t)")
AgeFits(retro0,rec_age=1,case_label="retro 0 assessment")
AgeFits(retro3,rec_age=1,case_label="retro 3 assessment")

#-M-----------------------------------------------------------------------------
df <- data.frame(Year=1977:2016,mod14.1$M)
ggplot(df,aes(x=Year,y=X1))+geom_line(size=2) + mytheme + ylab("Natural mortality") + ylim(c(0,.4))



#============================================================kkkkkkkkkk
# Model selected for 2013

# used in doc to compare this year w/ last...

rec_age=1
p.rec.hist(mod1,main="Model 1",ylab="Recruitment to age",fy=1930,ly=2015,plotmean="T")
# Survey fit
IndexFit(mod1,yf=1980,yl=2015,f=1,main="Model 1", ,ylab="Survey biomass (t)")
p.catch.fit(mod1,f=1,ylim=c(0,10000))
plot(mod1$SSB,typ="b",ylab="SSB",xlab="Year",ylim=c(0,.65e6))
mod1$SSB

#pdf(paste(Figdir,"TotB.pdf",sep=""),width=9, height=7)
detach(dat)
p.biom.pol(mod1,typ="SSB",main="Model 1",new=F,fy=1977,ly=2016)
lines(mod1$TotBiom[,1],mod1$TotBiom[,2],col="black",lty=2,lwd=3)
p.biom.pol(mod1,typ="TB",main="Model 1",new=F,fy=1977,ly=2014)
#dev.off()

M <- mod16.1
ylab="xx"
xlab="xx"
main="xx"
plt_srr<-function(M,xlab="Spawning biomass",ylab="Recruits",main="Model 16.0")
{
  df <- data.frame( yrs  = M$Stock_Rec[,1] , ssb  = M$Stock_Rec[,2], rec  = M$Stock_Rec[,4])
  df2 <- data.frame(stock = M$stock_Rec_Curve[,1], rec = M$stock_Rec_Curve[,2])
  ggplot() + geom_path(data=df,aes(x=ssb,y=rec) ) +geom_text(data=df,aes(x=ssb,y=rec, label=yrs)) + 
  xlab(xlab) + ylab(ylab) + ggtitle(main)  + mytheme + geom_line(data=df2,aes(x=stock, y=rec),col="salmon",size=2) +
  theme(legend.position="none")
}
plt_srr(mod16.0)

p.stock.rec(mod16.1,main="Atka mackerel")
p.stock.rec(mod1,main="Atka mackerel")
k
p.eff.n(mod6,typ="F",f=1,main="Model")
p.eff.n(mod1,typ="S",f=1,main="Model 1")
AgeFits(mod1,rec_age=1,case_label="2015 assessment")
AgeFitsSrv(mod1,rec_age=1,case_label="2015 assessment")
#dev.off()
detach(dat)
# spawning biomass and last year's estimates 
p.biom.pol(mod1,typ="SSB",new=F,fy=1977,ly=2015)
lines(1976:2016,mod5$SSB[12:52,2],lwd=2)
#names(mod1)

# Rec
p.biom.stk(mod1,typ="R")
# Numbers at age
p.bub.age(mod1,siz=1000)
#dev.off()

p.rec.hist(mod1,fy=1977,ly=2015,main="Model 1")
lines(modsigmaR$R[,1],modsigmaR$R[,2],col="purple",lwd=2)
#dev.off()

# Likelihood table
tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
tab
rt <- get_params(lstOuts)
rt
names(mod1)
do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)}))
do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)}))

do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][40,2]),0)}))
?do.call
names(mod1)



mod1$Like_Comp
mod1$Like_Comp_names
SSB_Lastyr=read.table("clipboard")
names(mod1)

p.biom.pol(mod0,typ="SSB",new=F)
p.biom.pol(mod1,typ="SSB",new=F)
p.biom.pol(mod2,typ="SSB",new=F,main="Model 2".ly=2013)
lines(mod0$SSB[,1],mod0$SSB[,2],col="red")
lines(mod1$SSB[,1],mod1$SSB[,2],col="green")

#++++SSB CV figure=========================
plot(d3$SSB[,1],d3$SSB[,3]/d3$SSB[,2],typ="l",lty=2,ylim=c(0,.4),ylab="CV on spawning biomass",xlab="Year",cex.lab=1.4)
lines(d3$SSB[,1],d1$SSB[,3]/d1$SSB[,2],lwd=2)
lines(d3$SSB[,1],d2$SSB[,3]/d2$SSB[,2],lty=1)
lines(d3$SSB[,1],d7$SSB[,3]/d7$SSB[,2],lty=3)
legend(1968,.4, c("sigma_d=0.1","sigma_d=0.2","sigma_d=0.3", "sigma_d=1.0"),lty=c(1,1,2,3),lwd=c(2,1,1,1))

lines(modvsel$SSB[,1],modvsel$SSB[,2],col="red")
lines(mod2$SSB[,1],mod1$SSB[,2],col="purple",lwd=2)
lines(mod2$SSB[,1],mod2$SSB[,2],col="green",lwd=2)
lines(mod2$SSB[,1],mod3$SSB[,2],col="pink",lwd=2)
lines(mod2$SSB[,1],mod4$SSB[,2],col="black",lwd=2)

lines(modestM$SSB[,1],modestM$SSB[,2],col="salmon",lwd=2)
lines(modsigmaR$SSB[,1],modsigmaR$SSB[,2],col="purple",lwd=2)
modestM$Index_Q_1
names(mod1)

#++++Selectivity figure=========================
SelLastYr=c(0.00228954, 0.0333104,	0.337153,	0.879725,	1.13117,	1.32736,	1.619,	1.61509,	1.46665,	1.29413,	1.29413)
SelLastYr[[1]][1:11]
d1=readList(paste(outdir,"arc\\ds.1_R.rep",sep=""))
d2=readList(paste(outdir,"arc\\ds.2_R.rep",sep=""))
d3=readList(paste(outdir,"arc\\ds.3_R.rep",sep=""))
d7=readList(paste(outdir,"arc\\ds1.0_R.rep",sep=""))
q1.3=readList(paste(outdir,"arc\\q1.3_R.rep",sep=""))
plot(d1$N[36,3:12],typ="p",pch=19)
lines(d7$N[36,3:12])
k=36
plot(1:11,d3$sel_fsh_1[k,3:13]/max(d3$sel_fsh_1[k,3:13]),typ="l", ylab="Selectivity",xlab="Age",lwd=3, cex.lab=1.8)
lines(d1$sel_fsh_1[k,3:13]/max(d1$sel_fsh_1[k,3:13]),lty=2)
lines(d7$sel_fsh_1[k,3:13]/max(d7$sel_fsh_1[k,3:13]),lty=3)
lines(q1.3$sel_fsh_1[k,3:13]/max(q1.3$sel_fsh_1[k,3:13]),lty=3)
lines(SelLastYr/max(SelLastYr),lty=4)
SelLastYr
abline(h=.5)
legend(1,.95, c("sigma_d=0.3","sigma_d=0.1","sigma_d=1.0", "2011 Assessment"),lty=1:4,lwd=c(3,1,1,1,1))
#END ofSelectivity figure=========================


lines(mod1$TotBiom[,1],mod1$TotBiom[,2],lty=2,lwd=2)
pdf("Atka_2013.pdf",width=9, height=7)
Mntns(mod0,"Model 0")
Mntns(mod1,"Model 1")
Mntns(mod2,"Model 2")
mod1$Q_Survey_1
mod2$Q_Survey_1
dev.off()

pdf(paste(Figdir,"Selectivity.pdf",sep=""),width=7, height=11)
par(mfrow=c(1,2))
sel.age.mountain(mod1, f=1, new="F",typ="F", xvec=c(1:11),main="Model 1")
sel.age.mountain(mod2, f=1, new="F",typ="F", xvec=c(1:11),main="Model 2")
sel.age.mountain(mod2.2, f=1, new="F",typ="F", xvec=c(1:11),main="Model 2.2")
sel.age.mountain(mod1.2, f=1, new="F",typ="F", xvec=c(1:11),main="Model 2.2")
dev.off()
par(mfcol=c(1,1),mar=c(5,5,4,2) + 0.1)  

p.catch.fit(mod1,f=1,ylab="Catch biomass (t)",ylim=c(0,120000))
names(mod1)
                 
mod1$Fshry_names="Trawl"                 
spwn_ratio(mod1,main="Model 1")
cont.f.age.res(mod1, typ = "F", f = 1, lage = 1, hage = 11, cl ="COL")
p.bub.age(mod1,lage=1,hage=11,fy=1977,ly=2011,siz=100)
detach(dat)
Plot_Phase()
AgeFits(mod2,f=1,case_label="2013 assessment",rec_age=1)
AgeFits(mod1,f=1,case_label="2013 assessment",rec_age=1)
                 
modsigmaR$R
                 xlab="Age",ylab="Year",zscale=2.5,new=F,cex.yax=1.,fy=1980)
detach(dat)
IndexFit(mod,yf=1980,yl=2010,f=2,main=main)

par(mfcol=c(1,1))

# show fit to catch biomass
CatchFit(mod2.2)

# show spawning biomass relative to population with no fishing
spwn_ratio(mod1,fy=1962,ly=2015) 
fix(spwn_ratio) 

detach(dat)
# example of writing multiple plots to pdf file:
#pdf("figs/agefits.pdf",width=9, height=7)
  #AgeFits(am1,f=1)
#dev.off()

# another example of writing multiple plots to pdf file:
#pdf("figs/survey_fit.pdf",width=7, height=9)
#Mntns(am1,"Model 1")
    # Indices(am1,"Model 1")
par(mfrow=c(1,1))
IndexFit(mod1,main="Model 1",yf=1990,ylab="Survey biomass (t)")
detach(dat) 
Plot_Fspr()
source("R/Plot_Atkas.R")
#spwn_ratio(am1) 
#plt_proj(am1)
#dev.off()

dev.off()

p.biom.pol(retro2,typ="SSB",main="Model 1",new=F,fy=1977,ly=2013)
for (i in 1:9) {
  rn=paste0("retro1",i);
  lines(get(rn)$SSB[,1],get(rn)$SSB[,2],col=i)
}
pdf(paste(Figdir,"Retro_Mod2.pdf",sep=""),width=7, height=9)
par(mfrow=c(2,1))
p.biom.pol(retro0,typ="SSB",main="Model 2",new=F,fy=1977,ly=2012)
#plot(retro0$SSB[,1],retro0$SSB[,2],ylim=c(0,550),      ylab="Spawning biomass (kt)",type="l",lwd=2,xlab="",lty=2)
ssb=1966:2013
retro0$R
names(retro1)
rrr=1977:2012
for (i in 0:10) {
  rn=paste("retro",i,sep="");
  lines(get(rn)$SSB[,1],get(rn)$SSB[,2],col=i)
  ssb=rbind(ssb,c(t(get(rn)$SSB[,2]),rep(NA,i)))
  rrr=rbind(rrr,c(t(get(rn)$R[,2]),rep(NA,i)))
  }
write.csv(ssb,"Atka_SSB.csv")
write.csv(rrr,"Atka_rec.csv")

system("atka_ssb.csv")
?write.csv
rrr
rep(NA,2)
plot(retro0$SSB[,1],rep(0,length(retro0$SSB[,1])),
     ylim=c(-.7,.7),
     ylab="Relative difference from terminal year",
     type="l",xlab="Year",lty=2,lwd=2)
for (i in 1:10) {
  rn=paste("retro",i,sep="");
  lines(get(rn)$SSB[,1],
        (get(rn)$SSB[,2]/retro0$SSB[1:(48-i),2])-1,col=i)
  }


dev.off()
p.biom.pol(mod2,typ="SSB",n.mod=1,main="Model 2",new=F,fy=1977,ly=2013)

p.biom.pol(mod2,typ="SSB",n.mod=2,mod1,main="Model 1",new=F,fy=1977,ly=2013)
lines(SSB_Lastyr[,1],SSB_Lastyr[,2]/2,lwd=2,col="red")

# Plot selectivity in multiple crappy panels
p.select.hist(mod1,typ="F",h="T",f=1,lage=1,hage=11,fy=1985,ly=2000)
# Plot selectivity in multiple crappy color grayscale 
c.select(mod1)
dev.off()  
