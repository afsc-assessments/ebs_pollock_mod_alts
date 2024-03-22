R
rm(list=ls())
ls()
source("R/prelims.R")
#-------------------------------------------------------------------------------
# Visual compare runs
#-------------------------------------------------------------------------------
library(ggridges)

#source("../../R/compareRuns.r")
source("../R/compareRuns.r")
mytheme <- mytheme + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank() )



#--Projections---------
pdf("proj.pdf")
pdt <- data.table(read.table("proj/bigfile.out",header=TRUE)) 
pdt <- data.table(read.table("arc/mod16.0b_N50_proj.rep",header=TRUE)) 
pdt$N <- 50
pd <- data.table(read.table("arc/mod16.0b_N100_proj.rep",header=TRUE)) 
pd$N <- 100
pdt <- rbind(pdt,pd)
pd$N <- 200
pdt <- rbind(pdt,pd)
pdt$Alternative <- as.factor(pdt$Alternative)
pdt$N <- as.factor(pdt$N)
names(pdt)

pt <- pdt[Yr>2018&Alternative==1,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB),B35=mean(B35),B0=mean(B0),B40=mean(B40),lb=quantile(SSB,.2) ,ub=quantile(SSB,.8) ),.(Yr)]

s1 <- ggplot(pt,aes(x=Yr,y=SSB)) + geom_line() + mytheme + ylim(c(0,100)) + geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.25) + labs(y="Spawning biomass (kt)",x="Year") + scale_x_continuous(breaks=seq(2015,2032,2)) + geom_line(aes(y=B40),color="blue",linetype="dashed") 
s1
pt <- pdt[Yr>2018&Alternative==1,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(Catch,.2) ,ub=quantile(Catch,.8) ),.(Yr)]
c1 <- ggplot(pt,aes(x=Yr,y=Catch)) + mytheme + labs(y="Catch (kt)",x="Year") + scale_x_continuous(breaks=seq(2015,2032,2))  + ylim(c(0,120)) +
 geom_ribbon(aes(ymin=lb,ymax=ub),fill="salmon",color="red",alpha=0.25) + geom_line(aes(x=Yr,y=ABC),size=1)
 c1

pt <- pdt[Yr>2016&Alternative==1,.(B35=mean(B35),B0=mean(B0),B40=mean(B40),SSB=median(SSB) ,lb=quantile(SSB,.2) ,ub=quantile(SSB,.8) ),.(Yr,N)]
s2 <- ggplot(pt,aes(x=Yr,y=SSB,linetype=N,color=N)) + geom_line(aes(linetype=N),size=1.2) + mytheme + ylim(c(0,100)) + geom_line(aes(y=B40,linetype=N)) + 
      labs(y="Spawning biomass (kt)",x="Year") + scale_x_continuous(breaks=seq(2015,2032,2)) 
s2

pt <- pdt[Yr>2016&Alternative==1,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(Catch,.2) ,ub=quantile(Catch,.8) ),.(Yr,N)]
c1 <- ggplot(pt,aes(x=Yr,y=Catch,color=N)) + mytheme + labs(y="Catch (kt)",x="Year") + scale_x_continuous(breaks=seq(2015,2032,2))  + ylim(c(0,120)) +
 geom_ribbon(aes(ymin=lb,ymax=ub,fill=N),alpha=0.25) #+ geom_line(aes(x=Yr,y=ABC),size=1)
 c1
pt <- pdt[Yr>2016&Alternative==1,.(F=mean(F),FABC=mean(FABC),FOFL=mean(FOFL),SSB=median(SSB) ,lb=quantile(F,.2) ,ub=quantile(F,.8) ),.(Yr,N)]
f1 <- ggplot(pt,aes(x=Yr,y=F,color=N)) + mytheme + labs(y="Fishing mortality ",x="Year") + scale_x_continuous(breaks=seq(2015,2032,2))  + ylim(c(0,.5)) +
 geom_line() # + geom_line(aes(x=Yr,y=FABC),size=1)+ geom_line(aes(x=Yr,y=FOFL),size=1)
 geom_ribbon(aes(ymin=lb,ymax=ub,fill=N),alpha=0.25) # + geom_line(aes(x=Yr,y=FABC),size=1)+ geom_line(aes(x=Yr,y=FOFL),size=1)
 f1
c1 <- c1 + geom_line(aes(x=Yr,y=OFL),size=2)
#c1 <- c1 + geom_line(data=pt[as.numeric(Alternative)==2,.(Yr,ABC)],aes(x=Yr,y=ABC))
c1
pt[as.numeric(Alternative)==2,.(Yr,ABC)]
pt <- pdt[Yr>2016,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL)),.(Yr,Alternative)] 
pt <- pdt[Yr>2016&Alternative==1,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(SSB,.2) ,ub=quantile(SSB,.8) ),.(Yr,N)]
pt
ggplot(pt,aes(x=Yr,y=OFL,color=Alternative)) + geom_line() + mytheme
pdt
pdx <-rbind(pdt2,pdt)
setkey(pdx,Yr,Alternative)
pt <- pdx[.(Yr>2016,(Alternative)==1),.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB)  ),.(Yr,config)]
pt <- pdx[Yr>2016&Alternative==1,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(Catch,.1) ,ub=quantile(Catch,.9)  ),.(Yr,config)]
ggplot(pt,aes(x=Yr,y=Catch,color=config,shape=config)) + geom_line() + geom_point() + mytheme + geom_ribbon(aes(ymin=lb,ymax=ub,fill=config),alpha=0.25) + labs(x="Year")+ scale_x_continuous(breaks=seq(2015,2032,2))  + ylim(c(0,130))
dev.off()

nll <- data.frame()
dtmp <- data.table(read.table("arc/mod16.0b_N50_NLL.rep")) 
dtmp$N <- (50) 
nll <- rbind(nll,dtmp)
dtmp <- data.table(read.table("arc/mod16.0b_N100_NLL.rep")) 
dtmp$N <- 100 
nll <- rbind(nll,dtmp)
dtmp <- data.table(read.table("arc/mod16.0b_N200_NLL.rep")) 
dtmp$N <- 200 
nll <- rbind(nll,dtmp)

names(dtmp) <- c("Year","Survey","Fishery_age","Survey_age")
ggplot(nll,aes(x=Year,y=Survey,color=Terminal_year)) + geom_line() + geom_point()+ mytheme + xlim(c(1990,2017))+ theme(legend.position="none") + ylab("Survey index -log likelihood")
ggplot(nll,aes(x=Year,y=Survey_age,color=Terminal_year)) + geom_line() + geom_point() + mytheme+ theme(legend.position="none") + ylab("Survey age -log likelihood")

