source("../../../R/prelims.R")
getwd()
source("../../R/plot_bts.R")
library(ggthemes)
library(tidyverse)
library(r4ak)
.THEME=theme_few()
#-------------------------------------------------------------------------------
# Visual compare runs
#-------------------------------------------------------------------------------
library(ggridges)

source("../../R/compareRuns.r")

#        Read in the output of the assessment
mod    <- read_admb("base/amak")
A = mod
#mod16.0b <- readList("mod16.0b/For_R.rep")
#lstOut1  <- list( "2017.Assessment"= modlyr, "Current.assessment"= mod16.0b)
M <- list(mod)
plot_bts(M,idx=2)
#tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
#tab       <- cbind(lstOuts[[1]]$,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
#tt <- tibble( srv_sdnr = lstOuts %>% map("sdnr_ind_1"),
  #fsha_sdnr =
mytheme = mytheme + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
compareTime(lstOut1,"R",SD=T,Sum=NULL,legendPos="top",startYear=1980)
rdf <- cbind(data.table(mod$R,"Model_16.0"))
names(rdf) <- c("yr","R","se","lb","ub","case")
rdf  <- rdf[yr>1990&yr<2017,.(yr,R=R,lb=lb,ub=ub,case)]
mnR <- mean(mod$R)
dodge <- position_dodge(width=0.8)
ggplot(rdf,aes(x=yr-1,y=R,fill=case)) + xlab("Year class") + ylab("Age 1 recruits (thousands)") + #ylim(c(0,18000)) +
       geom_bar(width=0.75,position="dodge",stat="identity",color="black") + 
       scale_x_continuous(breaks=seq(1990,2017,2)) +
       geom_errorbar(aes(ymin=lb,ymax=ub),width=.3,colour="blue",position=dodge) + mytheme + geom_hline(aes(yintercept=mnR))
#-------------------------------------------------------------------------------
# Fit to survey data
mdf <- .get_bts_df(M,idx=2)
mdf2 <- mdf[!is.na(mdf$obs),]
ggplot(mdf,aes(x=year,y=pre)) + geom_line(width=2,color="blue") + .THEME + geom_point(data=mdf2,aes(x=year,y=obs),size=2,color="red") + expand_limits(y = 0) + 
                            geom_errorbar(data=mdf2,aes(x=year,ymax=ub,ymin=lb),width=.5) + ylab("Survey biomass (t)")+ xlab("Year") + 
                            scale_x_continuous(breaks=seq(1970,2019,2), limits=c(1980,2019)) 
AgeFitsSrv(mod,rec_age=1,case_label="ebs pollock")

#-------------------------------------------------------------------------------
# Fit to fishery data
AgeFits(mod,rec_age=1,case_label="Fishery")


#-------------------------------------------------------------------------------
#Selectivity
df <- data.frame(mod$sel_fsh_1[,2:13] ); names(df) <- c("yr",1:11)
sdf <- gather(df,age,sel,2:12) %>% filter(yr>1976) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .2)+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)") #+ theme_ridges()  

mytheme <- mytheme + theme(text=element_text(size=14)) + theme(axis.title.x=element_text(size=24) ,axis.title.y=element_text(size=24))
df <- data.frame(mod$sel_fsh_1[,2:13] ); names(df) <- c("yr",1:11)
(M$sel_fsh_1)
names(M)
sdf <- gather(df,age,sel,2:12) %>% filter(yr>1976) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
head(sdf)
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .9,color="blue",fill="yellow",size=.5)+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr))))
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .4,color="blue",fill="yellow",size=.5)+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr))))
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .4,color="black",fill="orange",size=.5)+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr))))
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .2,fill="orange")+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)")+ xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr)))) 

rec_age=2
# Survey fit
IndexFit(mod,yf=1980,yl=2018,f=1,main="Model 1", ,ylab="Survey biomass (t)")
p.catch.fit(mod,f=1,ylim=c(0,1800))
p.biom.pol(mod,typ="SSB",main="Model x",new=F,fy=1977,ly=2018)
lines(1976:2017,mod$SSB[12:53,2],lwd=2,lty=2)
p.biom.pol(mod,typ="TB",main="Model 16.0b",new=F,fy=1977,ly=2018)


M <- mod
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
plt_srr(mod,main="Model 19")
