rm(list=ls())
ls()
source("R/prelims.R")
.THEME=mytheme
#-------------------------------------------------------------------------------
# Visual compare runs
#-------------------------------------------------------------------------------
library(ggridges)
library(ggthemes)
library(scales)

source("R/compareRuns.r")

#        Read in the output of the assessment
modlyr   <- readList("2019_Final/For_R.rep")
mod16.0b <- readList("mod16b/For_R.rep")
mod20 <- readList("m20/For_R.rep")
mod3par <- readList("dbl_log/For_R.rep")
lstOut1  <- list( "Non-parametric"= mod16.0b, "3-parameter double logistic"= mod3par)
M <- lstOuts <- lstOut1 
#tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
#tab       <- cbind(lstOuts[[1]]$,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
#tt <- tibble( srv_sdnr = lstOuts %>% map("sdnr_ind_1"),
  #fsha_sdnr = lstOuts %>% map("sdnr_age_fsh_1"),
  #srva_sdnr = lstOuts %>% map("sdnr_age_ind_1")
    #)
#map_df(data.table(lstOut1),"[","sumbiom")
caplyr <- "2019 assessment"
captyr <- "Current assessment "

#=====================
#   SSB 
#=====================
df  <- data.table(Model = "16.0b", mod16.0b$SSB )
df <- rbind(df, data.table(Model="lastyr",modlyr$SSB))

df  <- data.table(Model = "Non-parametric", mod16.0b$SSB )
df <- rbind(df, data.table(Model="3-parameter",mod3par$SSB))
# for (i in 2:3) df <- rbind(df, data.table(Model=paste0("Model ",i),lstOuts[[i]]$SSB))
names(df) <- c("Model","yr","SSB","SE","lb","ub")
bdf <- filter(df,yr>1980,yr<=2021) %>% arrange(yr);bdf

p1 <- ggplot(data=bdf,aes(x=yr,y=SSB,alpha=0.3,color=Model)) + scale_y_continuous(limits=c(0,1000000)) + ylab("Spawning biomass") + 
          xlab("Year") +  theme_few(base_size=19) + geom_line(data=bdf,aes(x=yr,y=SSB,color=Model)) +
          geom_ribbon(data=bdf ,aes(x=yr,y=SSB,ymin=lb,ymax=ub,fill=Model),alpha=.3)  + guides(alpha=FALSE,col=FALSE) ;p1
p1
#=====================
#   Rec 
#=====================
rdf <- cbind(data.table(mod16.0b$R,"Model_16.0b"))
rdf <- rbind(rdf,cbind(data.table(modlyr$R,"2019 Assessment")))
names(rdf) <- c("yr","R","se","lb","ub","case")
mytheme = mytheme + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
rdf  <- rdf[yr>1990&yr<2020,.(yr,R=R,lb=lb,ub=ub,case)]
mnR <- 580 #mean(mod16.0b$R)
dodge <- position_dodge(width=0.8)
ggplot(rdf,aes(x=yr-1,y=R,fill=case)) + xlab("Year class") + ylab("Age 1 recruits (thousands)") + #ylim(c(0,18000)) +
       geom_bar(width=0.75,position="dodge",stat="identity",color="black") + 
       scale_x_continuous(breaks=seq(1990,2018,2)) +
       geom_errorbar(aes(ymin=lb,ymax=ub),width=.3,colour="blue",position=dodge) + mytheme + geom_hline(aes(yintercept=mnR))

#=====================
#    Fishing mortality
#=====================
     library(data.table)
M <- mod16.0b
fdf <- data.table(M$F_fsh_1)
fdf$C_B <- M$Obs_catch_1/M$TotBiom[1:length(M$Obs_catch_1),2]
names(fdf) <- c("Year","MeanF","F","C_B")
fdf.m <-melt(fdf,id="Year",value="Estimate",variable="Type")
fdf.m
ggplot(fdf.m,aes(x=Year,y=Estimate,color=Type)) + xlab("Year") +  
       geom_line(size=2.) + scale_x_continuous(breaks=seq(1977,2019,3)) + mytheme  
       geom_hline(aes(yintercept=mnR),linetype="dashed") + theme(legend.position="none")


#=====================
# Fit to survey data
#=====================
AgeFitsSrv(mod16.0b,rec_age=1,case_label=captyr)

mdf <- .get_bts_df(lstOut1[2])
mdf2 <- mdf %>% filter(!is.na(obs))
ggplot(mdf,aes(x=year,y=pre)) + geom_line(width=2,color="blue") + .THEME + geom_point(data=mdf2,aes(x=year,y=obs),size=2,color="red") + expand_limits(y = 0) + 
                            geom_errorbar(data=mdf2,aes(x=year,ymax=ub,ymin=lb),width=.5) + ylab("Survey biomass (t)")+ xlab("Year") + 
                            scale_x_continuous(breaks=seq(1970,2019,2), limits=c(1990,2019)) 
IndexFit(mod16.0b,yf=1990,yl=2019,f=1,main=captyr,ylab="Survey biomass (t)")
AgeFitsSrv(modlyr,rec_age=1,case_label=caplyr)

AgeFitsSrv(mod16.0b,rec_age=1,case_label="16.0b")

#=====================
# Fit to fishery data
#=====================
AgeFits(mod16.0b,rec_age=1,case_label=captyr)
#AgeFits(mod20,rec_age=1,case_label="Model 20.0")
#AgeFits(modlyr,rec_age=1,case_label=caplyr)

dt <- data.table(modlyr$Obs_Survey_1[,1:4],mod16.0b$Obs_Survey_1[,3])
names(dt) <- c("Year","Observed","Model_lastyr","sd","Model_16.0b")
dt <- data.table(dt)[Year>1986,.(Year, Observed, Model_lastyr,Model_16.0b, lb=Observed-1.96*sd,ub=Observed+1.96*sd),]
dtp <-melt(dt,id="Year")
dtp
ggplot(dt,aes(x=Year,y=Model_16.0b) ) + geom_point(data=dtp[variable=="Observed"],size=3) + geom_line(data=dtp[variable!="Observed"]) + ylim(c(0,1800000)) + labs(x="Year",y="Survey biomass") + mytheme #+
,aes(x=Age,y=value,colour=variable)) 
dt[variable!="ub"&variable!="lb"] %>% ggplot(aes(x=Year,y=value,col=variable)) + geom_point(data=dt[variable=="Observed"],size=3) + geom_line(data=dt[variable!="Observed"]) + ylim(c(0,1800000)) + labs(x="Year",y="Survey biomass") + mytheme #+
            geom_line(data=)

#=====================
#Selectivity
#=====================
tt <- data.frame(mod16.0b$sel_fsh_1[,3:13]) %>% mutate(mak = do.call(pmax, (.))) 
tt
df <- data.frame(mod16.0b$sel_fsh_1[,2], mod16.0b$sel_fsh_1[,3:13] ) ; names(df) <- c("yr",1:11)
df %>% mutate(max=max(c_across(2:12))) %>% mutate(2:12/max)

df <- data.frame(mod3par$sel_fsh_1[,2:13] ); names(df) <- c("yr",1:11)
df <- data.frame(mod3par$sel_fsh_1[,2:13] ); names(df) <- c("yr",1:11)
df
sdf <- gather(df,age,sel,2:12) %>% filter(yr>1976) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
max(sdf$sel )
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .9,color="blue",fill="yellow",size=.5)+ xlim(c(1,11))+ ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr)))) + theme_few()
#ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .4,color="blue",fill="yellow",size=.5)+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr))))
#ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .4,color="black",fill="orange",size=.5)+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr))))
#ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .2,fill="orange")+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)")+ xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr)))) 

#=====================
# SRR
#=====================
M <- mod16.0b
xlab="Spawning biomass"
ylab="Recruits age 1 (thousands)"
library(scales)
main="Model 16.0b"
plt_srr <-function(M,xlab="Spawning biomass",ylab="Recruits",main="Model 16.0b") {
  df <- data.frame( yrs  = M$Stock_Rec[,1] , ssb  = M$Stock_Rec[,2], rec  = M$Stock_Rec[,4])
  df2 <- data.frame(stock = M$stock_Rec_Curve[,1], rec = M$stock_Rec_Curve[,2])
  ggplot() + geom_path(data=df,aes(x=ssb,y=rec) ) +geom_text(data=df,aes(x=ssb,y=rec, label=yrs)) + 
  xlab(xlab) + ylab(ylab) + ggtitle("Atka mackerel")  + theme_few(base_size=16) + geom_line(data=df2,aes(x=stock, y=rec),col="salmon",size=2) + theme_few() + 
   scale_y_continuous(labels = comma) + scale_x_continuous(labels = comma) 
}
plt_srr(mod16.0b,main="Model 16.0b")

#-------------------------------------------------------------------------------
p.eff.n(mod16.0b,typ="F",f=1,main="Model")
rec_age=1
p.rec.hist(mod16.0b,fy=1977,ly=2018,main="Atka mackerel, Model 16.0b")

p.eff.n(mod16.0b,typ="S",f=1,main="Model")
detach(dat)

compareTime(lstOut1,"SSB",SD=T,Sum=NULL,legendPos="top",startYear=1980)
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,startYear=2000,ylim=c(0,3.4e5))
compareTime(lstOuts,"R",SD=F,ylim=c(0,1500),startYear=2000)
compareTime(lstOuts,"R",SD=T,cex.lab=2)
compareTime(lstOuts,"TotBiom",SD=T)
compareMatrix(lstOuts,"TotF",SD=TRUE,Sum=NULL,YrInd=modlyr$Yr,Apply=mean,legendPos="right")
#-------------------------------------------------------------------------------


#Survey selectivity==============================================
srv_sel1 <- modlyr$sel_ind_1[1,3:13]
srv_sel1 <- data.frame(Age=1:11,Selectivity=srv_sel1/max(srv_sel1))
srv_sel1$Model <-caplyr

#=====================
# Retro
# get all retrospectives
#=====================
# Read in retro results
for (i in 0:15) {
  rn=paste0("mod16b/retro/r_",i,".rep")
  mn=paste0("retro",i)
  assign(mn,readList(rn))
  print(rn)
}
retouts <- list()
retouts <- list( R0=retro0, R1= retro1, R2= retro2, R3= retro3, R4= retro4 , R5= retro5 , R6= retro6 , R7= retro7 , R8= retro8 , R9= retro9 , R10= retro10 )  
df  <- data.table(Model = "16.0b", mod16.0b$SSB )
names(df) <- c("Model","yr","SSB","SE","lb","ub")
bdf <- filter(df,yr>1977) %>% arrange(yr)
bdf
ssb <- bdf %>% transmute(SSB)
ssb
mytheme = mytheme + theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0))
p1 <- ggplot() + scale_y_continuous(limits=c(0,530000)) + ylab("Spawning biomass (t)") + xlab("Year") +  mytheme + geom_line(data=bdf,aes(x=yr,y=SSB),size=4) +
        scale_y_continuous(labels = comma) +
           geom_ribbon(data=bdf ,aes(x=yr,ymin=lb,ymax=ub),fill="gold",col="grey",alpha=.6)  + guides(fill=FALSE,alpha=FALSE,col=FALSE) + scale_x_continuous(breaks=seq(1978,2020,2)) 
           p1
           i=8
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
# Relative... 
p1 <- ggplot() + scale_y_continuous(limits=c(0,2)) + ylab("Relative spawning biomass") + xlab("Year") +  
mytheme + geom_line(data=bdf,aes(x=yr,y=1. ),size=2) 
           p1
           i=8
for (i in 1:10) {
  rn=paste("retro",i,sep="");
  tmp= (data.frame(get(rn)$SSB)) 
  names(tmp) <- c("yr","SSB","SE","lb","ub"); 
  tmp= tmp %>%filter(yr>1977)
  tmp$SSB <- tmp$SSB/bdf$SSB[1:length(tmp$SSB)]
  p1 <- p1 + geom_line(data=tmp,aes(x=yr,y=SSB),col=i,linetype=i,size=1.25)
  #p1 <- p1 + geom_segment(data=tdf,aes(x=yr,xend=yr,yend=SSB,y=SSB),arrow=arrow(angle=90,length=unit(.2,"cm")),size=2,col=i)
  tdf <- tmp[dim(tmp)[1],]
  p1 <- p1 + geom_point(data=tdf,aes(x=yr,y=SSB),size=4,col=i)
  #p1 <- p1 + geom_point(get(rn)$SSB[lr,1],get(rn)$SSB[lr,2],pch=19,col=i)
}
p1

tab = list(R0=retro0$SSB)
for (i in 1:10) { tab=cbind(tab,retouts[[i]]$SSB)}
tab
#pdf(paste(Figdir,"Retro_Mods.pdf",sep=""),width=9, height=6)
bdf <- bdf[,-1]
bdf
# Color blind palette
# cb_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# p1 <-    p1 + scale_fill_manual(values=cb_palette)
df <- data.frame(mod16.0b$SSB )
names(df) <- c("yr","SSB","SE","lb","ub")
bdf2<- filter(df,yr>1977)
bdf <- filter(df,yr>1977,yr<=2020)
bdf
bdft <- bdf
for (i in 1:10) bdft <- cbind(bdft,rep(NA,4))
  dim(bdft)
  head(bdft)
  tail(bdft)
names(bdft)[6:15] <- paste0("SSB_",2020:2005)
i=1
get(paste0("retro",i))$SSB[14:(55-i),1:2]
for (i in 1:10) bdft[1:(44-i),i+5] <- get(paste0("retro",i))$SSB[14:(55-i),2]
bdft
bdf$SSB
for (i in 1:10) bdft[,i+5] <- bdft[,i+5]/bdft[,2]
bdft
p2 <- ggplot(bdft,aes(x=yr,y=SSB)) + scale_x_continuous(breaks=seq(1978,2021,2)) + 
      scale_y_continuous(limits=c(.25,1.75)) + ylab("Relative spawning biomass") + xlab("Year") +  mytheme  
for (i in 1:10) {
  tdf <- data.frame(cbind(bdft[1],SSB=bdft[5+i]))[1:(43-i),]
  names(tdf) <- c("yr", "SSB")
  p2 <- p2 + geom_line(data=tdf, aes(x=yr,y=SSB),col=i,size=1.5)
  tdf <- tdf[dim(tdf)[1],]
  p2 <- p2 + geom_point(data=tdf,aes(x=yr,y=SSB),size=4,col=i)
}    
#p2 <-p2 +geom_hline(aes(yintercept=1),size=3,linetype="dotted")
p2 <-p2 +geom_hline(aes(yintercept=1),size=1,col="grey")
p2
p1
grid.arrange(p1,p2,nrow=2)
# Mohn's rho
ref    = retro0$SSB[,2]
rho    = 0
for (i in 1:10){
  dtmp   = get(paste0("retro",i))$SSB[,2]
  tip    = length(dtmp)
  rho    = rho + (dtmp[tip] -ref[tip])/ref[tip]
  print(paste(i,rho/i))
}    

AgeFits(mod16.0b,rec_age=1,case_label="Model 16.0b, 2018 assessment")
(retro3$Q_Survey_1)
(retro0$Q_Survey_1)
(retro0$Sigmar)
(retro3$Sigmar)
(retro3$m_sigmar)
(retro0$m_sigmar)
(mod1$Q_Survey_1)
AgeFits(retro0,rec_age=1,case_label="retro 0 assessment")
AgeFits(retro3,rec_age=1,case_label="retro 3 assessment")

#-M-----------------------------------------------------------------------------
df <- data.frame(Year=1977:2016,mod14.1$M)
ggplot(df,aes(x=Year,y=X1))+geom_line(size=2) + mytheme + ylab("Natural mortality") + ylim(c(0,.4))



#============================================================
# Model selected for 2013

# used in doc to compare this year w/ last...

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

lines(modsigmaR$R[,1],modsigmaR$R[,2],col="purple",lwd=2)
#dev.off()

# Likelihood table
lstOuts <- lstOut1
tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)}))
do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)}))
tab  <- rbind(
cbind("q",do.call(cbind,lapply(lstOuts,function(x){round(x[["q_1"]][1,2],2)}))),
cbind("Npars",do.call(cbind,lapply(lstOuts,function(x){round(x[["Num_parameters_Est"]],2)}))),
cbind("M",do.call(cbind,lapply(lstOuts,function(x){round(x[["Mest"]][2],2)}))),
cbind("SigmaR",do.call(cbind,lapply(lstOuts,function(x){round(x[["Sigmar"]][2],2)}))),
cbind("EffN_Fish",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["EffN_Fsh_1"]][,2]),2)}))),
cbind("EffN_Surv",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["EffN_Survey_1"]][,2]),2)}))),
tab,
cbind("F2016",do.call(cbind,lapply(lstOuts,function(x){round((x[["F_fsh_1"]][40,2]),2)}))),
cbind("F2016/F40%",do.call(cbind,lapply(lstOuts,function(x){round((x[["F_fsh_1"]][40,2])/(x[["F40_est"]]),2)}))),
cbind("B 1977",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][1,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][1,3])/(x[["TotBiom"]][1,2])*100,0)}))),
cbind("B 2016",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][40,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][40,3])/(x[["TotBiom"]][38,2])*100,0)}))),
cbind("2006 YC",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][31,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][31,3])/(x[["R"]][31,2])*100,0)}))),
cbind("2012 YC",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][37,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][37,3])/(x[["R"]][37,2])*100,0)})))
)
tab
i=2
names(lstOuts)
for (i in 1:3){
  df <- data.table( lstOuts[[i]]$Obs_Survey_1[-1,] )
  names(df) <- c("Year","obs","pre","sd","pr","lnpr")
  df$pre <- as.numeric(df$pre)
  df$obs <- as.numeric(df$obs)
  df <- df[df$obs!="NA",]
  print( df %>% summarize(sqrt(mean(log(obs/pre)^2))) )
}
  df %>% filter(!is.na(obs)) 
tab
  lstOuts %>% map(
ttt <- 


tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
tab
cbind(lstOut1[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOut1,function(x){round(x[["Like_Comp"]],2)})))
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
#lstOuts   <- list( Model_14.1= mod14.1, Model_16.0= mod16.0, Model_16.1= mod16.1, Model_16.2= mod16.2)
#lstOuts   <- list( Model_14.1= mod14.1, Model_16.0= mod16.0, Model_16.1= mod16.1)
#lstOuts   <- list( Model_0.0=mod0.0, Model_0.1=mod0.1 , Model_0.2= mod0.2, Model_0.3= mod3)

tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
tab
write.csv(tab,file=file.path(inputPath,"LikelihoodTable2017.csv"),row.names=F)
system(paste("open ",inputPath,"/LikelihoodTable2015.csv",sep=""))
#----parameter table-----------------
rt <- get_params(lstOuts)
rt
get_params

#--Projections---------
pdt <- data.table(read.table(paste0("mod16b/proj/bigfile.out"),header=TRUE)) 
names(pdt)
#[1] "Alternative" "SpNo"        "Spp"         "Yr"          "ABC"        
# [6] "OFL"         "Catch"       "SSB"         "F"           "Tot_biom"   
#[11] "SPR_Implied" "Ntot"        "SexRatio"    "FABC"        "FOFL"       
#[16] "B0"          "B40"         "B35"        
library(ggthemes)
library(scales)
head(pdt)
nsims=1000
sim=rep(rep(1:14,nsims),7)
length(sim)
dim(pdt)
pdt$Sim=sim
tt <- (sample(nsims, size=30))
mtmp <- pdt %>% filter(Sim %in% tt) %>% select(Yr,sim,SSB,Alternative) 
#p1 <- p1 + facet_grid(Area~.,scales="free") + geom_line(data=mtmp,aes(group=sim),size=.2,col=SSB) 

p1 <- pdt %>% filter(Alternative==1) %>% mutate(Year=Yr) %>% ggplot(aes(x=Year,y=Catch,fill=Alternative)) + 
        expand_limits(y = 0) + 
        scale_y_continuous(labels = comma) +
        ggplot2::stat_summary(fun.ymin = function(x) quantile(x, 0.05), fun.ymax = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25, colour = NA) +
        ggplot2::stat_summary(fun.y = function(x) quantile(x, 0.5), geom = "line", lwd = 1) +
        ggplot2::stat_summary(fun.y = function(x) quantile(x, 0.5), geom = "point", size = 1) +
        theme_few()+ labs(x="Year",y='SSB (t)')    #+ geom_line(data=pdt,aes(x=Yr,y=B40))
       p1 <- p1 + geom_line(data=mtmp,aes(x=Yr,y=SSB),size=.2,col="grey40") 
p1
p1 <- pdt %>% filter(Alternative==1) %>% mutate(Year=Yr) %>% ggplot(aes(x=Year,y=SSB,fill=Alternative)) + 
        expand_limits(y = 0) + 
        scale_y_continuous(labels = comma) +
        ggplot2::stat_summary(fun.ymin = function(x) quantile(x, 0.05), fun.ymax = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25, colour = NA) +
        ggplot2::stat_summary(fun.y = function(x) quantile(x, 0.5), geom = "line", lwd = 1) +
        ggplot2::stat_summary(fun.y = function(x) quantile(x, 0.5), geom = "point", size = 1) +
        theme_few()+ labs(x="Year",y='SSB (t)')    + geom_line(data=pdt,aes(x=Yr,y=B40))
p1

pdt1 <- data.table(read.table(paste0("proj/",pfn,"_out/catch_fixed_2017_18.out"),header=TRUE)) 
pdt1[Yr>2016&Yr<2019,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(SSB,.2) ,ub=quantile(SSB,.8) ),.(Yr,Alternative)]
pfn <- "10yr"
pdt2 <- data.table(read.table(paste0("proj/",pfn,"_out/bigfile.out"),header=TRUE)) 
pdt$Alternative <- as.factor(pdt$Alternative)
pdt2$Alternative <- as.factor(pdt2$Alternative)
pdt$config <-"5yr"
pdt2$config <-"10yr"

pt <- pdt[Yr>2019,.(Catch=mean(Catch),ABC=mean(ABC),OFL=mean(OFL),SSB=median(SSB) ,lb=quantile(SSB,.2) ,ub=quantile(SSB,.8) ),.(Yr,Alternative)]
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

  # write some sims out for easy use later...
  bfs <- pdt %>% filter(Sim<=30)
  head(bfs)
  write.csv(bfs,"data/proj.csv")
 # head(bfs)
  bfss <- bfs %>% filter(Alternative==2) %>% transmute(Alternative,Yr,Catch,SSB,Sim=as.factor(Sim)) 
  pf <- data.frame(read.table("mod16b/proj/percentdb.out",header=F) )
  names(pf) <- c("stock","Alt","Yr","variable","value") 
  thisyr=2020
  dev.off()
  p1 <- pf %>% filter(substr(variable,1,1)=="C",variable!="CStdn",Alt==2) %>% select(Yr,variable,value) %>% spread(variable,value) ;p1#%>%
    ggplot(p1,aes(x=Yr,y=CMean),width=1.2) + geom_ribbon(aes(ymax=CUCI,ymin=CLCI),fill="goldenrod",alpha=.5) + theme_few() + geom_line() +
    scale_x_continuous(breaks=seq(thisyr,thisyr+14,2))  +  xlab("Year") + ylim(0,200000) + ylab("Tier 3 ABC (kt)") + geom_point() + 
    geom_line(aes(y=Cabc)) + geom_line(aes(y=Cofl),linetype="dashed") + geom_line(data=bfss,aes(x=Yr,y=Catch,col=as.factor(Sim)))+ guides(size=FALSE,fill=FALSE,alpha=FALSE,col=FALSE) 
  p2 <- pf %>% filter(substr(variable,1,1)=="S",variable!="SSBStdn",Alt==2) %>% select(Yr,variable,value) %>% spread(variable,value) %>%
    ggplot(aes(x=Yr,y=SSBMean),width=1.2) + geom_ribbon(aes(ymax=SSBUCI,ymin=SSBLCI),fill="coral",alpha=.5) + theme_few() + geom_line() +
    scale_x_continuous(breaks=seq(thisyr,thisyr+14,2))  +  xlab("Year") + ylim(0,5000) + ylab("Tier 3 Spawning biomass (kt)") + geom_point() + 
    geom_line(aes(y=SSBFabc)) + geom_line(aes(y=SSBFofl),linetype="dashed")+ geom_line(data=bfss,aes(x=Yr,y=SSB,col=as.factor(Sim)))+ guides(size=FALSE,fill=FALSE,alpha=FALSE,col=FALSE) 
  #t3 <- grid.arrange(p1, p2, nrow=2)
    library(patchwork)
    p1
  t3 <- p1/ p2
  print(t3)
  ggsave("figs/tier3_proj.pdf",plot=t3,width=5.4,height=7,units="in")
}

ggplot(pt,aes(x=Yr,y=Catch,color=Alternative)) + geom_line() + geom_point() + mytheme + geom_ribbon(aes(ymin=lb,ymax=ub,fill=config),alpha=0.25) + labs(x="Year")+ scale_x_continuous(breaks=seq(2015,2032,2))  + ylim(c(0,130))
