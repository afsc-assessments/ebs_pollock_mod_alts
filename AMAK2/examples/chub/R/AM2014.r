rm(list=ls())

require(PBSadmb)

library(lattice)
library(RColorBrewer)
library(doBy)

reposDir    <- "/Users/jim/_mymods/amak/examples/atka/"
codePath    <- file.path(reposDir,"R")
inputPath   <- reposDir
outputPath  <- file.path(reposDir,"arc/")
resultPath  <- file.path(reposDir,"results/")
fndir<-"/users/jim/Dropbox/R_common/"
outdir="/Users/jim/_mymods/AMAK/examples/Atka/"
Figdir="/Users/jim/_mymods/AMAK/examples/Atka/figs"
ret_dir="/Users/jim/_mymods/AMAK/examples/Atka/retro"
setwd(codePath)
source(paste(fndir,"adfunctions.r",sep=""))
#source(paste(fndir,"adfunctions2.r",sep=""))
source("diagnostics_v2.r") 
source("ADMB2R_15102013.r")
source("compareRuns.r")

# Read in the output of the assessment
run_name="mod1"
dat_name="am2014.dat"
amak.in  <- read.dat(iFilename = paste(dat_name,sep=""),iPath=inputPath)
amak.out <- readList(file.path(inputPath,paste("arc/",run_name,"_R.rep",sep="")))
amak.ypr <- readYPR(file.path(inputPath,paste( "arc/",run_name,".yld",  sep="")))


mod1 <- readList(file.path(inputPath,paste("arc/mod1_R.rep",sep="")))
mod2 <- readList(file.path(inputPath,paste("arc/mod2_R.rep",sep="")))
mod3 <- readList(file.path(inputPath,paste("arc/mod3_R.rep",sep="")))
mod4 <- readList(file.path(inputPath,paste("arc/mod4_R.rep",sep="")))
# Model selected for 2013
mod2013 <- readList("../arc/mod2013_R.rep")

#diagnostics(mod1,amak.in,amak.ypr,what=c("input","fit","projections"))
#diagnostics(amak.out,amak.in,what=c("input","fit","projections"))
#diagnostics(amak.out,amak.in,amak.ypr,what=c("fit","projections"))
lstOuts   <- list( Model_1= mod1, Model_2= mod2, Model_3= mod3, Model_4= mod4 )  
lstOuts   <- list( Model_1= mod1, Model_4= mod4 )  
lstOuts   <- list( Model_1= mod1, Model_2= mod2, Model_3= mod3)
thislast  <- list( "Model 1 this year"= mod1, "Last year's model"= mod2013)

# used in doc to compare this year w/ last...
compareTime(thislast,"SSB",SD=T,Sum=NULL,legendPos="topleft",startYear=1977,ylim=c(0,500000),main="Spawning biomass")
compareTime(retouts,"SSB",SD=F,Sum=NULL,legendPos="topleft",startYear=1977,ylim=c(0,500000),main="Spawning biomass")
compareTime(retouts,"R",SD=T,legendPos="topright",startYear=1998)
#compareTime(thislast,"TotBiom",legendPos="topright",SD=T)
compareTime(thislast,"R",SD=T,legendPos="topright")
#compareMatrix(thislast,"TotF",SD=F,legendPos="left",Sum=NULL,YrInd=mod1$Yr,Apply=mean,ylim=c(0,0.3))

compareTime(lstOuts,"SSB_NoFishR",SD=T,Sum=NULL,startYear=1977,ylim=c(0,1.0))
compareTime(lstOuts,"SSB",SD=T,Sum=NULL,legendPos="topright",startYear=1977,ylim=c(0,500000))
compareTime(lstOuts,"SSB",SD=T,Sum=NULL,legendPos="topleft",startYear=1970)
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,legendPos="topleft",startYear=1970)
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,legendPos="topleft",startYear=1953)
compareTime(lstOuts,"SSB",SD=T,Sum=NULL,legendPos="topleft",startYear=1970,ylim=c(0,6e5))
compareTime(lstOuts,"R",SD=T,legendPos="topright")
compareTime(lstOuts,"R",legendPos="left",SD=F)
compareTime(lstOuts,"TotBiom",legendPos="topright",SD=T)
compareMatrix(lstOuts,"TotF",SD=F,legendPos="left",Sum=NULL,YrInd=jjm1.2$Yr,Apply=mean,ylim=c(0,0.8))
compareMatrix(lstOuts,"N",   SD=F,legendPos="left",         YrInd=mod1$Yr,Apply=sum)

#pdf(paste(Figdir,"SurvFit.pdf",sep=""),width=9, height=7)
rec_age=1
p.rec.hist(mod1,main="Model 1",ylab="Recruitment",fy=1977,ly=2014,plotmean="T")
# Survey fit
IndexFit(mod1,yf=1989,yl=2014,f=1,main="Model 1", ,ylab="Survey biomass (t)")
p.sur.stk(mod1,f=1)
p.catch.fit(mod1,f=1,ylim=c(0,110000))
plot(mod1$SSB,typ="b",ylab="SSB",xlab="Year",ylim=c(0,.85e6))

#pdf(paste(Figdir,"TotB.pdf",sep=""),width=9, height=7)
p.biom.pol(mod1,typ="SSB",main="Model 1",new=F,fy=1977,ly=2014)
lines(mod1$TotBiom[,1],mod1$TotBiom[,2],col="black",lty=2,lwd=3)
p.biom.pol(mod1,typ="TB",main="Model 1",new=F,fy=1977,ly=2014)
#dev.off()

styr=1977
p.stock.rec(mod1,main="Atka mackerel, model 1")
p.eff.n(mod1,typ="F",f=1,main="Model 1")
p.eff.n(mod1,typ="S",f=1,main="Model 1")
AgeFits(mod1,rec_age=1,case_label="2014 1 assessment")
AgeFitsSrv(mod1,rec_age=1,case_label="2014 1 assessment")
#dev.off()
#detach(dat)
# spawning biomass and last year's estimates 
p.biom.pol(mod1,typ="TB",new=F)
lines(1977:2013,mod1$F_fsh_1[,3],lwd=2)
#names(mod1)

# Rec
p.biom.stk(mod1,typ="R")
# Numbers at age
p.bub.age(mod1,siz=200)
#dev.off()

p.rec.hist(mod2,fy=1977,ly=2014,main="Model 2")
lines(modsigmaR$R[,1],modsigmaR$R[,2],col="purple",lwd=2)
#dev.off()

# Likelihood table
tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
tab
do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)}))
do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)}))
rbind(
cbind("q",do.call(cbind,lapply(lstOuts,function(x){round(x[["q_1"]][1,2],2)}))),
cbind("Npars",do.call(cbind,lapply(lstOuts,function(x){round(x[["Num_parameters_Est"]],2)}))),
cbind("M",do.call(cbind,lapply(lstOuts,function(x){round(x[["Mest"]][2],2)}))),
cbind("SigmaR",do.call(cbind,lapply(lstOuts,function(x){round(x[["Sigmar"]][2],2)}))),
cbind("EffN_Fish",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["EffN_Fsh_1"]][,2]),2)}))),
cbind("EffN_Surv",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["EffN_Survey_1"]][,2]),2)}))),
tab,
cbind("F2014",do.call(cbind,lapply(lstOuts,function(x){round((x[["F_fsh_1"]][38,2]),2)}))),
cbind("F2014/F40%",do.call(cbind,lapply(lstOuts,function(x){round((x[["F_fsh_1"]][38,2])/(x[["F40_est"]]),2)}))),
cbind("B 1977",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][1,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][1,3])/(x[["TotBiom"]][1,2])*100,0)}))),
cbind("B 2014",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][38,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][38,3])/(x[["TotBiom"]][38,2])*100,0)}))),
cbind("2001 YC",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][26,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][26,3])/(x[["R"]][26,2])*100,0)}))),
cbind("2006 YC",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][31,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][31,3])/(x[["R"]][31,2])*100,0)})))
)

?do.call
names(mod1)
write.csv(tab,file=file.path(inputPath,"LikelihoodTable2014.csv"),row.names=F)
system(paste("open ",inputPath,"/Likelihoodtable2014.csv",sep=""))



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


IndexFit(mod2,yf=1989,yl=2012,f=1,main="Model 2",ylab="Survey biomass (t)")
IndexFit(mod1,yf=1989,yl=2012,f=1,main="Model 1",ylab="Survey biomass (t)")
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
Plot_Fspr()
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
spwn_ratio(mod1,fy=1977,ly=2013) 
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

#spwn_ratio(am1) 
#plt_proj(am1)
#dev.off()

#=====================
# Retro
# get all retrospectives
#=====================
setwd(outdir)
outdir
retouts <- list()
retouts <- list( R0=mod1, R1= retro1, R2= retro2, 
  R3= retro3, 
  R4= retro4 ,
  R5= retro5 ,
  R6= retro6 ,
  R7= retro7 ,
  R8= retro8 ,
  R9= retro9 ,
  R10= retro10 
  )  
tab = list(R0=mod1$SSB)

for (i in 1:10) { tab=cbind(tab,retouts[[i]]$SSB)}

for (i in 0:15) {
  rn=paste("retro/r_",i,".rep",sep="")
  mn=paste("retro",i,sep="")
  assign(mn,readList(rn))
  print(rn)
}
#pdf(paste(Figdir,"Retro_Mods.pdf",sep=""),width=9, height=6)
p.biom.pol(retro0,typ="SSB",main="Model 1",new=F,fy=1977,ly=2014)
for (i in 1:15) {
  rn=paste("retro",i,sep="");
  lines(get(rn)$SSB[,1],get(rn)$SSB[,2],col=i)
  lr=length(get(rn)$SSB[,1])
  points(get(rn)$SSB[lr,1],get(rn)$SSB[lr,2],pch=19,col=i)
}

plot(retro0$SSB[,1],rep(0,length(retro0$SSB[,1])), ylim=c(-.7,.7), ylab="Relative difference from terminal year",
type="l",xlab="Year",lty=2,lwd=2)
 for (i in 1:15) {
  rn=paste("retro",i,sep="");
  lr=length(get(rn)$SSB[,1])
  lines(get(rn)$SSB[,1],get(rn)$SSB[,2]/retro0$SSB[1:lr,2]-1,col=i)
  #points(get(rn)$SSB[lr,1],get(rn)$SSB[lr,2],pch=19,col=i)
}    
# Mohn's rho
rc = retro0$SSB[,2]
ntmp=0
rho=0
for (i in 1:15) {
  dtmp=get(paste("retro",i,sep=""))$SSB
  lr=length(dtmp[,1])
  ntmp= ntmp+(dtmp[lr,2] -rc[lr])/rc[lr]
  #rho = rho + (-(ALL[i,2]-ALL[*tsyrs-i,2+i]))/ALL[(j)*tsyrs-i,2]
  rho = rho + (-(dtmp[i,2]-rc[i]))/rc[lr]
  print(paste(i,ntmp/i,rho))
}    
       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
  #points(get(rn)$SSB[lr,1],get(rn)$SSB[lr,2],pch=19,col=i)
  #print(dtmp)
  #print(paste(dtmp[,1],rc[lr]))
 for(i in 1:15) {

dev.off()

p.biom.pol(retro2,typ="SSB",main="Model 1",new=F,fy=1977,ly=2013)
for (i in 1:9) {
  rn=paste("retro1",i,sep="");
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