#-------------------------------------------------------------------------------
# Script to run the assessment of Jack Mackerel and look at outputs
#-------------------------------------------------------------------------------

setwd("N:/Projecten/SouthPacific/2011/Assessment2/")

rm(list=ls())
memory.size(4000)

  # Set libraries & source code
library(lattice)
require(PBSadmb)
library(RColorBrewer)
source("R/ADMB2R.r")

  # Define population characteristics
popSettings     <- scan("mod1_Proj.dat",what='numeric', quiet=TRUE,sep="\n",comment.char="#",allowEscapes=T)
popSettingsMods <- scan("Mod_Schedules_ForNiels.txt",what='numeric',quiet=TRUE,sep="\n",comment.char="#",allowEscapes=T)
recrYrs         <- 5

ages            <- 2:12
natmort         <- na.omit(an(unlist(strsplit(popSettings[13]," "))))
stockwt         <- na.omit(an(unlist(strsplit(popSettings[15]," "))))
mat             <- na.omit(an(unlist(strsplit(popSettings[14]," "))))/2
spwn            <- na.omit(an(popSettings[10]))/12
landwt          <- na.omit(an(unlist(strsplit(popSettings[16]," "))))
fpattern        <- na.omit(an(unlist(strsplit(popSettings[17]," "))))

fpattern        <- matrix(c(c(0.047435075,0.210153162,0.460193999,0.722022158,0.813072429,0.771646581,0.947260551,1,0.928219578,0.928219578,0.928219578),
                            c(0.063793696,0.283687898,0.57459911,0.907039853,1,0.896024784,0.935015955,0.994236322,0.915103449,0.915103449,0.915103449),
                            c(0.046168417,0.208166215,0.458854538,0.732492445,0.848439077,0.806886618,0.947484642,1,0.953838298,0.953838298,0.953838298)),
                          nrow=3,ncol=length(ages),dimnames=list(model=1:3,ages=ages),byrow=T)
N2011           <- matrix(c(c(12488800,4555270,3626240,1026310,398459,215514,290787,166396,   70561,  55213.4,   60371.1),
                            c(15662100,5933790,5588200,1760890,755634,396340,520966,331232,  163681,   136339,  158328),
                            c(8632520, 3152440,2445910,706608, 272997,136430,167424,93548.2,40340.2,  32978.3,  36713.5)),
                          nrow=3,ncol=length(ages),dimnames=list(model=1:3,ages=ages),byrow=T)
recr            <- matrix(c(c(2786050,2423700,3501710,7428160,6043840),
                            c(3440060,3499960,5000230,10419100,7721940),
                            c(2588660,2159060,2891570,5532510,4275660)),
                          nrow=3,ncol=5,dimnames=list(model=1:3,years=2006:2010),byrow=T)
  # Get projection settings
projSettings    <- scan("mod1.prj",what='numeric', quiet=TRUE,sep="\n",comment.char="#",allowEscapes=T)
projYrs         <- an(projSettings[7])
simNum          <- an(projSettings[8])
scnNum          <- 5
ImYCatch        <- an(projSettings[13])



  # Store variables
SSB             <- array(NA,dim=c(1,            projYrs,simNum,scnNum,3),dimnames=list(ages=1,   years=1:projYrs,iter=1:simNum,scenario=1:scnNum,model=1:3))
F               <- array(NA,dim=c(length(ages), projYrs,simNum,scnNum,3),dimnames=list(ages=ages,years=1:projYrs,iter=1:simNum,scenario=1:scnNum,model=1:3))
Fmult           <- array(NA,dim=c(1,            projYrs,simNum,scnNum,3),dimnames=list(ages=1,   years=1:projYrs,iter=1:simNum,scenario=1:scnNum,model=1:3))
yield           <- array(NA,dim=c(length(ages), projYrs,simNum,scnNum,3),dimnames=list(ages=ages,years=1:projYrs,iter=1:simNum,scenario=1:scnNum,model=1:3))
N               <- array(NA,dim=c(length(ages), projYrs,simNum,scnNum,3),dimnames=list(ages=ages,years=1:projYrs,iter=1:simNum,scenario=1:scnNum,model=1:3))

  # Random structures
for(iMod in 1:3){
  meanR         <- mean(recr[iMod,])
  HmeanR        <- 1 / mean( 1 / recr[iMod,])
  gamm            <- meanR / HmeanR
  gi_beta         <- meanR
  delta           <- 1 / (gamm - 1)
  cvrec           <- sqrt(1 / delta)
  psi             <- rnorm(projYrs*simNum*scnNum,0,1)^2
  omega           <- gi_beta * (1 + (psi -sqrt(4 * delta * psi + psi^2)) / (2 * delta))
  zeta            <- gi_beta * (1 + (psi +sqrt(4 * delta * psi + psi^2)) / (2 * delta))
  gtheta          <- gi_beta / (gi_beta + omega)

    runifs        <- runif(projYrs*simNum*scnNum,0,1)
  N[1,,,,iMod][which(runifs <= gtheta)]     <- omega[which(runifs <= gtheta)]
  N[1,,,,iMod][which(runifs >  gtheta)]     <- zeta[which(runifs >  gtheta)]
}
for(iMod in 1:3) N[,1,,,iMod]               <- rep(N2011[iMod,],simNum*scnNum)

  # Functions needed for projection calculation
calcCatch   <-  function(x,f,m,n,fwt,tac){
                  catch <- sum(f*x / (f*x + m) * n * (1 - exp(-(f*x+m))) *fwt)
                return(sqrt((tac-catch)^2))}

# Simulate projections
for(iYr in 1:(projYrs-1)){
  for(iScn in 1:scnNum){
    if(iScn == 1) target <- 520000
    if(iScn == 2) target <- 390000
    if(iScn == 3) target <- 260000
    if(iScn == 4) target <- 130000
    if(iScn == 5) target <-   5000
    if(iYr == 1)  target <-         ImYCatch
    for(iTer in 1:simNum){
      for(iMod in 1:3){
        Fmult[,iYr,iTer,iScn,iMod] <- optim(1,fn=calcCatch,lower=0,upper=3,method="L-BFGS-B",f=fpattern[iMod,],m=natmort,n=N[,iYr,iTer,iScn,iMod],fwt=landwt,tac=target)$par
      }
    }
  }
  for(iMod in 1:3){
    F[,iYr,,,iMod]              <- outer(fpattern[iMod,],Fmult[,iYr,,,iMod],"*")
  }
  survivors                     <- N[,iYr,,,] * exp(-sweep(F[,iYr,,,],1,natmort,"+"))
  N[2:length(ages),iYr+1,,,]    <- survivors[-dim(survivors)[1],,,]
  # Plusgroup
  N[length(ages),iYr+1,,,]      <- N[length(ages),iYr+1,,,] + survivors[dim(survivors)[1],,,]
  SSB[,iYr,,,]                  <- apply(sweep(N[,iYr,,,],1,stockwt * mat,"*") * exp(-sweep(F[,iYr,,,]*spwn,1,natmort*spwn,"+")),2:4,sum)
}
save.image("projections.RData")



SSBmed  <- apply(SSB,  c(2,4:5),median,na.rm=T)
Fmean   <- apply(F,    c(2:5),  mean,na.rm=T)
Fmed    <- apply(Fmean,c(1,3:4),median,na.rm=T)

projRes <- as.data.frame(cbind(year=rep(2011:2024,5*3),F=Fmed,SSB=SSBmed,scenario=rep(rep(1:5,each=14),3),model=rep(1:3,each=14*5)))


#- Figures
  # Legend to figures
  ikey        <- simpleKey(text=c("99%","75%","50%","25%","0%"),
                           points=T,lines=T,columns = 3)
  ikey$points$pch <- c(3,8,15,16,17)
  ikey$points$col <- c("darkblue","red","green","purple","cyan")
  ikey$lines$lwd  <- 2
  ikey$lines$col  <- c("darkblue","red","green","purple","cyan")

Fmedrange <- c(0,0.4)#c(0,range(Fmed,na.rm=T)[2])
xyplot(F~year,data=subset(projRes,year %in% 2011:2021 & model == 1),
       groups=scenario,xlab="Years",ylab="Fishing mortality",main="Model 1",
       type="b",col=c("darkblue","red","green","purple","cyan"),
       pch=c(3,8,15,16,17),
       key=ikey,prepanel=function(...) {list(ylim=c(0,Fmedrange))})
xyplot(F~year,data=subset(projRes,year %in% 2011:2021 & model == 2),
       groups=scenario,xlab="Years",ylab="Fishing mortality",main="Model 2",
       type="b",col=c("darkblue","red","green","purple","cyan"),
       pch=c(3,8,15,16,17),
       key=ikey,prepanel=function(...) {list(ylim=c(0,Fmedrange))})
xyplot(F~year,data=subset(projRes,year %in% 2011:2021 & model == 3),
       groups=scenario,xlab="Years",ylab="Fishing mortality",main="Model 3",
       type="b",col=c("darkblue","red","green","purple","cyan"),
       pch=c(3,8,15,16,17),
       key=ikey,prepanel=function(...) {list(ylim=c(0,Fmedrange))})
savePlot("./Results/Model1_F.png",type="png")

SSBmedrange <- c(0,8000000)#c(0,range(SSBmed,na.rm=T)[2])
xyplot(SSB~year,data=subset(projRes,year %in% 2011:2021 & model == 1),
       groups=scenario,xlab="Years",ylab="Spawning Stock Biomass",main="Model 1",
       type="b",col=c("darkblue","red","green","purple","cyan"),
       pch=c(3,8,15,16,17),
       key=ikey,prepanel=function(...) {list(ylim=c(0,SSBmedrange))})
xyplot(SSB~year,data=subset(projRes,year %in% 2011:2021 & model == 2),
       groups=scenario,xlab="Years",ylab="Spawning Stock Biomass",main="Model 2",
       type="b",col=c("darkblue","red","green","purple","cyan"),
       pch=c(3,8,15,16,17),
       key=ikey,prepanel=function(...) {list(ylim=c(0,SSBmedrange))})
xyplot(SSB~year,data=subset(projRes,year %in% 2011:2021 & model == 3),
       groups=scenario,xlab="Years",ylab="Spawning Stock Biomass",main="Model 3",
       type="b",col=c("darkblue","red","green","purple","cyan"),
       pch=c(3,8,15,16,17),
       key=ikey,prepanel=function(...) {list(ylim=c(0,SSBmedrange))})

savePlot("./Results/Model1_B.png",type="png")