# ------------------------------------------------------------------------
# Demo script ------------------------------------------------------------
# jjmTools: Graphics and diagnostics libraries for SPRFMO's JJM mo -------
# ------------------------------------------------------------------------

library(devtools)
install_github(repo="imarpe/jjmTools")
library(jjmTools)
rm(list=ls())
# Set parameters ----------------------------------------------------------
# Path of JJM repository (from current working directory)
setwd("/Users/jim/_mymods/jjm/Code/R")
source("/Users/Jim/_mymods/jjmTools/Code/R/*")
.LIB      = "/Users/Jim/_mymods/jjmTools/R/"
.RFILES   = list.files(.LIB,pattern="\\.[Rr]$")
for(nm in .RFILES) source(file.path(.LIB, nm), echo=FALSE)
reposDir =  "../../../AMAK/examples/chub"
# Name models in a list
Mod2List = c("mod2","mod1" )
Mod0List = paste0("mod0.", 0:1)
Mod1List = paste0("mod1.", c(0:5,9,10,11))
Mod1Listn = paste0("mod1.", c(0:5,9,10,11),"n")
Mod2List
Mod1Listn
# Last year's model....
Mod2013 = "mod4.1"
# Run models --------------------------------------------------------------
#runJJM(modelName = Mod0List, path = reposDir)
#runJJM(modelName = Mod1List, path = reposDir)
mod2 = readJJM(model  = Mod0List[1], path = reposDir)
mod2 = readJJM(model  = Mod2List[1], path = reposDir)
mod2 = readJJM(model  = "mod2",    path = reposDir)
plot(mod2)
# One stop shopping for all annex figures
source("annexfigs.r")
kobe(mod2)
kobe(mod4.1,add=T,col="grey")
kobe(mod2.2,add=T,col="grey",lwd=2)
kobe(mod2.1,add=T,col="grey",lwd=2)
kobe(mod2.1)
kobe(mod2.3)
dp2.0=diagnostics(mod2.0)
dp2.1=diagnostics(mod2.1)
dp2.2=diagnostics(mod2.2)
dp2.3=diagnostics(mod2.3)
par(oma=rep(2,4))
plot(dp2.3, what="projections",var="ssbPrediction",ylab="SSB (kt)",main="Model 2.3", ylim=c(0, 15000), xlim=c(2000, 2025))
legend(2020,15000,1:3)
plot(dp2.0, what = "fit",var="summarySheet")
plot(dp2.1, what="projections",var="ssbPrediction",ylab="SSB (kt)",main="Model 2.1")
plot(dp2.2, what="projections",var="ssbPrediction",ylab="SSB (kt)",main="Model 2.2",ylim=c(0,15000))
plot(dp2.3, what="projections",var="ssbPrediction",ylab="SSB (kt)",main="Model 2.3",ylim=c(0,15000))

with(mod2.3$output$output,msy_mt)
	{ cbind( Catch_fut_1[,1:2], Catch_fut_2[,2], Catch_fut_3[,2], Catch_fut_4[,2])})
names(mod2.0$output$output)
cfut0 = with(mod2.0$output$output,{ cbind( Catch_fut_1[,1:2], Catch_fut_2[,2], Catch_fut_3[,2], Catch_fut_4[,2])})
cfut1 = with(mod2.1$output$output,{ cbind( Catch_fut_1[,1:2], Catch_fut_2[,2], Catch_fut_3[,2], Catch_fut_4[,2])})
cfut2 = with(mod2.2$output$output,{ cbind( Catch_fut_1[,1:2], Catch_fut_2[,2], Catch_fut_3[,2], Catch_fut_4[,2])})
cfut3 = with(mod2.3$output$output,{ cbind( Catch_fut_1[,1:2], Catch_fut_2[,2], Catch_fut_3[,2], Catch_fut_4[,2])})
rbind( cfut0[1:2,], cfut1[1:2,], cfut2[1:2,], cfut3[1:2,] )

cfut = with(mod2.1$output$output,{ cbind( Catch_fut_1[,1:2], Catch_fut_2[,2], Catch_fut_3[,2], Catch_fut_4[,2])})
cfut = with(mod2.2$output$output,{ cbind( Catch_fut_1[,1:2], Catch_fut_2[,2], Catch_fut_3[,2], Catch_fut_4[,2])})
with(mod2.3$output$output,{ cbind( Catch_fut_1[,1:2], Catch_fut_2[,2], Catch_fut_3[,2], Catch_fut_4[,2])})
ssbfut0= with(mod2.0$output$output,{ cbind(  SSB_fut_1[,1:2],  SSB_fut_2[,2],  SSB_fut_3[,2],  SSB_fut_4[,2],SSB_fut_5[,2])})
ssbfut2= with(mod2.2$output$output,{ cbind(  SSB_fut_1[,1:2],  SSB_fut_2[,2],  SSB_fut_3[,2],  SSB_fut_4[,2],SSB_fut_5[,2])})
ssbfut= with(mod2.1$output$output,{ cbind(  SSB_fut_1[,1:2],  SSB_fut_2[,2],  SSB_fut_3[,2],  SSB_fut_4[,2],SSB_fut_5[,2])})
ssbfut= with(mod2.2$output$output,{ cbind(  SSB_fut_1[,1:2],  SSB_fut_2[,2],  SSB_fut_3[,2],  SSB_fut_4[,2],SSB_fut_5[,2])})
with(mod2.3$output$output,{ cbind(  SSB_fut_1[,1:2],  SSB_fut_2[,2],  SSB_fut_3[,2],  SSB_fut_4[,2],SSB_fut_5[,2])})
cfut0

par(mfrow=c(1,1))
plot(cfut[,2],ssbfut[,2],ylim=c(0,8000),xlim=c(0,1500),ylab="SSB (kt)",xlab="Catch")
text(cfut[,2],ssbfut[,2],cfut[,1])
text(cfut2[,2],ssbfut2[,2],cfut[,1],col="red")
text(cfut2[,3],ssbfut2[,3],cfut[,1],col="red")
text(cfut2[,4],ssbfut2[,4],cfut[,1],col="red")
text(cfut2[,5],ssbfut2[,5],cfut[,1],col="red")
text(cfut2[,6],ssbfut2[,6],cfut[,1],col="red")
text(cfut0[,2],ssbfut0[,2],cfut[,1],col="red")

par(mfrow=c(2,1))
plot(ssbfut[,1],ssbfut[,2],ylim=c(0,8000),ylab="SSB (kt)",xlab="Year")
lines(ssbfut[,1],ssbfut[,3])
lines(ssbfut[,1],ssbfut[,4])
lines(ssbfut[,1],ssbfut[,5])
lines(ssbfut[,1],ssbfut[,6])
plot(cfut[,1],cfut[,2],ylim=c(0,1000),ylab="Catch (kt)",xlab="Year")
lines(cfut[,1],cfut[,3])
lines(cfut[,1],cfut[,4])
lines(cfut[,1],cfut[,5])

,ylim=c(0,1000),ylab="Catch (kt)",xlab="Year")
with(mod2.0$output$output,{ cbind(  SSB_fut_1[,1:2],  SSB_fut_2[,2],  SSB_fut_3[,2],  SSB_fut_4[,2],SSB_fut_5[,2])})

mod2.1$output$output$Catch_fut_1
mod2.1$output$output$Catch_fut_1
mod2.1$output$output$Catch_fut_1

mod1234 = combineModels(mod1.0,mod1.1, mod1.2, mod1.3, mod1.4,mod1.5,mod1.9,mod1.10,mod1.11)
mod01 = combineModels(mod1.0,mod1.1)
mod02 = combineModels(mod1.0,mod1.2)
mod03 = combineModels(mod1.0,mod1.3)
mod04 = combineModels(mod1.0,mod1.4)
mod05 = combineModels(mod1.0,mod1.5)
mod09 = combineModels(mod1.0,mod1.9)
mod010= combineModels(mod1.0,mod1.10)
mod011= combineModels(mod1.0,mod1.11)
mod02013 = combineModels(mod1.0,mod4.1)
mod002013 = combineModels(mod0.0,mod4.1)
mod2.3.2013 = combineModels(mod2.3,mod4.1)
mod2.0.2013 = combineModels(mod2.0,mod4.1)
mod0001 = combineModels(mod0.0,mod0.1)
mod0002 = combineModels(mod0.0,mod0.2)
plot(mod2.0.2013)
plot(mod0001)
plot(mod0002)
logLik(mod1234,details=T)
summary(mod1234)$lik
plot(mod1234)
plot(mod02013)
plot(mod002013)
plot(mod01)
plot(mod02)
plot(mod03)
plot(mod05)
plot(mod09)
plot(mod010)
plot(mod011)
names(mod1234$combined$outputs)
,file=file.path(reposDir,"likTable.csv")
logLik(mod1234)
#Survey lik
mods  <- mod1234
surfLik   <- do.call(cbind,lapply(mods$combined$outputs,function(x){do.call(rbind,x[grep("Survey_Index_",names(x))])}))
colnames(surfLik) <- mods$info

#Age survey lik
mods  <- mod123
ageSurvLik   <- do.call(cbind,lapply(mods$combined$outputs,function(x){do.call(rbind,x[grep("Age_Survey_",names(x))])}))
colnames(ageSurvLik) <- mods$info

rbind(surfLik,ageSurvLik)
Print_Figs <- function(Mod,outname="Modx.pdf"){
	pdf(outname,height=9,width=7)
  diagPlots = diagnostics(outputObject = Mod)
  plot(diagPlots, what = "input")
  plot(diagPlots, what = "fit")
  plot(diagPlots, what = "projections")
  dev.off()
}
Print_Figs(mod1n.11,"../figs/Mod1.11n.pdf")
Print_Figs(mod1.1,"../figs/Mod1.1.pdf")
Print_Figs(mod1.0,"../figs/Mod1.0.pdf")
Print_Figs(mod1.4,"../figs/Mod1.4.pdf")
Print_Figs(mod1.2,"../figs/Mod1.2.pdf")
Print_Figs(mod1.3,"../figs/Mod1.3.pdf")
Print_Figs(mod1.5,"../figs/Mod1.5.pdf")
Print_Figs(mod1.10,"../figs/Mod1.10.pdf")

Print_Figs(mod2.2,"../figs/Mod2.2.pdf")
Print_Figs(mod2.3,"../figs/Mod2.3.pdf")

# Reading -----------------------------------------------------------------
# OUTPUT Object
model = readJJM(model= modelName, path = reposDir)

# DIAG object
diagPlots = diagnostics(outputObject = model)

# Combine models ----------------------------------------------------------
plot(mod1234)
plot(mod1234,xlim=c(2000,2015))
plot(mod1234,xlim=c(2010, 2015), ylim=c(2,8), legendPos = "topleft")
mod012 = combineModels(mod0,mod1, mod2kod3)
mod012 = combineModels(mod0,mod1, mod2kod3)
plot(mod012)

# Integrating models ------------------------------------------------------
mod12 = combineStocks(mod1, mod2, model = "mod2s_12")

# Print -------------------------------------------------------------------

# Output object
print(model)

# List of outputs 
print(mod1234)

# Diagnostics object
print(diagPlots)

# Get and print summaries -------------------------------------------------
# Output object
sumModel = summary(model)
print(sumModel)

# List of outputs object
sumList = summary(mod012)
print(sumList)
# Diagnostics object
sumPlots = summary(diagPlots)
sumPlots


# Get and print plots -----------------------------------------------------
plot(diagPlots, what = "input")
plot(diagPlots, what = "fit")
pdf()
plot(dp, what = "fit")
plot(diagPlots, what = "projections")
plot(diagPlots, what = "ypr")


# Run models in parallel---------------------------------------------------
Mod0List = paste0("mod0.", 0:2)
nCores = 6
cl = makeCluster(rep("localhost",times=nCores), type = "SOCK")
registerDoSNOW(cl)
runJJM(models = Mod0List, path = reposDir, parallel=TRUE)
stopCluster(cl)

