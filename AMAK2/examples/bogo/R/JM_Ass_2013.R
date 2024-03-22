#-------------------------------------------------------------------------------
# Script to run the assessment of Jack Mackerel and look at outputs
#
# By: Niels Hintzen
# Upate: 31 Aug 2011
#-------------------------------------------------------------------------------
rm(list=ls())
memory.size(4000)

# Set libraries & source code
library(lattice)
require(PBSadmb)
library(RColorBrewer)
library(doBy)

#-------------------------------------------------------------------------------
# Set paths
#-------------------------------------------------------------------------------
reposDir    <- "D:/Repository/JackMackerel/"
reposDir    <- "c:/users/jim/documents/_mymods/jjm/"
codePath    <- file.path(reposDir,"Code/R/")
inputPath   <- file.path(reposDir,"Code/admb/")
outputPath  <- file.path(reposDir,"Code/admb/arc/")
resultPath  <- file.path(reposDir,"Results/Assessment/")
setwd(codePath)

# Specify control file
controlFile <- "mod4.1.ctrl"

# Run the assessment
source("diagnostics_v2.r")
source("ADMB2R_15102013.r")
source("compareRuns.r")
system(paste('"jjm.exe"','-ind',paste(controlFile,".ctr",sep=""),'-nox'), wait = TRUE)

# Read in the output of the assessment
run_name="mod1.4"
run_name="mod1.2"
dat_name="mod0.5.dat"
run_name="mod1.4"
run_name="mod1.5"
run_name="mod1.7"
run_name="mod3.1"
run_name="mod1.4"
run_name="mod3.2"
run_name="mod4.1"
run_name="mod4.2"
run_name="mod4.3"
run_name="mod4.4"



dat_name="mod3.dat"

jjm.in  <- read.dat(iFilename = paste(dat_name,sep=""),iPath=inputPath)
jjm.out <- readList(file.path(inputPath,paste("arc/",run_name,"_r.rep",sep="")))
jjm.ypr <- readYPR(file.path(inputPath,paste( "arc/",run_name,".yld",  sep="")))

#-------------------------------------------------------------------------------
# Create diagnostics
#-------------------------------------------------------------------------------
pdf(paste(resultPath,"summary_",run_name,".pdf",sep=""),height=29.7/2.54,width=21/2.54,pointsize = 16, bg = "white")
diagnostics(jjm.out,jjm.in,jjm.ypr,what=c("input","fit","projections","ypr"))
diagnostics(jjm.out,jjm.in,jjm.ypr,what=c("fit","projections"))
dev.off()
#Landscape
pdf(paste(resultPath,"summary_LS_",run_name,".pdf",sep=""),width=29.7/2.54,height=21/2.54,pointsize = 16, bg = "white")
diagnostics(jjm.out,jjm.in,jjm.ypr,what=c("input","fit","projections","ypr"))
diagnostics(jjm.out,jjm.in,jjm.ypr,what=c("fit","projections"))
dev.off()

#Write output to file
writeList(setOutputNames(jjm.out),fname=paste(controlFile,"_out.txt",sep=""),format="P")


#-------------------------------------------------------------------------------
# Visual compare runs
#-------------------------------------------------------------------------------
source("compareRuns.r")
getwd()
jjm0.0 <- readList(file.path(inputPath,paste("arc/Mod0.0_r.rep",sep="")))
jjm0.1 <- readList(file.path(inputPath,paste("arc/Mod0.1_r.rep",sep="")))
jjm0.2 <- readList(file.path(inputPath,paste("arc/Mod0.2_r.rep",sep="")))
jjm0.3 <- readList(file.path(inputPath,paste("arc/Mod0.3_r.rep",sep="")))
jjm0.4 <- readList(file.path(inputPath,paste("arc/Mod0.4_r.rep",sep="")))
jjm0.5 <- readList(file.path(inputPath,paste("arc/Mod0.5_r.rep",sep="")))
jjm0.6 <- readList(file.path(inputPath,paste("arc/Mod0.6_r.rep",sep="")))
jjm1.1 <- readList(file.path(inputPath,paste("arc/Mod1.1_r.rep",sep="")))
jjm1.2 <- readList(file.path(inputPath,paste("arc/Mod1.2_r.rep",sep="")))
jjm1.3 <- readList(file.path(inputPath,paste("arc/Mod1.3_r.rep",sep="")))
jjm1.4 <- readList(file.path(inputPath,paste("arc/Mod1.4_r.rep",sep="")))
jjm1.5 <- readList(file.path(inputPath,paste("arc/Mod1.5_r.rep",sep="")))
jjm1.6 <- readList(file.path(inputPath,paste("arc/Mod1.6_r.rep",sep="")))
jjm1.7 <- readList(file.path(inputPath,paste("arc/Mod1.7_r.rep",sep="")))
jjm1.8 <- readList(file.path(inputPath,paste("arc/Mod1.8_r.rep",sep="")))
jjm1.9 <- readList(file.path(inputPath,paste("arc/Mod1.9_r.rep",sep="")))
jjm1.10<- readList(file.path(inputPath,paste("arc/Mod1.10_r.rep",sep="")))
jjm3.1<- readList(file.path(inputPath,paste("arc/Mod3.1_r.rep",sep="")))
jjm3.2<- readList(file.path(inputPath,paste("arc/Mod3.2_r.rep",sep="")))
jjm4.1<- readList(file.path(inputPath,paste("arc/Mod3.1_r.rep",sep="")))
jjm4.2<- readList(file.path(inputPath,paste("arc/Mod3.2_r.rep",sep="")))


lstOuts <- list()
for (i in 1:9){
  ttt<-  paste(inputPath,"retro/r",i-1,"_r.rep",sep="")
  lstOuts[[i]] <- readList(ttt)
  names(lstOuts)[i] <- paste("retro_",i-1,sep="")
}

compareTime(lstOuts,"SSB",SD=T,Sum=NULL,startYear=1970,legendPos="left")
compareTime(lstOuts,"R",SD=T,Sum=NULL,startYear=1970)

lstOuts   <- list(
  Model_1.4= jjm1.4,
  Model_1.3= jjm1.3,
  Model_1.2= jjm1.2
)  
lstOuts   <- list(
  Model_3.1= jjm3.1,
  Model_1.4= jjm1.4,
  Model_3.2= jjm3.2
)  

lstOuts   <- list(
  Model_0.4= jjm0.4,
  Model_1.1= jjm1.1
)  
lstOuts   <- list(
  Model_1.6= jjm1.6,
  Model_1.2= jjm1.2
)  
lstOuts   <- list(
  Model_1.6= jjm1.6,
  Model_1.7= jjm1.7,
  Model_1.2= jjm1.2
)  
lstOuts   <- list(
  Model_1.9= jjm1.9,
  Model_1.2= jjm1.2
)  

lstOuts   <- list(
  Model_0.0= jjm0.0,
  Model_0.1= jjm0.1,
  Model_0.2= jjm0.2,
  Model_0.3= jjm0.3,
  Model_0.4= jjm0.4
  Model_1.1= jjm1.1,
  Model_1.2= jjm1.2,
  Model_1.3= jjm1.3,
  Model_1.4= jjm1.4,
  Model_1.5= jjm1.5,
  Model_1.6= jjm1.6,
  Model_1.7= jjm1.7,
  Model_1.8= jjm1.8,
  Model_1.9= jjm1.9
                    )
lstOuts   <- list(
  Model_0.5= jjm0.5,
  Model_0.6= jjm0.6
)
lstOuts   <- list(
  Model_1.2= jjm1.2,
  Model_0.6= jjm0.6
)

pdf(paste(resultPath,"Compare_0.pdf",sep=""),width=29.7/2.54,height=21/2.54,pointsize = 24, bg = "white")
pdf(paste(resultPath,"Compare_1.pdf",sep=""),width=29.7/2.54,height=21/2.54,pointsize = 24, bg = "white")

pdf(paste(resultPath,"Compare1_2_w_1_4_",".pdf",sep=""),height=7,width=10,pointsize = 16, bg = "white")
pdf(paste(resultPath,"Compare1_4_w_3_",".pdf",sep=""),height=7,width=10,pointsize = 16, bg = "white")
pdf(paste(resultPath,"Compare1_2_w_1_9",".pdf",sep=""),height=7,width=10,pointsize = 16, bg = "white")

compareTime(lstOuts,"SSB_NoFishR",SD=T,Sum=NULL,startYear=1990,ylim=c(0,.4))
compareTime(lstOuts,"SSB",SD=T,Sum=NULL,legendPos="left",startYear=1990,ylim=c(0,7000))
compareTime(lstOuts,"SSB",SD=T,Sum=NULL,legendPos="left",startYear=1970)
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,legendPos="left",startYear=1970)
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,legendPos="left",startYear=1953)
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,legendPos="left",startYear=1970)
compareTime(lstOuts,"R",SD=T,legendPos="left")
compareTime(lstOuts,"R",legendPos="left",SD=F)
compareTime(lstOuts,"TotBiom",legendPos="left",SD=T)
compareMatrix(lstOuts,"TotF",SD=F,legendPos="left",Sum=NULL,YrInd=jjm1.2$Yr,Apply=mean,ylim=c(0,0.8))
compareMatrix(lstOuts,"N",   SD=F,legendPos="left",         YrInd=jjm1.2$Yr,Apply=sum)
dev.off()
jjm.n2    <- readList(file.path(outputPath,"n3_r.rep"))
lstOuts   <- list(Model_N1=jjm.n1,Model_N2=jjm.n2)
jjm.s1    <- readList(file.path(outputPath,"s1_r.rep"))
jjm.s2    <- readList(file.path(outputPath,"s2_r.rep"))
lstOuts   <- list(Model_S1=jjm.s1,Model_S2=jjm.s2)


pdf(paste(resultPath,"Compare_0.pdf",sep=""),height=29.7/2.54,width=21/2.54,pointsize = 16, bg = "white")
#-------------------------------------------------------------------------------
# Numerical compare runs
#-------------------------------------------------------------------------------
# Likelihood table
tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
write.csv(tab,file=file.path(inputPath,"LikelihoodTable2013.csv"),row.names=F)
system(paste(inputPath,"/Likelihoodtable2013.csv",sep=""))
getwd()
# Risk tables
jjm.1.2  <- readList(file.path(inputPath,"mod1.1_r.rep"))gi
jjm.1.8  <- readList(file.path(inputPath,"mod1.8_r.rep"))
lstOuts   <- list(Model_1.2=jjm.1.2,Model_1.8=jjm.1.8)

# Get the future SSBs and SDs together in one file
fut       <- do.call(rbind,lapply(lstOuts,function(x){
                     do.call(rbind,lapply(x[grep("SSB_fut_",names(x))],
                                          function(y){return(y[,1:3])}))}))
fut       <- as.data.frame(fut,stringsAsFactors=F)
colnames(fut) <- c("year","SSB","SD")
fut$modelscenario <- paste(rep(names(lstOuts),each=nrow(lstOuts[[1]]$SSB_fut_1) *
                                                   length(grep("SSB_fut_",names(lstOuts[[1]])))),
                           paste("Scen",
                                 rep(1:length(grep("SSB_fut_",names(lstOuts[[1]]))),each=nrow(lstOuts[[1]]$SSB_fut_1)),
                                 sep="_"),
                           sep="_")
#Get the 2012 SSB and SDs together
fut
ass       <- do.call(rbind,lapply(lstOuts,function(x){
                                  x$SSB[which(x$SSB[,1] == (x$SSB_fut_1[1,1]-1)),1:3]}))
ass       <- as.data.frame(ass,stringsAsFactors=F)
colnames(ass) <- c("year","SSB","SD")
ass$modelscenario <- names(lstOuts)

#Do the risk calculation:
# Each final year, what is the risk of SSB(final year) < SSB(ass year)
rsktable    <- matrix(NA,nrow=length(lstOuts),ncol=length(grep("SSB_fut_",names(lstOuts[[1]]))),
                      dimnames=list(names(lstOuts),1:length(grep("SSB_fut_",names(lstOuts[[1]])))))
ratiotable  <- matrix(NA,nrow=length(lstOuts),ncol=length(grep("SSB_fut_",names(lstOuts[[1]]))),
                      dimnames=list(names(lstOuts),1:length(grep("SSB_fut_",names(lstOuts[[1]])))))

for(i in names(lstOuts)){
  futdat <-subset(fut,year==max(fut$year) &
                  paste("Model_",unlist(strsplit(fut$modelscenario,"_"))[seq(2,nrow(fut)*4,4)],sep="") == i)
  assdat <- subset(ass,modelscenario == i)
  rsktable[i,] <- dnorm(assdat$SSB,futdat$SSB,futdat$SD)
  ratiotable[i,] <- round((futdat$SSB / assdat$SSB),2 )# / futdat$SSB *100,1)
}
rsktable
ratiotable

#Average catch over entire timeseries
meancatch <- do.call(rbind,lapply(lstOuts,function(x){unlist(lapply(x[grep("Catch_fut_",names(x))],function(y){mean(y[,2])}))}))
meancatch
