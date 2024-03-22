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
controlFile <- "mod6a.ctrl"

# Run the assessment
source("diagnostics_v2.r")
source("ADMB2R_15102013.r")
system(paste('"jjm.exe"','-ind',paste(controlFile,".ctr",sep=""),'-nox'), wait = TRUE)

# Read in the output of the assessment
run_name="mod0.5"
dat_name="mod0.5.dat"

jjm.in  <- read.dat(iFilename = paste(dat_name,sep=""),iPath=inputPath)
jjm.out <- readList(file.path(inputPath,paste("arc/",run_name,"_r.rep",sep="")))
jjm.ypr <- readYPR(file.path(inputPath,paste( "arc/",run_name,".yld",  sep="")))

#-------------------------------------------------------------------------------
# Create diagnostics
#-------------------------------------------------------------------------------
pdf(paste(resultPath,"summary_",run_name,".pdf",sep=""),height=29.7/2.54,width=21/2.54,pointsize = 16, bg = "white")
diagnostics(jjm.out,jjm.in,jjm.ypr,what=c("input","fit","projections","ypr"))
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

lstOuts   <- list(
  Model_0.0= jjm0.0,
  Model_0.1= jjm0.1,
  Model_0.2= jjm0.2,
  Model_0.3= jjm0.3,
  Model_0.4= jjm0.4,
  Model_0.5= jjm0.5
                  )
lstOuts   <- list(
  Model_0.5= jjm0.5,
  Model_0.6= jjm0.6
)

pdf(paste(outputPath,"Compare_1_4.pdf",sep=""),height=29.7/2.54,width=21/2.54,pointsize = 24, bg = "white")

jjm.n2    <- readList(file.path(outputPath,"n3_r.rep"))
lstOuts   <- list(Model_N1=jjm.n1,Model_N2=jjm.n2)

jjm.s1    <- readList(file.path(outputPath,"s1_r.rep"))
jjm.s2    <- readList(file.path(outputPath,"s2_r.rep"))
lstOuts   <- list(Model_S1=jjm.s1,Model_S2=jjm.s2)

compareTime(lstOuts,"SSB",SD=T,Sum=NULL,legendPos="top")
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,startYear=1953)
compareTime(lstOuts,"R",SD=T)
compareTime(lstOuts,"TotBiom",SD=T)

compareMatrix(lstOuts,"TotF",SD=F,Sum=NULL,YrInd=jjm.s2$Yr,Apply=mean)
compareMatrix(lstOuts,"N",   SD=F,         YrInd=jjm.s2$Yr,Apply=sum)
dev.off()

pdf(paste(resultPath,"Compare_0",".pdf",sep=""),height=29.7/2.54,width=21/2.54,pointsize = 16, bg = "white")
#-------------------------------------------------------------------------------
# Numerical compare runs
#-------------------------------------------------------------------------------

# Likelihood table
tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
write.csv(tab,file=file.path(inputPath,"LikelihoodTable2013.csv"),row.names=F)
system("Likelihoodtable2013.csv")
tab
inputPath
# Risk tables
jjm.mod7  <- readList(file.path(inputPath,"mod7_r.rep"))
jjm.mod6  <- readList(file.path(inputPath,"mod6_r.rep"))
lstOuts   <- list(Model_6=jjm.mod6,Model_7=jjm.mod7)

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
