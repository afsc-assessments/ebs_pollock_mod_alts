#-------------------------------------------------------------------------------
# Script to run the assessment of Jack Mackerel and look at outputs
#
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
reposDir    <- "~/OneDrive/Models/Atka/2016/"
codePath    <- file.path(reposDir,"R/")
inputPath   <- file.path(reposDir,"")
outputPath  <- file.path(reposDir,"")
resultPath  <- file.path(reposDir,"Results/Assessment/")
setwd(codePath)

# Specify control file
controlFile <- "mod1/amak.dat"

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

lstOuts   <- list( Model_1= mod1, Model_2= mod2, Model_3= mod3, Model_4= mod4, Model_5= mod5, Model_6= mod6, Model_7= mod7 )
lstOuts   <- list( Model_1= mod1, Model_2= mod2, Model_3= mod3, Model_4= mod4)
tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
tab
lstOuts   <- list( Model_1= mod1, Model_5= mod5, Model_6= mod6, Model_7= mod7 )
lstOuts   <- list( Model_1= mod1, Model_5= mod5, Model_7= mod7 )
lstOuts   <- list( Model_1= mod1, Model_7= mod7 )
lstOuts   <- list( Model_1= mod1, Model_7= mod7 )
names(mod1)
mod1$M
mod2$M

pdf(paste(outputPath,"Compare_1_4.pdf",sep=""),height=6,width=10,pointsize = 24, bg = "white")

compareTime(lstOuts,"SSB",SD=T,Sum=NULL,legendPos="right",startYear=1980)
compareTime(lstOuts,"SSB",SD=F,Sum=NULL,startYear=1953)
compareTime(lstOuts,"R",SD=T)
compareTime(lstOuts,"TotBiom",SD=T)
compareTime(lstOuts,"TotBiom",SD=F)
library(data.table)
df <-data.frame(Age=1:11,Model_1=mod1$M[1,],
  Model_2=mod2$M[1,],
  Model_3=mod3$M[1,],
  Model_4=mod4$M[1,],Maturity=mod1$mature)
df <- data.table(df)
df
dt<-melt(df,id="Age")
mod1$mature
names(mod1)
ggplot(dt,aes(x=Age,y=value,colour=variable)) + geom_line(size=2) + labs(x="Age",y="M or proportion mature)") + mytheme
compareMatrix(lstOuts,"TotF",SD=TRUE,Sum=NULL,YrInd=mod1$Yr,Apply=mean,)
compareMatrix(lstOuts,"N",   SD=F,         YrInd=mod1$Yr,Apply=sum)
dev.off()

bdf(paste(resultPath,"Compare_0",".pdf",sep=""),height=29.7/2.54,width=21/2.54,pointsize = 16, bg = "white")
#-------------------------------------------------------------------------------
# Numerical compare runs
#-------------------------------------------------------------------------------

# Likelihood table
tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
write.csv(tab,file=file.path(inputPath,"LikelihoodTableM.csv"),row.names=F)
system("open LikelihoodTableM.csv")

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
#' Plot natural mortality
#'
#' @param replist List object created by read_admb function
#' @return Plot natural mortality over time and size
#' @export
plot_naturalmortality <- function(replist){
  A    <- replist
  df   <- data.frame(A$M)
  colnames(df) <- A$mid_points
  nrow   <- dim(A$M)[1]
  # Always saves for both sexes???
  df$sex <- c(rep(1,length=nrow/2),rep(2,length=nrow/2))
  df$Year <- A$mod_yrs
  mdf    <- melt(df,id=c("sex","Year"))
  
  p <- ggplot(mdf,aes(x=Year,y=as.double(variable),z=value))
  p <- p + geom_tile(aes(fill = value)) 
  p <- p + stat_contour(geom="polygon", aes(fill=(value)))
  p <- p + labs(x="Year",y="size bin",fill="M")
  p <- p + facet_wrap(~sex,scale="free")
  p2 <- ggplot(mdf,aes(x=Year,y=value))
  p2 <- p2 + geom_line() + ggtheme + labs(y="Natural mortality")
  plot_multiple(p2,p)
}

plot_survey <- function(replist){
  A <- replist
  df <- as.data.frame(A$Obs_Survey_1)
  colnames(df) <- c("year","obs","pred","stderr","res1","res2")
  sd <- df$stderr
  df$lb <- exp(log(df$obs)-1.96*log(sd))
  df$ub <- exp(log(df$obs)+1.96*log(sd))

  p  <- ggplot(df,aes(year,obs))
# p  <- p + geom_point(aes(col=sex))
  p  <- p + geom_pointrange(aes(year,obs,ymax=obs,ymin=lb))
  p  <- p + labs(x="Year",y="Survey")
# Fitted CPUE
  pCPUEfit <- p + geom_line(data=df,aes(year,pred))
  return(pCPUEfit)
}
