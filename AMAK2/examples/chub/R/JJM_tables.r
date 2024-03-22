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

  # Set paths
codePath    <- "N:/Projecten/SouthPacific/2012/SWG/R/"
inputPath   <- "D:/Repository/Dropbox/jurel/Codes/Work/"
outputPath  <- "N:/Projecten/SouthPacific/2012/SWG/Results/"
setwd(codePath)

  # Read in output data to create tables
jjm.comb  <- readList(file.path(inputPath,"arc/comb_r.rep"))
jjm.north <- readList(file.path(inputPath,"arc/north_r.rep"))
jjm.south <- readList(file.path(inputPath,"arc/south_r.rep"))

  # Likelihood table
tab       <- as.data.frame(cbind(jjm.comb$Like_Comp_names,round(jjm.comb$Like_Comp,1)))
tab       <- cbind(tab,round(jjm.north$Like_Comp,1),round(jjm.south$Like_Comp,1))
colnames(tab) <- c("name","comb","north","south")
tab$name  <- c("Catch biomass","Fishery age compositions","Fishery length composition",
               "Fishery selectivity"," ","Indices age composition","Indices selectivity",
               "Stock-recruitment","F penalty","Priors indices Q","priors","Residual","Total")

write.csv(tab,file=file.path(inputPath,"LikelihoodTable2012.csv"),row.names=F)
