library(PBSadmb)
# library(lattice)
library(plotrix)
library(ggplot2)
library(RColorBrewer)
library(doBy)
library(tidyverse)
library(data.table)
#install.packages("gridExtra")
#install.packages("PBSadmb")
library(gridExtra)
source("plot_ind.R")
source("~/Dropbox/R_Common/ADFunctions.r")
#source(paste(fndir,"adfunctions2.r",sep=""))
#source(paste0(here::here("R"),"/diagnostics_v2.r") )
#source("R/ADMB2R_15102013.r")
#source("R/compareRuns.r")
mytheme <- theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank() )
mytheme <- mytheme + theme(text=element_text(size=18)) + theme(axis.title.x=element_text(size=24) ,axis.title.y=element_text(size=24))
mytheme <- mytheme + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_line(colour="grey60", linetype="dashed"), panel.grid.major.y = element_blank() )
mytheme <- mytheme + theme( panel.background = element_rect(fill="white"), panel.border = element_rect(colour="black", fill=NA, size=1))
get_params <- function(lst){
	return(	data.table(rbind(
		cbind("q",do.call(cbind,lapply(lstOuts,function(x){round(x[["q_1"]][1,2],2)}))),
		cbind("Npars",do.call(cbind,lapply(lstOuts,function(x){round(x[["Num_parameters_Est"]],2)}))),
		cbind("M",do.call(cbind,lapply(lstOuts,function(x){round(x[["Mest"]][2],2)}))),
		cbind("SigmaR",do.call(cbind,lapply(lstOuts,function(x){round(x[["Sigmar"]][2],2)}))),
		cbind("EffN_Fish",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["EffN_Fsh_1"]][,2]),2)}))),
		cbind("sdnr_Fish_Age",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["sdnr_age_fsh_1"]][,2]),2)}))),
		cbind("EffN_Surv",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["EffN_Survey_1"]][,2]),2)}))),
		cbind("sdnr_Surv_Age",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["sdnr_age_ind_1"]][,2]),2)}))),
		cbind("sdnr_Survey",do.call(cbind,lapply(lstOuts,function(x){round((x[["sdnr_ind_1"]]),2)}))), tab,
		cbind("F2016",do.call(cbind,lapply(lstOuts,function(x){round((x[["F_fsh_1"]][40,2]),2)}))),
		cbind("F2016/F40%",do.call(cbind,lapply(lstOuts,function(x){round((x[["F_fsh_1"]][40,2])/(x[["F40_est"]]),2)}))),
		cbind("B 1977",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][1,2]),0)}))),
		cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][1,3])/(x[["TotBiom"]][1,2])*100,0)}))),
		cbind("B 2016",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][40,2]),0)}))),
		cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][40,3])/(x[["TotBiom"]][40,2])*100,0)}))),
		cbind("2001 YC",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][26,2]),0)}))),
		cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][26,3])/(x[["R"]][26,2])*100,0)}))),
		cbind("2006 YC",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][31,2]),0)}))),
		cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][31,3])/(x[["R"]][31,2])*100,0)})))
		))
	)
}

na.pad <- function(x,len){
    x[1:len]
}
makePaddedDataFrame <- function(l,...){
    maxlen <- max(sapply(l,length))
    data.frame(lapply(l,na.pad,len=maxlen),...)
}
