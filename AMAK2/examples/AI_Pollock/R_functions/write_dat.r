write_dat<-function(data=test,data_file="test.dat"){
 options(list(scipen=6))
 T1<-noquote("# General Parameters")
  write(data$T1,paste(data_file))
  ##Start YEAR
  write(data$SY,paste(data_file),append = T)
  ##END YEAR
  write(data$EY,paste(data_file),append = T)
  ## First Age class
  write(data$FAC,paste(data_file),append = T)
  ##Last Age Class
  write(data$LAC,paste(data_file),append = T)
  ## numbers of length intervals																																									
   write("41",paste(data_file),append = T)	
  ##Length bin
   write(c(10:50),paste(data_file),append = T)
  ## 10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50
 #	=======================================		
 
  ## Number of Fisheries
  write(data$NF,paste(data_file),append = T)

 T2<-noquote("#Fshry name")
   write(T2,paste(data_file),append = T)
   write(data$F_NAME,paste(data_file),append = T)

 T3<-noquote("# Catch Biomass:")
  write(T3,paste(data_file),append = T)
  write.table(t(data$CATCH$TONS/1000),paste(data_file),col.names=F,row.names=F,append = T)

 T4<-noquote("#Yr specific cv on catch:")
  write(T4,paste(data_file),append = T)
  write.table(t(data$Catch_CV),paste(data_file),col.names=F,row.names=F,append = T)

 T5<- noquote("# Number Years With Fishery Age Data:")
  write(T5,paste(data_file),append = T)
  write.table(length(data$F_NAGE$YEAR),paste(data_file),col.names=F,row.names=F,append = T)
 
 T25<- noquote("# Number Years With Fishery Length Data:")
  write(T25,paste(data_file),append = T)
  write.table(0,paste(data_file),col.names=F,row.names=F,append = T)
 
 T6<-noquote("# Fishery 1 years")
  write(T6,paste(data_file),append = T)
  write.table(t(data$F_NAGE$YEAR),paste(data_file),col.names=F,row.names=F,append = T)

 T7<-noquote("# Fishery 1 multinomial sample size")
  write(T7,paste(data_file),append = T)
  write.table(t(data$F_MN_SIZE),paste(data_file),col.names=F,row.names=F,append = T)

 T8<- noquote("# Fishery 1 numbers at age proportion")
  write(T8,paste(data_file),append = T)
  write.table(round(data$F_NAGE[,2:15],3),paste(data_file),,col.names=F,row.names=F,append = T)

 T9<- noquote("# Fishery 1 Weight at age 1978-2011")
  write(T9,paste(data_file),append = T)
  write.table(round(data$F_WEIGHT$AGE_WEIGHT[,2:15]),paste(data_file),col.names=F,row.names=F,append = T)

 T10<- noquote("#  Number of Surveys")
  write(T10,paste(data_file),append = T)
  write(data$NS,paste(data_file),append = T)

 T11<- noquote("#	Name of Surveys")
  write(T11,paste(data_file),append = T)
   write(data$S_name,paste(data_file),append = T)
   write(length(data$S_BIO$YEAR),paste(data_file),append = T)

 T12<- noquote("# Survey 1 Years of survey biomass")
  write(T12,paste(data_file),append = T)
  write.table(t(matrix(data$S_BIO$YEAR)),paste(data_file),col.names=F,row.names=F,append = T)
  
 T13<- noquote("# Survey 1 Number of years")
  write(T13,paste(data_file),append = T)
  write(length(data$S_BIO$YEAR),paste(data_file),append = T)

 T14<- noquote("# Survey 1 Biomass values")
  write(T14,paste(data_file),append = T)
  write.table(t(matrix(round(data$S_BIO$BIO/1000,3))),paste(data_file),col.names=F,row.names=F,append = T)

 T15<- noquote("# Survey 1 Biomass stdev")
  write(T15,paste(data_file),append = T)
  write.table(t(matrix(round(data$S_BIO$B_STDEV/1000,3))),paste(data_file),col.names=F,row.names=F,append = T)

 T16<- noquote("# Number of Surveys with ages")
  write(T16,paste(data_file),append = T)
   write(length(data$S_NAGE$YEAR),paste(data_file),append = T)
   
  T161<- noquote("# Number of Surveys with lengths")
  write(T161,paste(data_file),append = T)
   write(0,paste(data_file),append = T) 

 T17<- noquote("# Survey 1 Years")
  write(T17,paste(data_file),append = T)
  write.table(t(matrix(data$S_NAGE$YEAR)),paste(data_file),col.names=F,row.names=F,append = T)

 T18<- noquote("# Survey 1 Multinomial sample size")
  write(T18,paste(data_file),append = T)
  write.table(t(matrix(data$S_MN_SIZE)),paste(data_file),col.names=F,row.names=F,append = T)

 T19<-noquote("# Survey 1 numbers at age proportions")
  write(T19,paste(data_file),append = T)
  write.table(round(data$S_NAGE[,2:15],3),paste(data_file),col.names=F,row.names=F,append = T)

 T20<- noquote("# Survey 1 weight at age")
  write(T20,paste(data_file),append = T)
  write.table(round(data$S_WEIGHT$AGE_WEIGHT[,2:15]),paste(data_file),col.names=F,row.names=F,append = T)

 T21<-noquote("# Population weight at Age")
  write(T21,paste(data_file),append = T)
  write.table(t(matrix(round(data$PP_WEIGHT_AGE,3))),paste(data_file),col.names=F,row.names=F,append = T)

 T22<-noquote("# Maturity at Age ")
  write(T22,paste(data_file),append = T)
  write.table(data$GMat*0.5,paste(data_file),col.names=F,row.names=F,append = T)

 T23<-noquote("# Month of Spawning")
  write(T23,paste(data_file),append = T)
  write(data$MSPWN,paste(data_file),append = T)

  T24<-noquote("# Ageing error")
  write(T24,paste(data_file),append = T)
  write.table(data$AERR,paste(data_file),col.names=F,row.names=F,append = T)
  write(" ",paste(data_file),append = T)
  write(" ",paste(data_file),append = T)
 }




