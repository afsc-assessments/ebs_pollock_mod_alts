library(RODBC)
library(mgcv)
setwd("Z:/Users/steve.barbeaux/My Documents/Private_Data/Stock_Assessments/Aleutian_Assessment_11/data")


AFSC=odbcConnect("AFSC","sbarb","sbarb$1011")

FAge=sqlQuery(AFSC,"SELECT NORPAC.FOREIGN_HAUL.YEAR,
  NORPAC.FOREIGN_HAUL.GENERIC_AREA,
  NORPAC.FOREIGN_AGE.SPECIES,
  NORPAC.FOREIGN_AGE.SEX,
  NORPAC.FOREIGN_AGE.LENGTH,
  NORPAC.FOREIGN_AGE.AGE,
  NORPAC.FOREIGN_AGE.INDIV_WEIGHT,
  NORPAC.FOREIGN_HAUL.LATITUDE,
  NORPAC.FOREIGN_HAUL.E_W,
  NORPAC.FOREIGN_HAUL.LONGITUDE
FROM NORPAC.FOREIGN_HAUL
INNER JOIN NORPAC.FOREIGN_AGE
ON NORPAC.FOREIGN_HAUL.HAUL_JOIN = NORPAC.FOREIGN_AGE.HAUL_JOIN
WHERE NORPAC.FOREIGN_HAUL.YEAR   > 1977
AND NORPAC.FOREIGN_HAUL.GENERIC_AREA BETWEEN 539 AND 544
AND NORPAC.FOREIGN_AGE.SPECIES = 201
ORDER BY NORPAC.FOREIGN_HAUL.YEAR")


FAge1<-subset(FAge,FAge$GENERIC_AREA==541&FAge$LATITUDE<5300)
FAge2<-subset(FAge,FAge$GENERIC_AREA==542&FAge$LATITUDE<5230)
FAge3<-subset(FAge,FAge$GENERIC_AREA==543&FAge$LATITUDE<5330)
FAge<-rbind(FAge1,FAge2,FAge3)
FAge<-FAge[order(FAge$YEAR,FAge$GENERIC_AREA),]
remove(FAge1,FAge2,FAge3)
FAge<-subset(FAge,select= -c(LATITUDE,LONGITUDE,E_W) )

FLength=sqlQuery(AFSC," SELECT NORPAC.FOREIGN_HAUL.YEAR,
  NORPAC.FOREIGN_HAUL.GENERIC_AREA,
  NORPAC.FOREIGN_LENGTH.SPECIES,
  NORPAC.FOREIGN_LENGTH.SEX,
  NORPAC.FOREIGN_LENGTH.SIZE_GROUP,
  NORPAC.FOREIGN_LENGTH.FREQUENCY,
  NORPAC.FOREIGN_HAUL.LATITUDE,
  NORPAC.FOREIGN_HAUL.E_W,
  NORPAC.FOREIGN_HAUL.LONGITUDE
FROM NORPAC.FOREIGN_HAUL
INNER JOIN NORPAC.FOREIGN_LENGTH
ON NORPAC.FOREIGN_HAUL.HAUL_JOIN = NORPAC.FOREIGN_LENGTH.HAUL_JOIN
WHERE NORPAC.FOREIGN_HAUL.YEAR   > 1977
AND NORPAC.FOREIGN_HAUL.GENERIC_AREA BETWEEN 539 AND 544
AND NORPAC.FOREIGN_LENGTH.SPECIES = 201
ORDER BY NORPAC.FOREIGN_HAUL.YEAR")

FLength1<-subset(FLength,FLength$GENERIC_AREA==541&FLength$LATITUDE<5300)
FLength2<-subset(FLength,FLength$GENERIC_AREA==542&FLength$LATITUDE<5230)
FLength3<-subset(FLength,FLength$GENERIC_AREA==543&FLength$LATITUDE<5330)
FLength<-rbind(FLength1,FLength2,FLength3)
FLength<-FLength[order(FLength$YEAR,FLength$GENERIC_AREA),]
remove(FLength1,FLength2,FLength3)
FLength<-subset(FLength,select= -c(LATITUDE,LONGITUDE,E_W) )


DAge=sqlQuery(AFSC,"SELECT OBSINT.DEBRIEFED_AGE.YEAR,
  OBSINT.DEBRIEFED_AGE.NMFS_AREA,
  OBSINT.DEBRIEFED_AGE.SPECIES,
  OBSINT.DEBRIEFED_AGE.SEX,
  OBSINT.DEBRIEFED_AGE.LENGTH,
  OBSINT.DEBRIEFED_AGE.AGE,
  OBSINT.DEBRIEFED_AGE.WEIGHT,
  OBSINT.DEBRIEFED_AGE.LATDD_END,
  OBSINT.DEBRIEFED_AGE.LONDD_END
FROM OBSINT.DEBRIEFED_AGE
WHERE OBSINT.DEBRIEFED_AGE.NMFS_AREA BETWEEN 539 AND 544
AND OBSINT.DEBRIEFED_AGE.SPECIES = 201
ORDER BY OBSINT.DEBRIEFED_AGE.YEAR" )

DAge1<-subset(DAge,DAge$NMFS_AREA==541&DAge$LATDD_END<53.00)
DAge2<-subset(DAge,DAge$NMFS_AREA==542&DAge$LATDD_END<52.50)
DAge3<-subset(DAge,DAge$NMFS_AREA==543&DAge$LATDD_END<53.50)
DAge<-rbind(DAge1,DAge2,DAge3)
DAge<-DAge[order(DAge$YEAR,DAge$NMFS_AREA),]
remove(DAge1,DAge2,DAge3)
DAge<-subset(DAge,select= -c(LATDD_END,LONDD_END) )


DLength=sqlQuery(AFSC,"SELECT OBSINT.DEBRIEFED_LENGTH.YEAR,
  OBSINT.DEBRIEFED_LENGTH.NMFS_AREA,
  OBSINT.DEBRIEFED_LENGTH.SPECIES,
  OBSINT.DEBRIEFED_LENGTH.SEX,
  OBSINT.DEBRIEFED_LENGTH.LENGTH,
  OBSINT.DEBRIEFED_LENGTH.FREQUENCY,
  OBSINT.DEBRIEFED_LENGTH.LATDD_END,
  OBSINT.DEBRIEFED_LENGTH.LONDD_END
FROM OBSINT.DEBRIEFED_LENGTH
WHERE OBSINT.DEBRIEFED_LENGTH.NMFS_AREA BETWEEN 539 AND 544
AND OBSINT.DEBRIEFED_LENGTH.SPECIES = 201
ORDER BY OBSINT.DEBRIEFED_LENGTH.YEAR")

DLength1<-subset(DLength,DLength$NMFS_AREA==541&DLength$LATDD_END<53.00)
DLength2<-subset(DLength,DLength$NMFS_AREA==542&DLength$LATDD_END<52.50)
DLength3<-subset(DLength,DLength$NMFS_AREA==543&DLength$LATDD_END<53.50)
DLength<-rbind(DLength1,DLength2,DLength3)
DLength<-DLength[order(DLength$YEAR,DLength$NMFS_AREA),]
remove(DLength1,DLength2,DLength3)
DLength<-subset(DLength,select= -c(LATDD_END,LONDD_END) )



AAge<-read.csv("Acoustic_survey_data.csv",header=T)
AAge$WEIGHT<-AAge$WEIGHT/1000
ALength<-aggregate(list(FREQ=AAge$LENGTH),by=list(YEAR=AAge$YEAR,NMFS_AREA=AAge$NMFS_AREA,SPECIES=AAge$SPECIES,SEX=AAge$SEX,LENGTH=AAge$LENGTH),FUN=length)
ALength<-ALength[order(ALength$YEAR,ALength$LENGTH),]

names(FLength)<-names(DLength)
names(FAge)<-names(DAge)
names(AAge)<-names(DAge)
names(ALength)<-names(DLength)

FAge<-rbind(FAge,DAge,AAge)
FLength<-rbind(FLength,DLength,ALength)
remove(DAge,DLength,ALength,AAge)

## Fishery AGE
FAge_w<-subset(FAge,!is.na(FAge$AGE))
FAge_w<-subset(FAge_w,!is.na(FAge_w$LENGTH))
FAge_w1<-subset(FAge_w,FAge_w$YEAR<1999&FAge_w$YEAR!=1988&FAge_w$YEAR!=1989&FAge_w$YEAR!=1990&FAge_w$YEAR!=1991&FAge_w$YEAR!=1992&FAge_w$YEAR!=1993&FAge_w$YEAR!=1997)
FAge_w2<-subset(FAge_w,FAge_w$YEAR>2005&FAge_w$YEAR<2009)
FAge_w<-rbind(FAge_w1,FAge_w2)



FLength_w1<-subset(FLength,FLength$YEAR<1999&FLength$YEAR!=1988&FLength$YEAR!=1989&FLength$YEAR!=1990&FLength$YEAR!=1991&FLength$YEAR!=1992&FLength$YEAR!=1993&FLength$YEAR!=1997)
FLength_w2<-subset(FLength,FLength$YEAR>2005&FLength$YEAR<2009)
FLength_w<-rbind(FLength_w1,FLength_w2)


FAge_w$SEX1<-FAge_w$SEX
FAge_w$SEX<-as.character(FAge_w$SEX)
FAge_w$SEX[FAge_w$SEX1=="M"]=1
FAge_w$SEX[FAge_w$SEX1=='F']=2
FAge_w$SEX[FAge_w$SEX1=='U']=3
FAge_w$LENGTH<-FAge_w$LENGTH*10
FAge_w$WEIGHT<-FAge_w$WEIGHT*1000

FLength_w$SEX1<-FLength_w$SEX
FLength_w$SEX<-as.character(FLength_w$SEX)
FLength_w$SEX[FLength_w$SEX1=='M']=1
FLength_w$SEX[FLength_w$SEX1=='F']=2
FLength_w$SEX[FLength_w$SEX1=='U']=3

FLength_w$LENGTH<-FLength_w$LENGTH*10


source("find_AL.r")
years<-sort(unique(FAge_w$YEAR))
z<-data.frame(matrix(ncol=4,nrow=1))
names(z)<-c("YEAR","SEX","AGE","LENGTH" )
for( i in 1:length(years)){
 x=find_AL(FAge_w,FLength_w,years[i])
 y1=data.frame(YEAR=years[i],SEX=1,AGE=x$len1$age,LENGTH=x$len1$tl)
 y2=data.frame(YEAR=years[i],SEX=2,AGE=x$len1$age,LENGTH=x$len1$tl)
 z=rbind(z,y1,y2)
 }
 F_length_age<-subset(z,is.na(z$YEAR)==F)

F_Total_AGE<-aggregate(list(FREQ=F_length_age$AGE),by=list(SEX=F_length_age$SEX,AGE=F_length_age$AGE,YEAR=F_length_age$YEAR),FUN=length)

 GRID<-expand.grid(SEX=c(1,2),YEAR=years,AGE=c(1:15))
 F_Total_AGE<-merge(GRID,F_Total_AGE, all=T)
 F_Total_AGE$FREQ[is.na(F_Total_AGE$FREQ)]=0
 Total<-aggregate(list(Total=F_length_age$AGE),by=list(YEAR=F_length_age$YEAR),FUN=length)
 F_Total_AGE<-merge(Total,F_Total_AGE, all=T)
 F_Total_AGE$PROB<-F_Total_AGE$FREQ/F_Total_AGE$Total
 F_Total_AGE<-subset(F_Total_AGE,select=-c(Total))
 F_Total_AGE<-F_Total_AGE[order(F_Total_AGE$YEAR,F_Total_AGE$SEX,F_Total_AGE$AGE),]



 F_Total_LENGTH<-aggregate(list(FREQ=F_length_age$AGE),by=list(SEX=F_length_age$SEX,LENGTH=F_length_age$LENGTH,YEAR=F_length_age$YEAR),FUN=length)

 GRID<-expand.grid(YEAR=years,LENGTH=seq(min(FLength_w$LENGTH),max(FLength_w$LENGTH),10))
 F_Total_LENGTH<-merge(GRID,F_Total_LENGTH, all=T)
 F_Total_LENGTH$FREQ[is.na(F_Total_LENGTH$FREQ)]=0
 Total<-aggregate(list(Total=F_length_age$AGE),by=list(YEAR=F_length_age$YEAR),FUN=length)
 F_Total_LENGTH<-merge(Total,F_Total_LENGTH, all=T)
 F_Total_LENGTH$PROB<-F_Total_LENGTH$FREQ/F_Total_LENGTH$Total
 F_Total_LENGTH<-subset(F_Total_LENGTH,select=-c(Total))
 F_Total_LENGTH<-F_Total_LENGTH[order(F_Total_LENGTH$YEAR,F_Total_LENGTH$SEX,F_Total_LENGTH$LENGTH),]
 remove(Total,GRID)


 ## Fishery WEIGHT at Age data
source("SWEIGHT.r")
F_WEIGHT<-SWEIGHT(data=FAge_w,Total_A=F_Total_AGE,fyear=2011)
F_Total_AGE<-merge(F_Total_AGE,F_WEIGHT[[1]])
F_Total_AGE<-F_Total_AGE[order(F_Total_AGE$YEAR,F_Total_AGE$SEX,F_Total_AGE$AGE),]





## Fishery Catch Data
source("get_catch.r")
CATCH<-GET_CATCH()


F_Total_WB<-merge(F_Total_AGE,CATCH[[1]],all.x=T)
F_Total_WB<-F_Total_WB[order(F_Total_WB$YEAR,F_Total_WB$AGE),]


F_Total_WB$N_AGE<-F_Total_WB$PROB*F_Total_WB$TONS/F_Total_WB$WEIGHT2*1000000
F_Total_WB$W_AGE<-F_Total_WB$PROB*F_Total_WB$TONS

F_total<-aggregate(list(FREQ=F_Total_AGE$FREQ),by=list(YEAR=F_Total_AGE$YEAR,AGE=F_Total_AGE$AGE),FUN=sum)
F_number<-aggregate(list(NUMBER=F_Total_AGE$FREQ),by=list(YEAR=F_Total_AGE$YEAR),FUN=sum)
F_total<-merge(F_total,F_number)
F_total$PROB<-F_total$FREQ/F_total$NUMBER
F_total<-F_total[order(F_total$YEAR,F_total$AGE),]



F_NAGE<-data.frame(matrix(ncol=15,nrow=length(years)))
  names(F_NAGE)<-c("YEAR","A2","A3" ,"A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15")
  F_NAGE$YEAR<-years

 for( i in 1:length(F_NAGE$YEAR)){
  F_NAGE[i,2:15]<-F_total$PROB[F_total$YEAR==F_NAGE$YEAR[i]][1:14]
  }




Catch_CV<-rep(0.05,length(unique(CATCH[[1]]$YEARS)))

##Total_FWB$C_LENGTH<-Total_WB$Prob*Total_FWB$TONS

CATCH[[1]]$TONS         ## CATCH BIOMASS
F_NAGE                  ## fishery numbers at age
F_WEIGHT$S_AGE_WEIGHT   ## fishery WEIGTH at age

a<-noquote("# General Parameters")
write.table(a,"test.dat",row.names=F,col.names=F)
