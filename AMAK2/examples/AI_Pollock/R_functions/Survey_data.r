library(RODBC)
library(mgcv)
setwd("Z:/Users/steve.barbeaux/My Documents/Private_Data/Stock_Assessments/Aleutian_Assessment_11/data")


AFSC=odbcConnect("AFSC","sbarb","sbarb$1011")

## Survey Age Data

Age=sqlQuery(AFSC, "SELECT RACEBASE.SPECIMEN.REGION,
  TO_CHAR(RACEBASE.HAUL.START_TIME, 'yyyy') AS YEAR,
  RACEBASE.SPECIMEN.CRUISE,
  RACEBASE.HAUL.VESSEL,
  RACEBASE.SPECIMEN.HAUL,
  RACEBASE.SPECIMEN.SPECIES_CODE,
  RACEBASE.SPECIMEN.LENGTH,
  RACEBASE.SPECIMEN.SEX,
  RACEBASE.SPECIMEN.WEIGHT,
  RACEBASE.SPECIMEN.MATURITY,
  RACEBASE.SPECIMEN.AGE,
  RACEBASE.HAUL.END_LONGITUDE,
  RACEBASE.HAUL.HAUL_TYPE,
  RACEBASE.HAUL.PERFORMANCE,
  RACEBASE.SPECIMEN.SPECIMEN_SAMPLE_TYPE,
  RACEBASE.SPECIMEN.SPECIMENID,
  RACEBASE.SPECIMEN.BIOSTRATUM
FROM RACEBASE.SPECIMEN
INNER JOIN RACEBASE.HAUL
ON RACEBASE.SPECIMEN.HAULJOIN                 = RACEBASE.HAUL.HAULJOIN
AND RACEBASE.HAUL.CRUISEJOIN                  = RACEBASE.SPECIMEN.CRUISEJOIN
WHERE RACEBASE.SPECIMEN.REGION                = 'AI'
AND RACEBASE.SPECIMEN.SPECIES_CODE            = 21740
AND RACEBASE.HAUL.HAUL_TYPE                   = 3
AND RACEBASE.HAUL.PERFORMANCE                >= 0
AND TO_CHAR(RACEBASE.HAUL.START_TIME, 'yyyy') > 1977
ORDER BY TO_CHAR(RACEBASE.HAUL.START_TIME, 'yyyy')")

Age_w<-subset(Age,Age$END_LONGITUDE > 0 | Age$END_LONGITUDE < -170)

Age_w<-subset(Age_w,is.na(Age_w$AGE)==F)
Age_w$AGE1<-Age_w$AGE
Age_w$AGE1[Age_w$AGE1>=15]=15
#Age_w$END_LONG<-Age_w$END_LONGITUDE
#Age_w$END_LONG[Age_w$END_LONGITUDE<0]<- -Age_w$END_LONGITUDE[Age_w$END_LONGITUDE<0]
#Age_w$END_LONG[Age_w$END_LONGITUDE>0]<- Age_w$END_LONGITUDE[Age_w$END_LONGITUDE>0]+ (180-Age_w$END_LONGITUDE[Age_w$END_LONGITUDE>0])
#Age_w<-subset(Age_w,Age_w$SEX!=3)


##Length Data

Length=sqlQuery(AFSC,"SELECT RACEBASE.LENGTH.REGION,
  TO_CHAR(RACEBASE.HAUL.START_TIME, 'yyyy') AS YEAR,
  RACEBASE.LENGTH.CRUISE,
  RACEBASE.LENGTH.VESSEL,
  RACEBASE.LENGTH.HAUL,
  RACEBASE.LENGTH.SPECIES_CODE,
  RACEBASE.LENGTH.SEX,
  RACEBASE.LENGTH.LENGTH,
  RACEBASE.LENGTH.FREQUENCY,
  RACEBASE.HAUL.END_LONGITUDE,
  RACEBASE.HAUL.GEAR,
  RACEBASE.HAUL.STRATUM,
  RACEBASE.HAUL.HAUL_TYPE
FROM RACEBASE.LENGTH
INNER JOIN RACEBASE.HAUL
ON RACEBASE.LENGTH.CRUISEJOIN                 = RACEBASE.HAUL.CRUISEJOIN
AND RACEBASE.LENGTH.HAULJOIN                  = RACEBASE.HAUL.HAULJOIN
WHERE RACEBASE.LENGTH.REGION                  = 'AI'
AND RACEBASE.LENGTH.SPECIES_CODE              = 21740
AND TO_CHAR(RACEBASE.HAUL.START_TIME, 'yyyy') >= 1980
AND RACEBASE.HAUL.HAUL_TYPE                   = 3 ")


Length_w<-subset(Length,Length$END_LONGITUDE > 0 | Length$END_LONGITUDE < -170)
Length_w<-subset(Length_w,Length_w$YEAR!=1982)
Length_w<-subset(Length_w,Length_w$YEAR!=2010)
Length_w<-subset(Length_w,Length_w$SEX!=3)



source("find_AL.r")
years<-sort(unique(Age_w$YEAR))
z<-data.frame(matrix(ncol=4,nrow=1))
names(z)<-c("YEAR","SEX","AGE","LENGTH" )
for( i in 1:length(years)){
 x=find_AL(Age_w,Length_w,years[i])
 y1=data.frame(YEAR=years[i],SEX=1,AGE=x$len1$age,LENGTH=x$len1$tl)
 y2=data.frame(YEAR=years[i],SEX=2,AGE=x$len2$age,LENGTH=x$len2$tl)
 z=rbind(z,y1,y2)
 }
 
 S_length_age<-subset(z,is.na(z$YEAR)==F)
 
 Total_AGE<-aggregate(list(FREQ=S_length_age$AGE),by=list(SEX=S_length_age$SEX,AGE=S_length_age$AGE,YEAR=S_length_age$YEAR),FUN=length)

 GRID<-expand.grid(SEX=c(1,2),YEAR=years,AGE=c(1:15))
 Total_AGE<-merge(GRID,Total_AGE, all=T)
 Total_AGE$FREQ[is.na(Total_AGE$FREQ)]=0
 Total<-aggregate(list(Total=S_length_age$AGE),by=list(YEAR=S_length_age$YEAR),FUN=length)
 Total_AGE<-merge(Total,Total_AGE, all=T)
 Total_AGE$PROB<-Total_AGE$FREQ/Total_AGE$Total
 Total_AGE<-subset(Total_AGE,select=-c(Total))
 Total_AGE<-Total_AGE[order(Total_AGE$YEAR,Total_AGE$SEX,Total_AGE$AGE),]



 Total_LENGTH<-aggregate(list(FREQ=S_length_age$AGE),by=list(SEX=S_length_age$SEX,LENGTH=S_length_age$LENGTH,YEAR=S_length_age$YEAR),FUN=length)

 GRID<-expand.grid(YEAR=years,LENGTH=seq(min(Length_w$LENGTH),max(Length_w$LENGTH),10))
 Total_LENGTH<-merge(GRID,Total_LENGTH, all=T)
 Total_LENGTH$FREQ[is.na(Total_LENGTH$FREQ)]=0
 Total<-aggregate(list(Total=S_length_age$AGE),by=list(YEAR=S_length_age$YEAR),FUN=length)
 Total_LENGTH<-merge(Total,Total_LENGTH, all=T)
 Total_LENGTH$PROB<-Total_LENGTH$FREQ/Total_LENGTH$Total
 Total_LENGTH<-subset(Total_LENGTH,select=-c(Total))
 Total_LENGTH<-Total_LENGTH[order(Total_LENGTH$YEAR,Total_LENGTH$SEX,Total_LENGTH$LENGTH),]
 remove(Total,GRID)



## Survey WEIGHT at Age data
source("SWEIGHT.r")
S_WEIGHT<-SWEIGHT(data=Age_w,fyear=2011)
Total_AGE<-merge(Total_AGE,S_WEIGHT[[1]])
Total_AGE<-Total_AGE[order(Total_AGE$YEAR,Total_AGE$SEX,Total_AGE$AGE),]

## Survey Biomass Data

BIOMASS=sqlQuery(AFSC,"SELECT AI.BIOMASS_INPFC.SURVEY,
  AI.BIOMASS_INPFC.YEAR,
  AI.BIOMASS_INPFC.SUMMARY_AREA,
  AI.BIOMASS_INPFC.SPECIES_CODE,
  AI.BIOMASS_INPFC.AREA_BIOMASS,
  AI.BIOMASS_INPFC.BIOMASS_VAR,
  AI.BIOMASS_INPFC.AREA_POP,
  AI.BIOMASS_INPFC.POP_VAR
FROM AI.BIOMASS_INPFC
WHERE AI.BIOMASS_INPFC.SUMMARY_AREA != 799
AND AI.BIOMASS_INPFC.SPECIES_CODE    = 21740
ORDER BY AI.BIOMASS_INPFC.YEAR,
  AI.BIOMASS_INPFC.SUMMARY_AREA ")

BIO<-data.frame(aggregate(list(BIO=BIOMASS$AREA_BIOMASS,POP=BIOMASS$AREA_POP),by=list(YEAR=BIOMASS$YEAR),FUN=sum)  )
VAR<-data.frame(aggregate(list(BIO=BIOMASS$BIOMASS_VAR,POP=BIOMASS$POP_VAR),by=list(YEAR=BIOMASS$YEAR),FUN=sum)  )
survey_biomass=data.frame(YEAR=BIO$YEAR,BIO=BIO$BIO,B_STDEV=sqrt(VAR$BIO))


## calculate numbers at age

BIO<-subset(BIO,YEAR!=2010)
Total_WB<-merge(Total_AGE,BIO,all=T)
Total_WB<-Total_WB[order(Total_WB$YEAR,Total_WB$AGE),]
Total_WB$N_AGE<-Total_WB$PROB*Total_WB$POP
Total_WB$W_AGE<-Total_WB$N_AGE*Total_WB$WEIGHT/1000000

NAGE=aggregate(list(NUMBER=Total_WB$N_AGE),by=list(AGE=Total_WB$AGE,YEAR=Total_WB$YEAR,POP=Total_WB$POP),FUN=sum)
NAGE$PROP<- NAGE$NUMBER/NAGE$POP
NAGE<-subset(NAGE,select=-c(POP))


S_NAGE<-data.frame(matrix(ncol=15,nrow=7))
  names(S_NAGE)<-c("YEAR","A2","A3" ,"A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15")
  S_NAGE$YEAR<-c(1991,1994,1997,2000,2002,2004,2006)

 for( i in 1:length(S_NAGE$YEAR)){
  S_NAGE[i,2:15]<-NAGE$PROP[NAGE$YEAR==S_NAGE$YEAR[i]][2:15]
  }


##FISHEY
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

odbcClose(AFSC)

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



## population weight at age
test<-gam(log(WEIGHT)~s(AGE1),gamma=1.4,data=subset(Age_w,Age_w$YEAR>2002))
PP_WEIGHT_AGE<-exp(predict(test,newdata=data.frame(AGE1=c(2:15))))


## Data set
  ##Start YEAR
  SY<-1978
  ##END YEAR
  EY<-2011
  ## First Age class
  FAC<-2
  ##Last Age Class
  LAC<-15
  ## Number of Fisheries
  NF<-1
  ##Number of Surveys
  NS<-1
  ## Month of spawning
  MSPWN<-3
  ## Fishery and survey names
  F1_name<-noquote("Aleutians")                                                     ## Fishery Name
  S1_name<-noquote("NMFS_summer_bottom_trawl")                                      ## Survey Name


 ## other data
Catch_CV<-rep(0.05,length(unique(CATCH[[1]]$YEAR) ))                                 ## Fishery CV
AERR<-read.csv("AGEING_ERROR.csv",header=F)                                          ## Aging error matrix
GMat<- read.csv("GOA_MAT.csv",header=F)                                              ## GOA maturity 1983-2006
F_MN_SIZE<-c(177,103,131,99,670,125,288,155,220,269,159,75,84,187,500,500,500)       ## Fishery Multinomial sample size
S_MN_SIZE<-c(1,100,100,100,100,100,100)                                              ## Survey Multinomial sample size



##CATCH[[1]]$TONS/1000         ## CATCH BIOMASS
#F_NAGE[,2:15]                  ## fishery numbers at age
#round(F_WEIGHT$AGE_WEIGHT[,2:15])   ## fishery WEIGTH at age

#survey_biomass$BIO/1000      ## SURVEY BIOMASS
#survey_biomass$B_STDEV  ## SURVEY BIOMASS STDEV
#S_NAGE[,2:15]                 ## Survey numbers at age
#round(S_WEIGHT$AGE_WEIGHT[,2:15])   ## SURVEY WEIGTH at age

source("write_dat.r")
write_dat("ai11_AEC.dat")