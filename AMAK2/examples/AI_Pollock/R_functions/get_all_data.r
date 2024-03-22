get_all_data<-function(fyear=2013, species=21740,region="'AI'",area="AI",MINAGE=100,filename="ai13_Nov.dat",WRITE=TRUE)  {

  ## Survey Age Data
  Age1=get_Sage(region=region,species=species,sy=1977,GEAR=172)
  Age2=get_Sage(region=region,species=species,sy=1977,GEAR=160)
  Age<-rbind(Age2,Age1)

  
  Age<-subset(Age,Age$YEAR<=fyear)
  Age_w<-subset(Age,Age$END_LONGITUDE > 0 | Age$END_LONGITUDE < -170)
  Age_w<-subset(Age_w,is.na(Age_w$AGE)==F)
  Age_w$AGE1<-Age_w$AGE
  Age_w$AGE1[Age_w$AGE1>=15]=15

##Survey Length Data
  Length1=get_Slength(region=region,species=21740,sy=1977,GEAR=172)
  Length2=get_Slength(region=region,species=21740,sy=1977,GEAR=160)
  Length<-rbind(Length2,Length1)
  
  Length<-subset(Length,Length$YEAR<=fyear)
  Length_w<-subset(Length,Length$END_LONGITUDE > 0 | Length$END_LONGITUDE < -170)
  Length_w<-subset(Length_w,Length_w$YEAR!=1982)
  Length_w<-subset(Length_w,Length_w$SEX!=3)


## Survey WEIGHT at Age
  TOTAL_AGE<-AGE_LENGTH(adata=Age_w,ldata=Length_w,fyear=fyear)  ##
  
## Survey Biomass Data
  BIOMASS<-get_AIbiom(species=21740)

  BIOMASS<-subset(BIOMASS,BIOMASS$YEAR<=fyear)

  BIO<-data.frame(aggregate(list(BIO=BIOMASS$AREA_BIOMASS,POP=BIOMASS$AREA_POP),by=list(YEAR=BIOMASS$YEAR),FUN=sum)  )
  VAR<-data.frame(aggregate(list(BIO=BIOMASS$BIOMASS_VAR,POP=BIOMASS$POP_VAR),by=list(YEAR=BIOMASS$YEAR),FUN=sum)  )
  survey_biomass=data.frame(YEAR=BIO$YEAR,BIO=BIO$BIO,B_STDEV=sqrt(VAR$BIO))

## calculate numbers at age

##BIO<-subset(BIO,YEAR!=2010)
  Total_WB<-merge(TOTAL_AGE$Total_AGE,BIO,all=T)
  Total_WB<-Total_WB[order(Total_WB$YEAR,Total_WB$AGE),]
  Total_WB$N_AGE<-Total_WB$PROB*Total_WB$POP
  Total_WB$W_AGE<-Total_WB$N_AGE*Total_WB$WEIGHT/1000000

  NAGE=aggregate(list(NUMBER=Total_WB$N_AGE),by=list(AGE=Total_WB$AGE,YEAR=Total_WB$YEAR,POP=Total_WB$POP),FUN=sum)
  NAGE$PROP<- NAGE$NUMBER/NAGE$POP
  NAGE<-subset(NAGE,select=-c(POP))

  years<-sort(unique(Age_w$YEAR))
  S_NAGE<-data.frame(matrix(ncol=15,nrow=length(years)))
  names(S_NAGE)<-c("YEAR","A2","A3" ,"A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15")
  S_NAGE$YEAR<-years

 for( i in 1:length(years)){
  S_NAGE[i,2:15]<-NAGE$PROP[NAGE$YEAR==S_NAGE$YEAR[i]][2:15]
  }

  AI_DATA<-list(S_NAGE=S_NAGE)
##FISHEY
  FAge=get_Fage(species=201,sy=1977,area=area)

  FAge1<-subset(FAge,FAge$GENERIC_AREA==541&FAge$LATITUDE<5300)
  FAge2<-subset(FAge,FAge$GENERIC_AREA==542&FAge$LATITUDE<5230)
  FAge3<-subset(FAge,FAge$GENERIC_AREA==543&FAge$LATITUDE<5330)
  FAge<-rbind(FAge1,FAge2,FAge3)
  FAge<-FAge[order(FAge$YEAR,FAge$GENERIC_AREA),]
  remove(FAge1,FAge2,FAge3)
  FAge<-subset(FAge,select= -c(LATITUDE,LONGITUDE,E_W) )

  FLength= get_Flength(species=201,sy=1977,area=area)

  FLength1<-subset(FLength,FLength$GENERIC_AREA==541&FLength$LATITUDE<5300)
  FLength2<-subset(FLength,FLength$GENERIC_AREA==542&FLength$LATITUDE<5230)
  FLength3<-subset(FLength,FLength$GENERIC_AREA==543&FLength$LATITUDE<5330)
  FLength<-rbind(FLength1,FLength2,FLength3)
  FLength<-FLength[order(FLength$YEAR,FLength$GENERIC_AREA),]
  remove(FLength1,FLength2,FLength3)
  FLength<-subset(FLength,select= -c(LATITUDE,LONGITUDE,E_W) )

  DAge=get_Dage(species=201,area=area)
  
  DAge<-subset(DAge,DAge$YEAR<=fyear)
  DAge1<-subset(DAge,DAge$NMFS_AREA==541&DAge$LATDD_END<53.00)
  DAge2<-subset(DAge,DAge$NMFS_AREA==542&DAge$LATDD_END<52.50)
  DAge3<-subset(DAge,DAge$NMFS_AREA==543&DAge$LATDD_END<53.50)
  DAge<-rbind(DAge1,DAge2,DAge3)
  DAge<-DAge[order(DAge$YEAR,DAge$NMFS_AREA),]
  remove(DAge1,DAge2,DAge3)
  DAge<-subset(DAge,select= -c(LATDD_END,LONDD_END) )
  DAge<-subset(DAge,DAge$YEAR<1999)

  DLength=get_Dlength(species=201,area=area)

  DLength<-subset(DLength,DLength$YEAR<=fyear)
  DLength1<-subset(DLength,DLength$NMFS_AREA==541&DLength$LATDD_END<53.00)
  DLength2<-subset(DLength,DLength$NMFS_AREA==542&DLength$LATDD_END<52.50)
  DLength3<-subset(DLength,DLength$NMFS_AREA==543&DLength$LATDD_END<53.50)
  DLength<-rbind(DLength1,DLength2,DLength3)
  DLength<-DLength[order(DLength$YEAR,DLength$NMFS_AREA),]
  remove(DLength1,DLength2,DLength3)
  DLength<-subset(DLength,select= -c(LATDD_END,LONDD_END) )
  names(FLength)<-names(DLength)
  names(FAge)<-names(DAge)

  FAge<-rbind(FAge,DAge)
  FLength<-rbind(FLength,DLength)

  if(fyear  > 2006){
   AAge<-read.csv("Acoustic_survey_data.csv",header=T)
   AAge<-subset(AAge,AAge$YEAR<=fyear)

   AAge$WEIGHT<-AAge$WEIGHT/1000
   ALength<-aggregate(list(FREQ=AAge$LENGTH),by=list(YEAR=AAge$YEAR,NMFS_AREA=AAge$NMFS_AREA,SPECIES=AAge$SPECIES,SEX=AAge$SEX,LENGTH=AAge$LENGTH),FUN=length)
   ALength<-ALength[order(ALength$YEAR,ALength$LENGTH),]
   AAge$GEAR<-1
   ALength$GEAR<-1
   names(AAge)<-names(FAge)
   names(ALength)<-names(FLength)
   FAge<-rbind(FAge,AAge)
   FLength<-rbind(FLength,ALength)
   }

## Fishery AGE
  FAge_w<-subset(FAge,!is.na(FAge$AGE))
  FAge_w<-subset(FAge_w,!is.na(FAge_w$LENGTH))

  ## removing years with too few datapoints
  x<-aggregate(list(NUMBER=FAge_w$AGE),by=list(YEAR=FAge_w$YEAR),FUN=sum)
  x<-subset(x,x$NUMBER >=MINAGE)
  x<-unique(x$YEAR)
  FAge_w<-subset(FAge_w,FAge_w$YEAR%in%x)
  FLength_w<-subset(FLength,FLength$YEAR%in%x)
  remove(x)


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
  
  F_Total_AGE<-AGE_LENGTH(adata=FAge_w,ldata=FLength_w,fyear=fyear,type="F")

## Fishery Catch Data
  source("get_catch.r")
  
  GET_CATCH(area="AI",species="'PLCK'",FYR=fyear,ADD_OLD=TRUE,OLD_FILE="OLD_CATCH.csv")
  
  
  CATCH<-GET_CATCH(FYR=fyear)
  CATCH<-subset(CATCH,CATCH$YEAR<=fyear)


  F_Total_WB<-merge(F_Total_AGE$Total_AGE,CATCH,all.x=T)
  F_Total_WB<-F_Total_WB[order(F_Total_WB$YEAR,F_Total_WB$AGE),]


  F_Total_WB$N_AGE<-F_Total_WB$PROB*F_Total_WB$TONS/F_Total_WB$WEIGHT2*1000000
  F_Total_WB$W_AGE<-F_Total_WB$PROB*F_Total_WB$TONS

  F_total<-aggregate(list(FREQ=F_Total_AGE$Total_AGE$FREQ),by=list(YEAR=F_Total_AGE$Total_AGE$YEAR,AGE=F_Total_AGE$Total_AGE$AGE),FUN=sum)
  F_number<-aggregate(list(NUMBER=F_Total_AGE$Total_AGE$FREQ),by=list(YEAR=F_Total_AGE$Total_AGE$YEAR),FUN=sum)
  F_total<-merge(F_total,F_number)
  F_total$PROB<-F_total$FREQ/F_total$NUMBER
  F_total<-F_total[order(F_total$YEAR,F_total$AGE),]

  years<-sort(unique(F_total$YEAR))

  F_NAGE<-data.frame(matrix(ncol=15,nrow=length(years)))
  names(F_NAGE)<-c("YEAR","A2","A3" ,"A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15")
  F_NAGE$YEAR<-years

 for( i in 1:length(F_NAGE$YEAR)){
  F_NAGE[i,2:15]<-F_total$PROB[F_total$YEAR==F_NAGE$YEAR[i]][1:14]
  }



## population weight at age
  test<-gam(log(WEIGHT)~s(AGE1),gamma=1.4,data=subset(Age_w,Age_w$YEAR>2002))
  PP_WEIGHT_AGE<-exp(predict(test,newdata=data.frame(AGE1=c(2:15))))

## get fishery length sample size
SS=get_FSS()
ss<-rbind(get_FSS(),get_DSS())
ss<-aggregate(list(SS=ss$SS),by=list(YEAR=ss$YEAR),FUN=sum)
ss<-subset(ss,ss$YEAR %in% F_NAGE$YEAR)
ss$SS[ss$YEAR>=2006]<-100




## Data set
  ##Start YEAR
  AI_DATA<-list(SY=1978,EY=fyear,FAC=2,LAC=15,NF=1,
  F_NAME=noquote("Aleutians"),
  CATCH=CATCH,
  Catch_CV=rep(0.05,length(unique(CATCH$YEAR) )),
  F_NAGE=F_NAGE,
  F_MN_SIZE=round(ss$SS),
  F_WEIGHT=F_Total_AGE$F_WEIGHT,
  NS=1,
  S_name=noquote("NMFS_summer_bottom_trawl"),
  S_BIO=survey_biomass,
  S_MN_SIZE=c(1,1,1,100,100,100,100,100,100,100),
  S_WEIGHT=TOTAL_AGE$S_WEIGHT,
  S_NAGE=subset(S_NAGE,S_NAGE > 1980),
  PP_WEIGHT_AGE=PP_WEIGHT_AGE,
  GMat= read.csv("GOA_MAT.csv",header=F),
  MSPWN=3,
  AERR=read.csv("AGEING_ERROR.csv",header=F),
  FILENAME=filename
  )
 if(WRITE==TRUE){
  write_dat(data=AI_DATA,data_file=filename)
  }
 AI_DATA
}