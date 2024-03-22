FWEIGHT<-function(data2=FAge_w, sy=1978, fy=fyear){
  require(mgcv)
  
  data2$PERIOD<-"F"
  data2$PERIOD[data2$YEAR>1990]<-"D1"
  data2$PERIOD[data2$YEAR>1995]<-"D2"
  data2$PERIOD<-as.factor(data2$PERIOD)
  
   test<-gam(log(WEIGHT)~PERIOD+as.factor(YEAR)+s(AGE,by=PERIOD,YEAR),gamma=1.4,data=data2)
   test1<-gam(log(WEIGHT)~PERIOD+s(AGE,by=as.factor(PERIOD)),gamma=1.4,data=data2)
  
  
  FWeight<-subset(data2,!is.na(data2$WEIGHT))
  FWeight$AGE1<-FWeight$AGE
  FWeight$AGE1[FWeight$AGE>15]<-15
  FWEIGHT<-aggregate(list(WEIGHT=FWeight$WEIGHT),by=list(AGE=FWeight$AGE1,YEAR=FWeight$YEAR, PERIOD=FWeight$PERIOD),FUN=mean)
  FNUMBER<-aggregate(list(NUMBER=FWeight$WEIGHT),by=list(AGE=FWeight$AGE1,YEAR=FWeight$YEAR,PERIOD=FWeight$PERIOD),FUN=length)
  FWEIGHT$NUMBER<-FNUMBER$NUMBER
  grid=expand.grid(YEAR=sort(unique(FWEIGHT$YEAR)),AGE=c(min(data2$AGE):15))
  FWEIGHT_A<-merge(grid,FWEIGHT,all=T)
  FWEIGHT_A$NUMBER[is.na(FWEIGHT_A$NUMBER)]<-0
  FWEIGHT_A$PERIOD<-"F"
  FWEIGHT_A$PERIOD[FWEIGHT_A$YEAR>1990]<-"D1"
  FWEIGHT_A$PERIOD[FWEIGHT_A$YEAR>1995]<-"D2"
  FWEIGHT_A$PERIOD<-as.factor(FWEIGHT_A$PERIOD)
 
  #test<-gam(log(WEIGHT)~PERIOD+as.factor(YEAR)+s(AGE,by=PERIOD,YEAR),gamma=1.4,data=subset(FWEIGHT_A,is.na(FWEIGHT_A$AGE)==F))
  

  FWEIGHT_A$pred<-exp(predict(test,newdata=FWEIGHT_A))
  FWEIGHT_A$WEIGHT2[FWEIGHT_A$NUMBER>=20]<-FWEIGHT_A$WEIGHT[FWEIGHT_A$NUMBER>=20]
  FWEIGHT_A$WEIGHT2<-FWEIGHT_A$pred
  #FWEIGHT_A$WEIGHT2<-FWEIGHT_A$WEIGHT
  #FWEIGHT_A$WEIGHT2[is.na(FWEIGHT_A$WEIGHT2)==T]<-FWEIGHT_A$pred[is.na(FWEIGHT_A$WEIGHT2)==T]



  GRID<-expand.grid(YEAR=c(sy:fy),AGE=c(2:15))
  FWEIGHT_B<-merge(GRID,FWEIGHT_A,all=T)
  FWEIGHT_B$PERIOD<-"F"
  FWEIGHT_B$PERIOD[FWEIGHT_B$YEAR>1990]<-"D1"
  FWEIGHT_B$PERIOD[FWEIGHT_B$YEAR>1995]<-"D2"
  FWEIGHT_B$PERIOD<-as.factor(FWEIGHT_B$PERIOD)


  FWEIGHT_B$pred<-exp(predict(test1,newdata=FWEIGHT_B))
  #FWEIGHT_B$WEIGHT2[is.na(FWEIGHT_B$WEIGHT)==F]<-FWEIGHT_B$WEIGHT[is.na(FWEIGHT_B$WEIGHT)==F]
  FWEIGHT_B$WEIGHT2[is.na(FWEIGHT_B$WEIGHT2)==T]<-FWEIGHT_B$pred[is.na(FWEIGHT_B$WEIGHT2)==T]
  


  F_AGE_WEIGHT<-data.frame(matrix(ncol=15,nrow=length(unique(FWEIGHT_B$YEAR))))
  names(F_AGE_WEIGHT)<-c("YEAR","A2","A3" ,"A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15")
  F_AGE_WEIGHT$YEAR<-c(sy:fy)

  Years1<-c(sy:fy)
  for( i in 1:length(Years1)){
  F_AGE_WEIGHT[i,2:15]<-FWEIGHT_B$WEIGHT2[FWEIGHT_B$YEAR==Years1[i]][1:14]
  }
  
  FWEIGHT<-list( FWEIGHT=FWEIGHT_B,AGE_WEIGHT=F_AGE_WEIGHT)
  FWEIGHT
}