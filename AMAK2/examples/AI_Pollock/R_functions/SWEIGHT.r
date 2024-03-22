  SWEIGHT<-function(data2=Age_w,TA=Total_AGE,fyear=fyear){
  require(mgcv)
  #data=subset(data,!is.na(data$WEIGHT))
  ##fill in the 1980 weights by length
  #test<-gam(log(WEIGHT)~as.factor(YEAR)+s(LENGTH,by=as.factor(YEAR)),gamma=1.4,data=data)
  #new_data<-subset(data,data$YEAR==1980)
  #new_data$YEAR=1983
  #data$WEIGHT[data$YEAR==1980]<-exp(predict(test,newdata=new_data))

  Weight<-subset(data2,!is.na(data2$WEIGHT))
  Weight$AGE1<-Weight$AGE
   Weight$AGE1[Weight$AGE>15]<-15
  

  WEIGHT<-aggregate(list(WEIGHT=Weight$WEIGHT),by=list(SEX=Weight$SEX,AGE=Weight$AGE1,YEAR=Weight$YEAR),FUN=mean)
  grid=expand.grid(YEAR=sort(unique(WEIGHT$YEAR)),SEX=c(1,2),AGE=c(min(Weight$AGE1):15))
  WEIGHT_A<-merge(grid,WEIGHT,all=T)
  test<-gam(log(WEIGHT)~as.factor(YEAR)+as.factor(SEX)+s(AGE,by=as.factor(YEAR)),gamma=1.4,data=subset(WEIGHT_A,is.na(WEIGHT_A$AGE)==F))
  WEIGHT_A$pred<-exp(predict(test,newdata=WEIGHT_A))
  WEIGHT_A$WEIGHT2<-WEIGHT_A$WEIGHT
  WEIGHT_A$WEIGHT2[is.na(WEIGHT_A$WEIGHT2)==T]<-WEIGHT_A$pred[is.na(WEIGHT_A$WEIGHT2)==T]

  WEIGHT_A2<-subset(WEIGHT_A,select=c(YEAR,SEX,AGE,WEIGHT2))

  Total_A<-merge(TA,WEIGHT_A2,all=T)
  Total_A<-Total_A[order(Total_A$YEAR,Total_A$SEX,Total_A$AGE),]

  Total<-subset(Total_A,!is.na(Total_A$WEIGHT2))
  Total$FREQ[is.na(Total$FREQ)]<-0
  Total$w<-Total$WEIGHT*Total$FREQ
  Total$w[Total$w==0]<-Total$WEIGHT2[Total$w==0]
  
  Total<-aggregate(list(Total=Total$FREQ,WW=Total$w),by=list(YEAR=Total$YEAR,AGE=Total$AGE),FUN=sum)
  Total$Total[Total$Total==0]<-2
  Total$AWEIGHT<-Total$WW/Total$Total
  Total<-Total[order(Total$YEAR,Total$AGE),]
  Total<-subset(Total,select=-c(Total,WW))

  GRID<-expand.grid(YEAR=c(1978:fyear),AGE=c(1:15))
  WEIGHT_B<-merge(GRID,Total,all=T)

  WEIGHT_B$PERIOD<-"F1"
  WEIGHT_B$PERIOD[WEIGHT_B$YEAR>1985]<-"F2"
  WEIGHT_B$PERIOD[WEIGHT_B$YEAR>1990]<-"D1"
  WEIGHT_B$PERIOD[WEIGHT_B$YEAR>1995]<-"D2"
  WEIGHT_B$PERIOD[WEIGHT_B$YEAR>2000]<-"D3"
  
  
  WEIGHT_B$PERIOD<-as.factor(WEIGHT_B$PERIOD)

  test1<-gam(log(AWEIGHT)~PERIOD+s(AGE,by=as.factor(PERIOD)),gamma=1.4,data=subset(WEIGHT_B,!is.na(WEIGHT_B$AGE)))

    
  WEIGHT_B$pred[!is.na(WEIGHT_B$AGE)]<-exp(predict(test1,newdata=WEIGHT_B))
  
  WEIGHT_B$AWEIGHT2[is.na(WEIGHT_B$AWEIGHT)]<-WEIGHT_B$pred[is.na(WEIGHT_B$AWEIGHT)]
  
  WEIGHT_B$AWEIGHT2[is.na(WEIGHT_B$AWEIGHT2)]<- WEIGHT_B$AWEIGHT[is.na(WEIGHT_B$AWEIGHT2)]


  S_AGE_WEIGHT<-data.frame(matrix(ncol=15,nrow=length(unique(WEIGHT_B$YEAR))))
  names(S_AGE_WEIGHT)<-c("YEAR","A2","A3" ,"A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15")
  S_AGE_WEIGHT$YEAR<-c(1978:fyear)

  Years1<-c(1978:fyear)
  for( i in 1:length(Years1)){
  S_AGE_WEIGHT[i,2:15]<-WEIGHT_B$AWEIGHT2[WEIGHT_B$YEAR==Years1[i]][2:15]
  }
  WEIGHT2<-subset(WEIGHT_A,select=c(YEAR,SEX,AGE,WEIGHT2))

  SWEIGHT<-list( WEIGHT=WEIGHT2,AGE_WEIGHT=S_AGE_WEIGHT)
  SWEIGHT
}