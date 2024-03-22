
AGE_LENGTH<-function(adata=Age_w,ldata=Length_w,fyear=2012,type="S"){
  years<-sort(unique(adata$YEAR))
  z<-data.frame(matrix(ncol=4,nrow=1))
  names(z)<-c("YEAR","SEX","AGE","LENGTH" )
  for( i in 1:length(years)){
    x=find_AL(adata,ldata,years[i])
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

  if(type=="S"){
  S_WEIGHT<-SWEIGHT(data2=adata,TA=Total_AGE,fyear=fyear)
  Total_AGE<-merge(Total_AGE,S_WEIGHT[[1]])
  Total_AGE<-Total_AGE[order(Total_AGE$YEAR,Total_AGE$SEX,Total_AGE$AGE),]
  TOTAL_AGE<-list(Total_AGE=Total_AGE,S_WEIGHT=S_WEIGHT)

  }
  
  if(type=="F"){
  S_WEIGHT=FWEIGHT(data2=adata, sy=1978, fy=fyear)
  Total_AGE<-merge(Total_AGE,S_WEIGHT[[1]])
  Total_AGE<-Total_AGE[order(Total_AGE$YEAR,Total_AGE$SEX,Total_AGE$AGE),]
  TOTAL_AGE<-list(Total_AGE=Total_AGE,F_WEIGHT=S_WEIGHT)

  }
   TOTAL_AGE
 }
