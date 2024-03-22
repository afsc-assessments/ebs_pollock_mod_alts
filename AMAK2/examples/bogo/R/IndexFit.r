IndexFit<-function(dat,f=1,main="Model 1",ylab="", yf=1970,yl=2014 ,Doleg=F,ylim=c(0,1))
{
  attach(dat)
  require(plotrix)
  x<-paste("Obs_Survey_",f,sep="")
  Data<-get(x)
  d1<-subset(Data,Data[,2]>0)
#  windows(8,8)
  plotCI( d1[,1], d1[,2], 2*d1[,4],
  new=F,
  main=main,
          #Index_names[f],
  cex.main=1.3,
  cex.lab=1.2,
  ylim=ylim, #c(0,(max(d1[,2])*1.9)),
  xlim=c(yf,yl),pch=19,cex=1,
  ylab=ylab,
  xlab="Year",col="blue",scol="black",lab=c(10,10,7)
)

  points(Data[,3]~Data[,1],type="l",col="red",lwd=2)
  points(Data[Data[,2]>0,1],Data[Data[,2]>0,2],type="p",pch=19,col="blue",cex=1)

  if (Doleg){
    text(min(Data[,1])+2,(max(d1[,2])*0.5),"Observed estimate",pos=4,cex=0.8)
    text(min(Data[,1])+2,(max(d1[,2])*0.4),"Model Fit",pos=4,cex=0.8)
    points(min(Data[,1])+2,(max(d1[,2])*0.5),pch=19,col="blue",cex=0.8)
    text(min(Data[,1])+2,(max(d1[,2])*0.4),"--",cex=1.5,col="red")
  }
  detach(dat)
}

Indices<-function(dat=mod1c,Title="Model 1c"){
  # Index fits
  # windows(,11, 6.5)
  par(mfrow=c(2,4), oma=c(0,1,3,1) )
  for(i in 1:length(dat$Index_names))
  {
    if (i==1) 
      Doleg=T
    else 
      Doleg=F
  
   IndexFit(dat,f=i,Doleg=Doleg)
  }
  # text(2000,  20, Title , cex=3, col="black",family="sans serif",font=1)
  mtext(Title,3,-0,outer=T, cex=2, col="blue",family="sans",font=2)
}

CatchFit<-function(dat=am1,Title="Model 1"){
 par(mfrow=c(2,2), oma=c(0,1,3,1) )
 for (i in 1:4){
   p.catch.fit(dat,f=i)
 }
 mtext(Title,3,-0,outer=T, cex=2, col="blue",family="sans",font=2)
 }
