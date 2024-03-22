length2<-subset(length, length$STRATUM<700)


length3<-aggregate(list(FREQ=length2$FREQUENCY),by=list(YEAR=length2$YEAR,LENGTH=length2$LENGTH),FUN=sum)
grid<-expand.grid(YEAR=c(2002,2004,2006,2010,2012),LENGTH=seq(1,100,1)

LENGTH=merge(grid,length3,all=T)

 LENGTH$FREQ[is.na(LENGTH$FREQ)==T]<-0
 
LENGTH$PROB[LENGTH$YEAR==2012]<-LENGTH$FREQ[LENGTH$YEAR==2012]/sum(LENGTH$FREQ[LENGTH$YEAR==2012])
 LENGTH$PROB[LENGTH$YEAR==2010]<-LENGTH$FREQ[LENGTH$YEAR==2010]/sum(LENGTH$FREQ[LENGTH$YEAR==2010])
 LENGTH$PROB[LENGTH$YEAR==2006]<-LENGTH$FREQ[LENGTH$YEAR==2006]/sum(LENGTH$FREQ[LENGTH$YEAR==2006])
 LENGTH$PROB[LENGTH$YEAR==2004]<-LENGTH$FREQ[LENGTH$YEAR==2004]/sum(LENGTH$FREQ[LENGTH$YEAR==2004])
 LENGTH$PROB[LENGTH$YEAR==2002]<-LENGTH$FREQ[LENGTH$YEAR==2002]/sum(LENGTH$FREQ[LENGTH$YEAR==2002])
 
 
 
HL12<-HL
HL12$counts<-LENGTH$FREQ[LENGTH$YEAR==2012]
HL12$density<-LENGTH$PROB[LENGTH$YEAR==2012]
plot(HL12,col="brown",freq=F,ylab="Proportion",xlab="Fork length(cm)",main=2012,cex.lab=1.2,cex.axis=1.2,cex.main=1.5)

HL10<-HL
HL10$counts<-LENGTH$FREQ[LENGTH$YEAR==2010]
HL10$density<-LENGTH$PROB[LENGTH$YEAR==2010]
plot(HL10,col="brown",freq=F,ylab="Proportion",xlab="Fork length(cm)",main=2010,cex.lab=1.2,cex.axis=1.2,cex.main=1.5)


HL06<-HL
HL06$counts<-LENGTH$FREQ[LENGTH$YEAR==2006]
HL06$density<-LENGTH$PROB[LENGTH$YEAR==2006]
plot(HL06,col="brown",freq=F,ylab="Proportion",xlab="Fork length(cm)",main=2006,cex.lab=1.2,cex.axis=1.2,cex.main=1.5)

HL04<-HL
HL04$counts<-LENGTH$FREQ[LENGTH$YEAR==2004]
HL04$density<-LENGTH$PROB[LENGTH$YEAR==2004]
plot(HL04,col="brown",freq=F,ylab="Proportion",xlab="Fork length(cm)",main=2004,cex.lab=1.2,cex.axis=1.2,cex.main=1.5)

HL02<-HL
HL02$counts<-LENGTH$FREQ[LENGTH$YEAR==2002]
HL02$density<-LENGTH$PROB[LENGTH$YEAR==2002]
plot(HL02,col="brown",freq=F,ylab="Proportion",xlab="Fork length(cm)",main=2002,cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
