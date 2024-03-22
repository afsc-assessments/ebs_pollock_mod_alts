plot(FWEIGHT_B$WEIGHT2[FWEIGHT_B$YEAR==1978]~FWEIGHT_B$AGE[FWEIGHT_B$YEAR==1978],type="l",ylim=c(0,2500))
for(i in 1979:2012){
 points(FWEIGHT_B$WEIGHT2[FWEIGHT_B$YEAR==i]~FWEIGHT_B$AGE[FWEIGHT_B$YEAR==i],col=(i-1978),type="l")
 }
 
 
plot(FWEIGHT_A$WEIGHT~FWEIGHT_A$AGE,col=as.numeric(FWEIGHT_B$PERIOD))
 #plot(FWEIGHT_A$pred[FWEIGHT_A$YEAR==1978]~FWEIGHT_A$AGE[FWEIGHT_A$YEAR==1978],type="l",ylim=c(0,2500))
for(i in 1978:2012){
 points(FWEIGHT_A$pred[FWEIGHT_A$YEAR==i]~FWEIGHT_A$AGE[FWEIGHT_A$YEAR==i],col=as.numeric(FWEIGHT_A$PERIOD[FWEIGHT_A$YEAR==i]),type="l")
 }
 
 for(i in 1978:2012){
 points(FWEIGHT_B$WEIGHT2[FWEIGHT_B$YEAR==i]~FWEIGHT_B$AGE[FWEIGHT_B$YEAR==i],col=as.numeric(FWEIGHT_B$PERIOD[FWEIGHT_B$YEAR==i]),type="l",lwd=1)
 }
 
 
 points(data2$WEIGHT~data2$AGE,col=as.numeric(data2$PERIOD)cex=0.5,pch=1)

