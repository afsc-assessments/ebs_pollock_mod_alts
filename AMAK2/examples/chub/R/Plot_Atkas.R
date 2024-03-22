Plot_Phase<- function (mod=mod1,main="Model 1",ylim=c(0,1.5),xlim=c(0,10)){
  plot(mod$msy_mt[,12],mod$msy_mt[,5],typ="l",main=main,ylim=ylim,xlim=xlim,xlab="B/Bmsy",ylab="F/Fmsy")
  text(mod$msy_mt[,12],mod$msy_mt[,5],mod$msy_mt[,1]-100*floor(mod$msy_mt[,1]/100))
  lastone=length(mod$msy_mt[,12])
  text(mod$msy_mt[lastone,12],mod$msy_mt[lastone,5],mod$msy_mt[lastone,1]-100*floor(mod$msy_mt[lastone,1]/100),col="blue",cex=1.4)
  abline(v=1,h=1,col="red",lty=2)
}

Plot_Fspr2<- function (mod=mod1,mod2=mod2,ylim=c(0,1)){
  plot(mod$msy_mt[,1],1-mod$msy_mt[,2],lwd=2,col="red",typ="l",ylim=ylim,xlab="Year",ylab="1-SPR%")
  lines(mod2$msy_mt[,1],1-mod2$msy_mt[,2],lwd=2,col="blue",lty=2,wid=2)
  legend("topleft",c("Model 1","Model 2"),lty=c(1,2),col=c("red","blue"),inset=.02)
  
}

Plot_Fspr<- function (mod=mod1,main="Model 1",ylim=c(0,1)){
  plot(mod$msy_mt[,1],1-mod$msy_mt[,2],lwd=2,col="blue",typ="l",main=main,ylim=ylim,xlab="Year",ylab="1-SPR%")
  abline(h=.65,col="red",lty=2)
}
