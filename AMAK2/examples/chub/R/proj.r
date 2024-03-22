plt_proj<-function(dat,Title="Model 1")    {
  #windows(,6.5,9)
  FutBiom(dat)
  #windows(,6.5,9)
  FutCatch(dat)
  # mtext(Title,3,0,outer=T, cex=2, col="blue",family="sans",font=2)
}

FutBiom <- function (dat, new=TRUE,ylab="SSB, both sexes (kt)")
{
  ltmp=c("F35","F=SQ","F=.5SQ","F=.25SQ", "F=0.0")
  attach(dat)
  Data = get("SSB")  
  xlim = c(min(Data[,1]+12),max(Data[,1])+10)
  ymax = 1.3*max(Data[13:length(Data[,2]),2])
  ylim = c(0,ymax)
  plot(Data[,1],Data[,2],type="n",lwd=2,
  xlab="Year",ylab=ylab,
  #lab=c(10,10,7),
  xlim=xlim, ylim=ylim, col="red")
  n=length(Data[,1])
  x=c(Data[13:n,1],Data[length(Data[,1]):13,1])
  y=c(Data[13:n,4],Data[length(Data[,1]):13,5])
  polygon(x,y,col="salmon",lty=3,border="grey")
  lines(Data[13:n,1],Data[13:n,2],lwd=2, col="darkred")
  Data = get("SSB_fut_5")  
  x=c(Data[,1],Data[length(Data[,1]):1,1])
  y=c(Data[,4],Data[length(Data[,1]):1,5])
  polygon(x, y,col=rainbow(1,alpha=.2,start=.8),
          lty=3,border="grey")
  Data = get("SSB_fut_2")  
  x=c(Data[,1],Data[length(Data[,1]):1,1])
  y=c(Data[,4],Data[length(Data[,1]):1,5])
  polygon(x, y,col=rainbow(1,alpha=.5,start=.12),
          lty=3,border="grey")
  for (i in 1:5)
  {
    Data = get(paste("SSB_fut_",i,sep=""))  
    lines(Data[,1],Data[,2],lwd=2, col=i)
    points(Data[,1],Data[,2],pch=i, col=i)
  }

  # legend("topright", ltmp, col=1:5,lwd=2, title="          ", inset = .02)
  legend("topright", ltmp, pch=1:5,lwd=2,col=1:5,       title="Future SSB", inset = .02)
  detach(dat)
}

FutCatch <- function (dat, new=TRUE,ylab="Total Catch (kt)",fy=1977,ly=2020)
{
  require(lattice)
  attach(dat)
 nfsh = length( dat$Fshry_names)
 catch_tot = Obs_catch_1
 for (i in 1:nfsh) catch_tot= catch_tot+get(paste("Obs_catch_",i,sep=""))
  xlim = c(fy,ly)
  ylim = c(0,1.5*max(catch_tot))
  plot(fy:(fy+length(catch_tot)-1),catch_tot,
    type="l",lwd=2,
    xlab="Year",ylab=ylab,xlim=xlim, ylim=ylim, col="red")
  x= Catch_fut_1[,1]
  y= Catch_fut_1[,2]
  ltmp=c("F35","F=SQ","F=.5SQ","F=.25SQ")
  print(ltmp)
  for (i in 1:4)
  {
    y = get(paste("Catch_fut_",i,sep=""))
    #lines(x,y[,2],lty=i,lwd=3)
    points(x,y[,2],pch=i      )
  }
  # legend("topright", ltmp, lty=1:4,lwd=3, title="Future catch", inset = .02)
  legend("topright", ltmp, pch=1:4,       title="Future catch", inset = .02)
  detach(dat)
}
