Mntns_srv<-function(dat,Title)    {
  # windows(,9,6.5)
  par(mfrow=c(1, length(dat$Fshry_names)),oma=c(0,1,3,1) )
  for(i in 1:length(dat$Fshry_names))
  {
  	sel.age.mountain(dat, f=i, typ='S', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2.5, new=F, 
  	     cex.yax=1.)
  }
  mtext(Title,3,0,outer=T, cex=2, col="blue",family="sans",font=2)
}
Mntns<-function(dat,Title)    {
  # windows(,9,6.5)
  par(mfrow=c(1, length(dat$Fshry_names)),oma=c(0,1,3,1) )
  for(i in 1:length(dat$Fshry_names))
  {
  	sel.age.mountain(dat, f=i, typ='F', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2.5, new=F, 
  	     cex.yax=1.)
  }
  mtext(Title,3,00,outer=T, cex=2, col="blue",family="sans",font=2)
}
