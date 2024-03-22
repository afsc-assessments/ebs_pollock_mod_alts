##biomass plots with polygons
## this function can handle the results of up to 2 separate models
# Still to do: 
#	make it able to handle more models (how many?)
#	make polygons transparent
p.biom.pol <- function (dat1, dat2, n.mod=1, typ="SSB",color=c("light blue", 'darkblue', 'light green', 'darkgreen'), new=TRUE)
{
  if(typ=="SSB")titl<-"Spawning Biomass (t)"
  if(typ=="TB"){typ="TotBiom"
                titl<-"Total Biomass (t)"
                }
  if(typ=="R"){typ="R"
                titl<-"Age 1 Recruites(N)"
                            }
  require(lattice)
  
  attach(dat1)
  Data = get(typ)  
  xlim1 = c(min(Data[,1]),max(Data[,1]))
  ymax1 = max(Data[,2])
  if(n.mod==1){
  	xlim = xlim1
  	ylim = c(0, ymax1*1.5)
  }
  detach(dat1)
  if(n.mod!=1){
  	attach(dat2)
	Data = get(typ)  
 	xlim2 = c(min(Data[,1]),max(Data[,1]))
 	ymax2 = max(Data[,2])
  	detach(dat2)
  	xlim = range(c(xlim1, xlim2))
  	ylim = c(0, max(ymax1, ymax2)*1.5)
  }
 
  attach(dat1)
  Data=get(typ)
  ##if(new==TRUE) windows(,6,6)
  plot(Data[,1],Data[,2],type="n",lwd=2,
  xlab="Year",ylab=paste(titl),lab=c(10,10,7),xlim=xlim, ylim=ylim, col=color[2])
  x=c(Data[,1],Data[length(Data[,1]):1,1])
  y=c(Data[,4],Data[length(Data[,1]):1,5])

  polygon(x,y,col=color[1],lty=3,border="grey")
  lines(Data[,1],Data[,2],lwd=2, col=color[2])
  detach(dat1)
  
  if(n.mod!=1){
  	attach(dat2)
  Data=get(typ)
  x=c(Data[,1],Data[length(Data[,1]):1,1])
  y=c(Data[,4],Data[length(Data[,1]):1,5])

  polygon(x,y,col=color[3],lty=3,border="grey")
  lines(Data[,1],Data[,2], lwd=2, col=color[4])
  detach(dat2)
  
  attach(dat1)
  Data=get(typ)
  lines(Data[,1],Data[,2], lwd=2, col=color[2])
  detach(dat1)	
  }
  
}

mountains <-
  function(zmat, xvec=NULL, yvec=NULL, zscale=3, nshades=100,
           xaxs='i', yaxs='i', xlab="", ylab="", las=1, addbox=FALSE, cex.xax=1, cex.yax=1,...)
{
  ## DESCRIPTION:
  # a function by Ian Taylor designed to look like the cool-looking Figure 7 in
  # Butterworth D.S., Ianelli J.N., Hilborn R. (2003) A statistical model for
  # stock assessment of southern bluefin tuna with temporal changes in selectivity.
  # South African Journal of Marine Science 25:331-362.

  errors <- FALSE
  for(icol in 1:ncol(zmat)){
    if(!is.numeric(zmat[,icol])){
      errors <- TRUE
      print(paste("error: column",icol,"of zmat is not numeric"))
    }
  }
  if(errors) return(invisible())

  # fill in vectors if not provided
  nrowz <- nrow(zmat)
  ncolz <- ncol(zmat)
  if(is.null(yvec)) yvec <- 1:nrowz
  if(is.null(xvec)) xvec <- 1:ncolz

  # define some limits
  xmin <- min(xvec)
  xmax <- max(xvec)
  zmax <- zscale*max(zmat)

  ny <- length(yvec)
  if(ny!=nrowz){
      print("length(yvec) must equal nrow(zmat)",quote=FALSE)
      return()
  }
  if(length(xvec)!=ncolz){
      print("length(xvec) must equal ncol(zmat)",quote=FALSE)
      return()
  }

  zseq <- seq(0, zmax, length=nshades)
  xvec2 <- c(xmin, xvec, xmax) # adding extra points for bottom corners of polygon

  # plot(0, type='n', xlim=c(xmin, xmax), ylim=c(0, 1.1*(ymax+ny)), xaxs='i', yaxs='i', ...)
  plot(0, type='n', xlim=c(xmin, xmax), ylim=c(min(yvec), (max(yvec) + 1.1*zmax)),
       xaxs=xaxs, yaxs=yaxs, xlab=xlab, ylab=ylab, axes=F, ...)

  for(iy in ny:1){
    zvec <- as.numeric(zmat[iy, ])
    zvec2 <- c(0, zscale*zvec, 0) # row from z matrix

    # calculate set of all intersections between polygon and the horizontal lines
    x3list <- list()
    for(iz in 1:nshades){
      z <- zseq[iz]
      x3 <- numeric()
      for(ix in 2:length(xvec2)){
          z1 <- zvec2[ix-1]
          z2 <- zvec2[ix]
          x1 <- xvec2[ix-1]
          x2 <- xvec2[ix]
          if(z >= min(z1, z2) & z < max(z1, z2)){
              x3 <- c(x3, (z-z1)*(x2-x1)/(z2-z1)+x1)
          }
      }
      x3list[[iz]] <- x3
    }
    # draw little polygons between each pair of horizontal lines
    for(iz in 2:length(x3list)){

      z2 <- zseq[iz]
      z1 <- zseq[iz-1]

      x3hi <- x3list[[iz]] # x-values of intersections along upper line
      x3lo <- x3list[[iz-1]]   # x-values of intersections along lower line

      npoly <- length(x3lo)/2

      for(ipoly in 1:npoly){
        xlo <- x3lo[ipoly*2 + -1:0] # lower line intersections for individual polygon
        xhi <- x3hi[x3hi>=xlo[1] & x3hi<=xlo[2]] # upper line intersections
        extra <- (zvec2 >= z1 & zvec2 <= z2 & xvec2>=xlo[1] & xvec2<=xlo[2]) # identifying extra points to add
        xhi2 <- c(xhi,xvec2[extra]) # adding extra points to vector of upper x-values
        zhi <- c(rep(z2,length(xhi)), zvec2[extra]) # add corresponding z-values
        zhi2 <- zhi[order(xhi2)] # put the z-values in order based on order of x-values
        xhi2 <- sort(xhi2) # now order the x-values

        # make polygon
        polygon(x = c(xlo[2:1],xhi2),
                y = yvec[iy]+ c(z1,z1,zhi2),
                col = grey(1-.9*z1/zmax),border=grey(1-.9*z1/zmax))
      }
    }
    # black polygon around the outside
    polygon(xvec2, yvec[iy]+zvec2)
  }

  # add axes
  axis(1,at=xvec, cex.axis=cex.xax)
  axis(2,at=yvec,las=las, cex.axis=cex.yax)
#  axis(4,at=yvec,las=las,labels=FALSE) # extra ticks on right hand side
  if(addbox) box() # add box if desired
}


## Fish Selectivity curve over ages and years, for fisheries or surveys
sel.age.mountain = function(dat, f=1, typ='F', xvec=NULL, yvec=NULL, zscale=3, nshades=100, xaxs='i', yaxs='i', xlab="", ylab="", las=1, new=TRUE, addbox=FALSE, cex.xax=1, cex.yax=1, ...)
{	
	attach(dat)
	if(typ=='F'){
		x = paste('sel_fsh_', f, sep='')
		mtitl=paste(Fshry_names[f], 'Selectivity')
	}
	if(typ=='S'){
		x = paste('sel_srv_', f, sep='')
		mtitl=paste(Index_names[f], 'Selectivity')
	}
	Select = get(x)
	zmat = Select[,-c(1:2)]
	yrs = Select[,2]
	#if(new==TRUE) windows(,5,6.5)
	mountains(zmat, xvec=xvec, yvec=yrs, zscale=zscale, nshades=nshades, xaxs=xaxs, yaxs=yaxs, xlab=xlab, ylab=ylab, las=1, addbox=addbox, cex.xax=cex.xax, cex.yax=cex.yax)
	title(main=mtitl)
	detach(dat)	
}
