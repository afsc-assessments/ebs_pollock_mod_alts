
############################################################################
#Fit of the model sbtmodxx.tpl to the data for jack mackerel
#Predicted length frequency for longline fishery 1 to the observed data.
#Required the output file _lab.rep.
#Requires library PBSmodelling, available from 
#Programmed by Trevor A. Branch 29 June 2009
#outputs from a _lab.rep file, assuming that the naming convention in the file has sections of 
#outputs separated by $name comments so that it can be easily parsed by readList() function of Jon Schnute 
#IMPORTANT: library PBSmodelling required, available at http://code.google.com/p/pbs-software/downloads/list
#then in R go to menu Packages, Install Package from Local zip file
#Run by:
###source("LengthsLL1.r")
##########################################################################
#LF_Fit <- function(labrep.file="r_report.rep", case_label="test1") {
   #library(PBSmodelling)  #for readList
   #library(gplots)        #for rich.colors()
   #subtle.color <- "gray40"
   #x <- readList(labrep.file)
#
   # length.list <- seq(from=x$lengths[1],by=x$lengths[2],length.out=x$lengths[3])
   #length.list <- seq(from=16,by=2,length.out=23)

   #obs.data <- x$p_len_obs  #extract the lengths where the fishery is =1, and exclude the fishery and year
   #pred.data <- x$p_len_est  
   # print(pred.data)
   years <- x$yrs+1900
   nyears <- length(x$yrs)
   nlengths <- length(length.list)
   mfcol <- c(ceiling(nyears/5),5)
   par(mfcol=mfcol,oma=c(3.8,4.5,3.5,1),mar=c(0,0,0,0))
   ylim <- c(0,1.05*max(obs.data,pred.data))
   cohort.color <- rich.colors(nlengths)   
   
   #these used for drawing outside line on histogram
   xvals <- vector(length=2*nlengths)
   yvals <- vector(length=2*nlengths)
   for (yr in 1:nyears) { 
      names.arg <- rep("",nlengths)
      #histogram lines and bars the same color and no gaps between them
      x <- barplot(obs.data[yr,],space=0,ylim=ylim,las=1,names.arg=names.arg, cex.names=0.5, xaxs="i",yaxs="i",border=cohort.color,
                  col=cohort.color,axes=F,ylab="",xlab="")
      
      #do some fancy footwork to get top line of histogram drawn but not lines in between bins
      bin.width <- x[2]-x[1]
      bw2 <- bin.width/2
      for (i in 1:nlengths) {
         xvals[2*i-1] <- x[i]-bw2            
         xvals[2*i] <- x[i]+bw2
         yvals[2*i-1] <- obs.data[yr,i]
         yvals[2*i] <- obs.data[yr,i]
      }
      lines(xvals,yvals,col=subtle.color,lwd=0.5)
      
      #now plot the axes
      if (yr %% mfcol[1] == 0 || yr==nyears) {
         axis(side=1,at=x,lab=length.list, line=-0.1,col.axis=subtle.color, col=subtle.color,lwd=0,lwd.ticks=0)  #just use for the labels, to allow more control than names.arg
      }
      if (yr <= mfcol[1]) {
        axis(2,las=1,at=c(0,0.15),col=subtle.color,col.axis=subtle.color,lwd=0.5)
      }
      par(new=T)
      par(xpd=NA)
      #print(pred.data[yr,])
      plot(x=x,y=pred.data[yr,],ylim=ylim, xlim=par("usr")[1:2], las=1,xaxs="i",yaxs="i",col="black",bg="white",pch=19,cex=0.6,axes=F,ylab="",xlab="")
      box(col=subtle.color,lwd=0.5)
      x.pos <- par("usr")[1] + 0.15*diff(par("usr")[1:2])   #par("usr") spits out the current coordinates of the plot window
      y.pos <- par("usr")[3] + 0.82*diff(par("usr")[3:4])   #par("usr") spits out the current coordinates of the plot window
      text(x=x.pos,y=y.pos,years[yr],cex=1.2, col=subtle.color)
      par(xpd=T)
   }
   mtext(side=1,outer=T,"Lengths",line=2.3)
   mtext(side=2,outer=T,"Proportions",line=3.2)
   mtext(side=3,outer=T,line=1.2,"Longline 1 length frequency")
   mtext(side=3,outer=T,line=0.2,paste("(",case_label,")",sep=""),cex=0.6)
}
pdf("lf1.pdf",width=8, height=11.5)
LF_Fit()
dev.off()
 # win.graph(width=7,height=11.5)
# LengthsLL1(labrep.file="sbtmod22_lab.rep",case_label="c1s1l1orig.5_h1m1M1O1C2a1")

#Required the output file _lab.rep.
#Requires library PBSmodelling, available from 
#Programmed by Trevor A. Branch 29 June 2009
#outputs from a _lab.rep file, assuming that the naming convention in the file has sections of 
#outputs separated by $name comments so that it can be easily parsed by readList() function of Jon Schnute 
#IMPORTANT: library PBSmodelling required, available at http://code.google.com/p/pbs-software/downloads/list
#then in R go to menu Packages, Install Package from Local zip file
#Run by:
###source("AussieAgeFits.r")
##########################################################################
#AgeFits <- function(labrep.file="mod1_r.rep", case_label="Mystuff") {
   #library(PBSmodelling)
   #subtle.color <- "gray40"
   #x <- readList(labrep.file)
   #ages <- 1:length(x$pobs_fsh_1)  #age range
   #obs.data <- x$pobs_fsh_1
   #pred.data <- x$phat_fsh_1
   #years <- 1977:2008
   
   #nyears <- length(years)
   #ages.list <- ages[1]:ages[2]
   #nages <- length(ages.list)
#
   #mfcol <- c(ceiling(nyears/3),3)
   #par(mfcol=mfcol,oma=c(3.5,4.5,3.5,1),mar=c(0,0,0,0))
   #cohort.color <- rainbow(mfcol[1]+2)[-c(1:2)]   #use hideous rainbow colors because they loop more gracefully than rich.colors
   #ncolors <- length(cohort.color)
   
   #ylim <- c(0,1.05*max(obs.data,pred.data))
   #for (yr in 1:nyears) { 
      #names.arg <- rep("",nages)
      #x <- barplot(obs.data[yr,],space=0.2,ylim=ylim,las=1,names.arg=names.arg, cex.names=0.5, xaxs="i",yaxs="i",border=subtle.color,
                  #col=cohort.color[1:nages],axes=F,ylab="",xlab="")
      #cohort.color <- c(cohort.color[ncolors],cohort.color[-1*ncolors])  #loop around colors
      #if (yr %% mfcol[1] == 0) {
         #axis(side=1,at=x,lab=ages.list, line=-0.1,col.axis=subtle.color, col=subtle.color,lwd=0,lwd.ticks=0)  #just use for the labels, to allow more control than names.arg
      #}
      #if (yr <= mfcol[1]) {
        #axis(2,las=1,at=c(0,0.5),col=subtle.color,col.axis=subtle.color,lwd=0.5)
      #}
      #par(new=T)
      #par(xpd=NA)
      #plot(x=x,y=pred.data[yr,],ylim=ylim, xlim=par("usr")[1:2], las=1,xaxs="i",yaxs="i",col="black",bg="white",pch=19,cex=1.3,axes=F,ylab="",xlab="")
      #box(col=subtle.color,lwd=0.5)
      #x.pos <- par("usr")[1] + 0.85*diff(par("usr")[1:2])   #par("usr") spits out the current coordinates of the plot window
      #y.pos <- par("usr")[3] + 0.75*diff(par("usr")[3:4])   #par("usr") spits out the current coordinates of the plot window
      #text(x=x.pos,y=y.pos,years[yr],cex=1.2, col=subtle.color)
      #par(xpd=T)
   #}
   #mtext(side=1,outer=T,"Age",line=2)
   #mtext(side=2,outer=T,"Proportion",line=3.2)
   #mtext(side=3,outer=T,line=1.2,"Australian age composition data")
   #mtext(side=3,outer=T,line=0.2,paste("(",case_label,")",sep=""),cex=0.6)
#}
#pdf("figs\\AussieAge.pdf",width=6, height=11.5)
#win.graph(width=8,height=11.5)
#AussieAgeFits(labrep.file="sbtmod22_lab.rep",case_label="c1s1l1orig.5_h1m1M1O1C2a1")
#dev.off()

