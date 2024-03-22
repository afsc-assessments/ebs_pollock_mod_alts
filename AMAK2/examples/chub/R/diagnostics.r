#-------------------------------------------------------------------------------
# Generic diagnostics code, developed for the joint-jack mackerel stock assessment
# By: Niels Hintzen
# Updated: 31 aug 2011
#-------------------------------------------------------------------------------

diagnostics <- function(jjm.out,jjm.in,jjm.ypr,what){

  #- Generic attributes of the stock assessment
  Nfleets   <- length(c(jjm.out$Fshry_names))
  Nsurveys  <- length(c(jjm.out$Index_names))
  ages      <- jjm.in$ages[1]:jjm.in$ages[2]
  lengths   <- jjm.in$lengths[1]:jjm.in$lengths[2]
  if(length(grep("pobs_fsh_",names(jjm.out)))>0){
    ageFleets <- unlist(strsplit(names(jjm.out)[grep("pobs_fsh_",names(jjm.out))],"_"))
    ageFleets <- ageFleets[seq(3,length(ageFleets),3)]
  } else { ageFleets <- 0}
  if(length(grep("pobs_len_fsh_",names(jjm.out)))>0){
    lgtFleets <- unlist(strsplit(names(jjm.out)[grep("pobs_len_fsh_",names(jjm.out))],"_"))
    lgtFleets <- lgtFleets[seq(4,length(lgtFleets),4)]
  } else { lgtFleets <- 0}
  
  #- Generic functions to be called in the diagnostics plot
  createDataFrame <- function(data,years,class){
                        dims  <- dim(data)
                        res   <- data.frame(year=rep(years,dims[2]),data=c(data),class=rep(class,each=dims[1]))
                     return(res)}

  plot.bubbles    <- function(x, xlab=" ", ylab=" ", main=" ", factx=0, facty=0, amplify=1){
                        my.col <- c("white"," ","black")
                        xval <- c(col(x)) + factx
                        yval <- c(row(x)) + facty
                        plot(x=xval, y=yval, cex=amplify*c(sqrt(abs(x/pi))), xlab=xlab, ylab=ylab, main=main, pch=19,xaxt="n",yaxt="n", col = my.col[2 + check.zero(c(x)/check.zero(abs(c(x))))]) ## area of bubble is proportional to value.
                        points(x=xval, y=yval, cex=amplify*c(sqrt(abs(x/pi)))) ## area of bubble is proportional to value.
                     }

  check.zero      <- function(x){
                        ## checks if there are zeros and replaces them with 1.
                        x[x == 0] <- 1
                        return(x)
                     }

  bubbles         <- function(x, data, bub.scale=2.5, col=c("black","black"),...){
                            	dots <- list(...)
                            	dots$data <- data
                            	dots$cex <- bub.scale*(abs(data$resids)/max(abs(data$resids),na.rm=T))+bub.scale*0.1
                            	dots$col <- ifelse(data$resids>0, col[1], col[2])
                            	dots$pch <- ifelse(data$resids>0,19,1)
                            	dots$panel <- function(x, y, ..., cex, subscripts){
                                panel.grid(h=-1, v= -1)
                                panel.xyplot(x, y, cex=cex[subscripts], ...)
                            	}
                            	call.list <- c(x=x, dots)
                            	ans <- do.call("xyplot", call.list)
                            	ans
                    }

                     
  ac              <- function(x){return(as.character(x))}
  an              <- function(x){return(as.numeric(x))}
  
  #-----------------------------------------------------------------------------
  #- Plots of the input data
  #-----------------------------------------------------------------------------
if("input" %in% what){
  pic <- plot(1,1,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",main="Input data",cex.main=3)
  print(pic)
  
  # 1: Weight in the fishery by fleet

  for(iFleet in 1:Nfleets){
    res <- createDataFrame(get("jjm.out")[[paste("wt_fsh_",iFleet,sep="")]][,-1],jjm.out$Yr,ages)
    if(iFleet == 1) tot <- cbind(res,c(jjm.out$Fshry_names[iFleet]))
    if(iFleet != 1) tot <- rbind(tot,cbind(res,c(jjm.out$Fshry_names[iFleet])))
  }
  colnames(tot) <- c("year","data","age","fleet"); res <- tot
  cols <-  brewer.pal(11, "Paired")
  ikey  <- simpleKey(text=as.character(ages),
                             points=F,lines=T,columns = 4,title="Ages")
  ikey$lines$col <- cols
  ikey$lines$lwd <- 4


  pic <- xyplot(data~year|fleet,data=res,groups=age,
         type="l",lwd=2,
         xlab="Years",ylab="Weight",main="Weight at age in the fishery",col=cols,
         key=ikey)

  print(pic)

  # 2: Weight at age in the survey
  for(iSurvey in 1:Nsurveys){
    res <- createDataFrame(get("jjm.out")[[paste("wt_ind_",iSurvey,sep="")]][,-1],jjm.out$Yr,ages)
    if(iSurvey == 1) tot <- cbind(res,c(jjm.out$Index_names[iSurvey]))
    if(iSurvey != 1) tot <- rbind(tot,cbind(res,c(jjm.out$Index_names[iSurvey])))
  }
  colnames(tot) <- c("year","data","age","survey"); res <- tot

  ikey  <- simpleKey(text=as.character(ages),
                             points=F,lines=T,columns = 4,title="Ages")
  ikey$lines$col <- cols
  ikey$lines$lwd <- 4


  pic <- xyplot(data~year|survey,data=res,groups=age,
         type="l",lwd=2,
         xlab="Years",ylab="Weight",main="Weight at age in the survey",col=cols,
         key=ikey)

  print(pic)

  # 3: Weight by cohort in the fleet
  for(iFleet in 1:Nfleets){
    res <- createDataFrame(get("jjm.out")[[paste("wt_fsh_",iFleet,sep="")]][,-1],jjm.out$Yr,ages)
    if(iFleet == 1) tot <- cbind(res,c(jjm.out$Fshry_names[iFleet]))
    if(iFleet != 1) tot <- rbind(tot,cbind(res,c(jjm.out$Fshry_names[iFleet])))
  }
  colnames(tot) <- c("year","data","age","fleet"); res <- tot
  res$cohort <- res$year - res$age
  cols <-  brewer.pal(12, "Paired")
  yrs <- sort(rev(jjm.out$Yr)[1:13])
  ikey  <- simpleKey(text=as.character(unique(subset(res,cohort%in%yrs)$cohort)),
                             points=F,lines=T,columns = 4,title="Cohort")
  ikey$lines$col <- rev(cols[1:length(unique(subset(res,cohort%in%yrs)$cohort))])
  ikey$lines$lwd <- 4

  pic <- xyplot(data~age|fleet,data=subset(res,cohort%in%yrs),groups=cohort,
         type="b",lwd=2,pch=19,cex=0.6,
         xlab="Age",ylab="Weight",main="Weight at age by cohort in the fleet",col=rev(cols[1:length(unique(subset(res,cohort%in%yrs)$cohort))]),
         key=ikey,
         panel=function(...){
          panel.grid(h=-1, v= -1)
          panel.xyplot(...)
         })

  print(pic)
  # 4: Weight by cohort in the survey

  for(iSurvey in 1:Nsurveys){
    res <- createDataFrame(get("jjm.out")[[paste("wt_ind_",iSurvey,sep="")]][,-1],jjm.out$Yr,ages)
    if(iSurvey == 1) tot <- cbind(res,c(jjm.out$Index_names[iSurvey]))
    if(iSurvey != 1) tot <- rbind(tot,cbind(res,c(jjm.out$Index_names[iSurvey])))
  }
  colnames(tot) <- c("year","data","age","survey"); res <- tot
  res$cohort <- res$year - res$age

  yrs <- sort(rev(jjm.out$Yr)[1:13])
  ikey  <- simpleKey(text=as.character(unique(subset(res,cohort%in%yrs)$cohort)),
                             points=F,lines=T,columns = 4,title="Cohort")
  ikey$lines$col <- rev(cols[1:length(unique(subset(res,cohort%in%yrs)$cohort))])
  ikey$lines$lwd <- 4


  pic <- xyplot(data~age|survey,data=subset(res,cohort%in%yrs),groups=cohort,
         type="b",lwd=2,pch=19,cex=0.6,
         xlab="Age",ylab="Weight",main="Weight at age by cohort in the survey",col=rev(cols[1:length(unique(subset(res,cohort%in%yrs)$cohort))]),
         key=ikey,
         panel=function(...){
          panel.grid(h=-1, v= -1)
          panel.xyplot(...)
         })

  print(pic)
  
  # 5: Age composition of the catch
  if(an(ageFleets)[1] != 0){
    for(iFleet in an(ageFleets)){
      res <- createDataFrame(sweep(jjm.in$Fagecomp[,,iFleet],1,apply(jjm.in$Fagecomp[,,iFleet],1,sum,na.rm=T),"/"),
                             jjm.in$years[1]:jjm.in$years[2],ages)
      res <- cbind(res,jjm.in$Fnames[iFleet])
      if(iFleet == an(ageFleets)[1]) tot <- res
      if(iFleet != an(ageFleets)[1]) tot <- rbind(tot,res)
    }
    colnames(tot) <- c("year","data","age","fleet"); res <- tot

    cols  <- brewer.pal(11, "PiYG")
    pic<- wireframe(data ~ age*year | fleet,data=res,
              scales = list(arrows = FALSE),
              col.regions=cols,cuts=10,
              main="Age composition in fleets",
              zlab="",ylab="",colorkey=T,drape=T,screen=list(z=-60,x=-60))
    print(pic)

    pic <- levelplot(data ~ age*year | fleet,data=res,
              scales = list(arrows = FALSE),
              col.regions=cols,cuts=10,
              main="Age composition in fleets",
              zlab="",ylab="",colorkey=T)
    print(pic)

    res$cohort <- (res$year - res$age) %% length(ages) + 1
    for(iFleet in unique(res$fleet)){
      pic <-   xyplot(data ~ age | as.factor(year),data=subset(res,fleet==iFleet),
                type="h",col=subset(res,fleet==iFleet)$cohort,lwd=5,
                main=paste("Age composition in fleet",iFleet))
      print(pic)
    }
  }
  
  # 6: Age composition of the survey
  if(length(which(jjm.in$Inumageyears>0))>0){
    for(iSurvey in which(jjm.in$Inumageyears>0)){
      res <- createDataFrame(sweep(jjm.in$Ipropage[,,iSurvey],1,apply(jjm.in$Ipropage[,,iSurvey],1,sum,na.rm=T),"/"),
                             jjm.in$years[1]:jjm.in$years[2],ages)
      res <- cbind(res,jjm.in$Inames[iSurvey])
      if(iSurvey == 1) tot <- res
      if(iSurvey != 1) tot <- rbind(tot,res)
    }
    colnames(tot) <- c("year","data","age","survey"); res <- tot

    pic <- levelplot(data ~ age*year | survey,data=res,
              scales = list(arrows = FALSE),
              col.regions=cols,cuts=10,
              main="Age composition in surveys",
              zlab="",ylab="",colorkey=T)
    print(pic)
  }

  # 7: Weight in the population
  res <- data.frame(year=1,data=jjm.in$Pwtatage,age=ages)
  pic <- xyplot(weight~age,data=res,
          xlab="Age",ylab="Weight",main="Weight in the stock",
          panel=function(...){
          panel.grid(h=-1, v= -1)
          panel.xyplot(...,type="b",pch=19,cex=0.6,lwd=2,col=1)
        })
  print(pic)
  
  # 8: Maturity at age in the population
  res <- data.frame(year=1,data=jjm.in$Pmatatage,age=ages)
  pic <- xyplot(maturity~age,data=res,
          xlab="Age",ylab="Maturity",main="Maturity in the stock",
          panel=function(...){
          panel.grid(h=-1, v= -1)
          panel.xyplot(...,type="b",pch=19,cex=0.6,lwd=2,col=1)
        })
  print(pic)
  
  # 9: Length composition of the catch
  if(an(lgtFleets)[1] != 0){
     for(iFleet in an(lgtFleets)){
      res <- createDataFrame(sweep(jjm.in$Flengthcomp[,,iFleet],1,apply(jjm.in$Flengthcomp[,,iFleet],1,sum,na.rm=T),"/"),
                             jjm.in$years[1]:jjm.in$years[2],lengths)
      res <- cbind(res,jjm.in$Fnames[iFleet])
      if(iFleet == an(lgtFleets)[1]) tot <- res
      if(iFleet != an(lgtFleets)[1]) tot <- rbind(tot,res)
    }
    colnames(tot) <- c("year","data","length","fleet"); res <- tot


    cols  <- brewer.pal(11, "PiYG")
    for(iFleet in unique(tot$fleet)){
      pic <- levelplot(data ~ length*year | fleet,data=subset(tot,fleet==iFleet),
                scales = list(arrows = FALSE),
                col.regions=cols,cuts=10,
                main=paste("Length composition in surveys",iFleet),
                zlab="",ylab="",colorkey=T)
      print(pic)
    }

    tot$cohort <- (tot$year - tot$length) %% length(lengths) + 1
    for(iFleet in unique(tot$fleet)){
      pic <-   xyplot(data ~ length | as.factor(year),data=subset(tot,fleet==iFleet),
                type="h",col=subset(tot,fleet==iFleet)$cohort,lwd=4,
                main=paste("Length composition in survey",iFleet),
                scales=list(alternating=1,y=list(relation="free",rot=0)))
      print(pic)
    }
    
  }
}

  #-----------------------------------------------------------------------------
  #- Plots of the fit of the catch data
  #-----------------------------------------------------------------------------
if("fit" %in% what){
  pic <- plot(1,1,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",main="Fit of catch data",cex.main=3)
  print(pic)

  # 9: Trend in catch
  for(iFleet in 1:Nfleets){
    res <- cbind(jjm.out$Yr,jjm.out[[paste("Obs_catch_",iFleet,sep="")]])
    if(Nfleets == 1) res <- cbind(res,jjm.out$Fshry_names[iFleet])
    if(Nfleets > 1)  res <- cbind(res,jjm.out$Fshry_names[iFleet,1])
    if(iFleet == 1) tot <- res
    if(iFleet != 1) tot <- rbind(tot,res)
  }
  colnames(tot) <- c("year","catch","fleet"); res <- data.frame(tot)
  res$catch <- as.numeric(as.character(res$catch))

  totcatch <- numeric()
  for(iYr in sort(unique(res$year))){
    totcatch[which(iYr==sort(unique(res$year)))] <- sum(subset(res,year==iYr)$catch)
  }

  pic <- xyplot(totcatch~jjm.out$Yr,
        xlab="Years",ylab="Catch in kt",main="Total catch",
        panel=function(...){
          panel.grid(h=-1, v= -1)
          panel.barchart(...,horizontal=F,origin=0,box.width=1,col="grey")
        })
  print(pic)

  #9b: trends in catch by fleet as polygon
  res$year <- an(ac(res$year))
  if(Nfleets > 1){
    for(iFleet in 2:Nfleets){
      idx <- which(res$fleet == jjm.out$Fshry_names[iFleet])
      res[idx,"catch"] <- subset(res,fleet == jjm.out$Fshry_names[iFleet])$catch + subset(res,fleet == jjm.out$Fshry_names[iFleet-1])$catch
    }

    ikey                  <- simpleKey(text=jjm.out$Fshry_names,
                                       points=F,lines=T,columns = 2)

    ikey$lines$col       <- 1:4
    ikey$lines$lwd       <- 2
    ikey$lines$lty       <- 1

    pic <- xyplot(catch ~ year,data=res,groups=fleet,
           xlab="Years",ylab="Catch by fleet in kt",main="Total catch by fleet",
           type="l",key=ikey,
           panel=function(...){
            lst <- list(...)
            idx     <- mapply(seq,from=seq(1,length(lst$y),length(lst$y)/4),to=seq(1,length(lst$y),length(lst$y)/4)+(length(lst$y)/4-1))
            panel.grid(h=-1,v=-1)
            panel.xyplot(...,col="white")
            for(iFleet in Nfleets:1){
              panel.polygon(x=c(lst$x[idx[,iFleet]],rev(lst$x[idx[,iFleet]])),
                            y=c(rep(0,length(lst$y[idx[,iFleet]])),rev(lst$y[idx[,iFleet]])),
                            col=iFleet,border=0)
            }

           })
    print(pic)
  }

  # 10: Log residual total catch by fleet
  for(iFleet in 1:Nfleets){
    res <- data.frame(year=jjm.out$Yr,obs=jjm.out[[paste("Obs_catch_",iFleet,sep="")]],model=jjm.out[[paste("Pred_catch_",iFleet,sep="")]],fleet=jjm.out$Fshry_names[iFleet])
    if(iFleet == 1) tot <- res
    if(iFleet != 1) tot <- rbind(tot,res)
  }
  tot$obs[tot$obs<=0]   <- NA;
  tot$obs[tot$model<=0] <- NA
  res                   <- tot

  resids                <- log(res$obs + 1) - log(res$model+1); res$resids <- resids
    scalar              <- 3/max(resids,na.rm=T)
  residRange            <- range(resids,na.rm=T)
  ikey                  <- simpleKey(text=as.character(round(seq(residRange[1],residRange[2],length.out=6),2)),
                                     points=T,lines=F,columns = 2)
  ikey$points$cex       <- abs(round(seq(residRange[1],residRange[2],length.out=6),2))*scalar
  ikey$points$col       <- 1
  ikey$points$pch       <- ifelse(round(seq(residRange[1],residRange[2],length.out=6),2)>0,19,1)


  pic <- xyplot(resids*scalar ~ year | as.factor(fleet),data=res,
         xlab="Years",ylab="Residuals",main="Catch residuals by fleet",
         prepanel=function(...) {list(ylim=c(1,1))},
         layout=c(1,Nfleets),
         type="p",lwd=3,
         cex.axis=1.2,font=2,
         key=ikey,
         scales=list(y=list(draw=F)),
         panel=function(x,y){
          panel.grid(h=-1, v= -1)
          panel.points(x,1,cex=abs(y),col=ifelse(y>0,"black","white"),pch=19)
          panel.points(x,1,cex=abs(y),col=1,pch=1)
         })
  print(pic)
  
  # 11: Absolute residual catch by fleet
  pic <- xyplot(model - obs ~ year | fleet,data=res,allow.multiple=T,
          xlab="Years",ylab="Absolute residual catch",main="Absolute residual catch by fleet",
         panel=function(...){
          panel.grid(h=-1,v=-1,lty=3)
          panel.barchart(...,horizontal=FALSE,origin=0,box.width=1,col="grey")
         })

  print(pic)

  # 12a: Proportions catch by age modelled and observed
  if(an(ageFleets)[1] != 0){
    for(iFleet in an(ageFleets)){
      obs <- createDataFrame(jjm.out[[paste("pobs_fsh_",iFleet,sep="")]][,-1],jjm.out[[paste("pobs_fsh_",iFleet,sep="")]][,1],ages)
      mod <- createDataFrame(jjm.out[[paste("phat_fsh_",iFleet,sep="")]][,-1],jjm.out[[paste("phat_fsh_",iFleet,sep="")]][,1],ages)
      if(iFleet == an(ageFleets)[1]) tot <- cbind(obs,mod$data,rep(jjm.out$Fshry_names[iFleet,1],nrow(obs)))
      if(iFleet != an(ageFleets)[1]) tot <- rbind(tot,cbind(obs,mod$data,rep(jjm.out$Fshry_names[iFleet,1],nrow(obs))))
    }
    colnames(tot) <- c("year","obs","age","model","fleet");
    res           <- tot

    resids        <- log(res$obs + 1) - log(res$model + 1)
    res$resids    <- resids
      scalar              <- 3/max(resids,na.rm=T)
    residRange            <- range(resids,na.rm=T)
    ikey                  <- simpleKey(text=as.character(round(seq(residRange[1],residRange[2],length.out=6),2)),
                                       points=T,lines=F,columns = 2)
    ikey$points$cex       <- abs(round(seq(residRange[1],residRange[2],length.out=6),2))*scalar
    ikey$points$col       <- 1
    ikey$points$pch       <- ifelse(round(seq(residRange[1],residRange[2],length.out=6),2)>0,19,1)

    pic <- bubbles(age~year|fleet,data=res,allow.multiple=T,
                   xlab="Years",ylab="Age",main="Residuals catch-at-age by fleet",key=ikey)
    print(pic)
  }

  # 12b: Proportions catch by length modelled and observed
  if(an(lgtFleets)[1] != 0){
    for(iFleet in an(lgtFleets)){
      obs <- createDataFrame(jjm.out[[paste("pobs_len_fsh_",iFleet,sep="")]][,-1],jjm.out[[paste("pobs_len_fsh_",iFleet,sep="")]][,1],lengths)
      mod <- createDataFrame(jjm.out[[paste("phat_len_fsh_",iFleet,sep="")]][,-1],jjm.out[[paste("phat_len_fsh_",iFleet,sep="")]][,1],lengths)
      if(Nfleets == 1){
        if(iFleet == an(lgtFleets)[1]) tot <- cbind(obs,mod$data,rep(jjm.out$Fshry_names[iFleet],nrow(obs)))
        if(iFleet != an(lgtFleets)[1]) tot <- rbind(tot,cbind(obs,mod$data,rep(jjm.out$Fshry_names[iFleet],nrow(obs))))
      }
      if(Nfleets != 1){
        if(iFleet == an(lgtFleets)[1]) tot <- cbind(obs,mod$data,rep(jjm.out$Fshry_names[iFleet,1],nrow(obs)))
        if(iFleet != an(lgtFleets)[1]) tot <- rbind(tot,cbind(obs,mod$data,rep(jjm.out$Fshry_names[iFleet,1],nrow(obs))))
      }
    }
    colnames(tot) <- c("year","obs","length","model","fleet");
    res           <- tot

    resids        <- log(res$obs + 1) - log(res$model + 1)
    res$resids    <- resids
      scalar              <- 3/max(resids,na.rm=T)
    residRange            <- range(resids,na.rm=T)
    ikey                  <- simpleKey(text=as.character(round(seq(residRange[1],residRange[2],length.out=6),2)),
                                       points=T,lines=F,columns = 2)
    ikey$points$cex       <- abs(round(seq(residRange[1],residRange[2],length.out=6),2))*scalar
    ikey$points$col       <- 1
    ikey$points$pch       <- ifelse(round(seq(residRange[1],residRange[2],length.out=6),2)>0,19,1)

    pic <- bubbles(length~year|fleet,data=res,allow.multiple=T,
                   xlab="Years",ylab="Length",main="Residuals catch-at-length by fleet",key=ikey)
    print(pic)
  }


  # 13a: Fitted age by year by fleet
  if(an(ageFleets)[1] != 0){
    for(iFleet in an(ageFleets)){
      obs <- createDataFrame(jjm.out[[paste("pobs_fsh_",iFleet,sep="")]][,-1],jjm.out[[paste("pobs_fsh_",iFleet,sep="")]][,1],ages)
      mod <- createDataFrame(jjm.out[[paste("phat_fsh_",iFleet,sep="")]][,-1],jjm.out[[paste("phat_fsh_",iFleet,sep="")]][,1],ages)
      if(iFleet == an(ageFleets)[1]){
        x <- cbind(obs,rep("obs",nrow(obs)),jjm.out$Fshry_names[iFleet]); colnames(x) <- c("year","data","age","class","fleet")
        y <- cbind(mod,rep("model",nrow(mod)),jjm.out$Fshry_names[iFleet]);colnames(y) <- c("year","data","age","class","fleet")
        tot <- rbind(x,y)
      }
      if(iFleet != an(ageFleets)[1]){
        x <- cbind(obs,rep("obs",nrow(obs)),jjm.out$Fshry_names[iFleet]); colnames(x) <- c("year","data","age","class","fleet")
        y <- cbind(mod,rep("model",nrow(mod)),jjm.out$Fshry_names[iFleet]);colnames(y) <- c("year","data","age","class","fleet")
        res <- rbind(x,y)
        tot <- rbind(tot,res)
      }
    }
    res <- tot
    res$cohort <- (res$year - res$age) %% length(ages) + 1

    ikey          <- simpleKey(text=c("Observed","Predicted"),
                               points=T,lines=F,rect=T,columns = 2)
    ikey$rectangles$alpha=c(1,0)
    ikey$rectangles$col="white"
    ikey$rectangles$lty=c(1,0)
    ikey$points$pch=c(-1,19)
    ikey$points$col=c("white","black")
    ikey$points$cex=c(0,1.1)

    cols  <- rainbow(length(ages))
    for(iFleet in c(jjm.out$Fshry_names)[an(ageFleets)]){
      tmpres  <- subset(res,fleet==iFleet)
    pic <- xyplot(data ~ age | factor(year),data=tmpres,
            groups=class,
            xlab="Age",ylab="Proportion at age",main=paste("Age fits",iFleet),
            key=ikey,
            as.table=TRUE,
            panel=function(x,y){
              idx     <- mapply(seq,from=seq(1,length(y),length(ages)),to=seq(1,length(y),length(ages))+(length(ages)-1))
              first   <- c(idx[,seq(1,dim(idx)[2],3)])
              second  <- c(idx[,seq(2,dim(idx)[2],3)])
              #cols <- tmpres$cohort[which(is.na(pmatch(tmpres$data,y))==F & is.na(pmatch(tmpres$age,x))==F)]
              yr      <- names(which.max(table(tmpres$year[which(tmpres$data %in% y)])))
              colidx  <- tmpres$cohort[which(tmpres$data %in% y & tmpres$year == an(yr))]
              panel.barchart(x[first],y[first],horizontal=F,origin=0,box.width=1,col=cols[colidx])
              panel.points(x[second],y[second],pch=19,col=1,cex=0.5)
              })
    print(pic)
    }
  }

  # 13b: Fitted length by year by fleet
  if(an(lgtFleets)[1] != 0){
    for(iFleet in an(lgtFleets)){
      obs <- createDataFrame(jjm.out[[paste("pobs_len_fsh_",iFleet,sep="")]][,-1],jjm.out[[paste("pobs_len_fsh_",iFleet,sep="")]][,1],lengths)
      mod <- createDataFrame(jjm.out[[paste("phat_len_fsh_",iFleet,sep="")]][,-1],jjm.out[[paste("phat_len_fsh_",iFleet,sep="")]][,1],lengths)
      if(iFleet == an(lgtFleets)[1]){
        x <- cbind(obs,rep("obs",nrow(obs)),jjm.out$Fshry_names[iFleet]); colnames(x) <- c("year","data","length","class","fleet")
        y <- cbind(mod,rep("model",nrow(mod)),jjm.out$Fshry_names[iFleet]);colnames(y) <- c("year","data","length","class","fleet")
        tot <- rbind(x,y)
      }
      if(iFleet != an(lgtFleets)[1]){
        x <- cbind(obs,rep("obs",nrow(obs)),jjm.out$Fshry_names[iFleet]); colnames(x) <- c("year","data","length","class","fleet")
        y <- cbind(mod,rep("model",nrow(mod)),jjm.out$Fshry_names[iFleet]);colnames(y) <- c("year","data","length","class","fleet")
        res <- rbind(x,y)
        tot <- rbind(tot,res)
      }
    }
    res <- tot
    res$cohort <- res$length

    ikey          <- simpleKey(text=c("Observed","Predicted"),
                               points=T,lines=F,rect=T,columns = 2)
    ikey$rectangles$alpha=c(1,0)
    ikey$rectangles$col="white"
    ikey$rectangles$lty=c(1,0)
    ikey$points$pch=c(-1,19)
    ikey$points$col=c("white","black")
    ikey$points$cex=c(0,1.1)

    cols  <- rainbow(length(lengths))
    for(iFleet in c(jjm.out$Fshry_names)[an(lgtFleets)]){
      tmpres  <- subset(res,fleet==iFleet)
    pic <- xyplot(data ~ length | factor(year),data=tmpres,
            groups=class,
            xlab="Length",ylab="Proportion at length",main=paste("Length fits",iFleet),
            key=ikey,
            as.table=TRUE,
            panel=function(x,y){
              idx     <- mapply(seq,from=seq(1,length(y),length(lengths)),to=(seq(1,length(y),length(lengths))+length(lengths)-1))
              first   <- c(idx[,seq(1,dim(idx)[2],3)])
              second  <- c(idx[,seq(2,dim(idx)[2],3)])
              #cols <- tmpres$cohort[which(is.na(pmatch(tmpres$data,y))==F & is.na(pmatch(tmpres$age,x))==F)]
              yr      <- names(which.max(table(tmpres$year[which(tmpres$data %in% y)])))
              colidx  <- tmpres$cohort[which(tmpres$data %in% y & tmpres$year == an(yr))]
              panel.barchart(x[first],y[first],horizontal=F,origin=0,box.width=1,col=cols[colidx])
              panel.points(x[second],y[second],col=1,pch=19,cex=0.25)
              },
              scales=list(alternating=1,y=list(relation="free",rot=0)))
    print(pic)
    }
  }


  
  cols  <- brewer.pal(11, "Paired")
  # 14: Absolute catch by fleet modelled and observed
  for(iFleet in 1:Nfleets){
    res <- data.frame(year=jjm.out$Yr,obs=jjm.out[[paste("Obs_catch_",iFleet,sep="")]],model=jjm.out[[paste("Pred_catch_",iFleet,sep="")]],fleet=jjm.out$Fshry_names[iFleet])
    if(iFleet == 1) tot <- rbind(cbind(res$year,res$obs,ac(res$fleet),rep("obs",nrow(res))),cbind(res$year,res$model,ac(res$fleet),rep("model",nrow(res))))
    if(iFleet != 1) tot <- rbind(tot,rbind(cbind(res$year,res$obs,ac(res$fleet),rep("obs",nrow(res))),cbind(res$year,res$model,ac(res$fleet),rep("model",nrow(res)))))
  }
  colnames(tot) <- c("year","data","fleet","classing")
  res <- data.frame(tot,stringsAsFactors=F); res$year <- an(ac(res$year)); res$data <- an(ac(res$data))

  ikey          <- simpleKey(text=c("Observed","Predicted"),
                             points=F,lines=T,rect=T,columns = 2)
  ikey$rectangles$alpha=c(1,0)
  ikey$rectangles$col="grey"
  ikey$rectangles$lty=c(1,0)
  ikey$lines$lty=c(0,1)
  ikey$lines$col=1
  ikey$lines$lwd=3

  pic <- xyplot(data ~ year | as.factor(fleet),data=res,
          groups=as.factor(classing),
          key=ikey,
          main="Predicted and observed catches by fleet",
          xlab="Years",ylab="Thousand tonnes",
          panel=function(x,y){
            first = 1:length(jjm.out$Yr)
            second = (length(jjm.out$Yr)+1):(length(jjm.out$Yr)*2)
            panel.grid(h=-1, v= -1)
            panel.barchart(x[first],y[second],horizontal=FALSE,origin=0,box.width=1,col="grey")
            panel.xyplot(x[first],y[second],type="l",lwd=5,col=1,lty=1)
          })

  print(pic)

  #-----------------------------------------------------------------------------
  #- Plots of the fit of the survey data
  #-----------------------------------------------------------------------------

  pic <- plot(1,1,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",main="Fit of survey data",cex.main=3)
  print(pic)


  # 15: Standardized indices observed with error and modelled
  rm(tot); rm(res)
  for(iSurvey in 1:Nsurveys){
    if(any(jjm.out$Yr %in% jjm.out[[paste("Obs_Survey_",iSurvey,sep="")]][,1]==T)){
      addToDF <- jjm.out$Yr[which(!jjm.out$Yr %in% jjm.out[[paste("Obs_Survey_",iSurvey,sep="")]][,1])]
      addToDF <- as.data.frame(rbind(cbind(addToDF,NA,"model"),cbind(addToDF,NA,"obs"),cbind(addToDF,NA,"sd")))
      colnames(addToDF) <- c("year","data","class"); addToDF$year <- an(ac(addToDF$year)); addToDF$data <- an(ac(addToDF$data)); addToDF$class <- as.character(addToDF$class)
    }
    res <- createDataFrame(jjm.out[[paste("Obs_Survey_",iSurvey,sep="")]][,-1],jjm.out[[paste("Obs_Survey_",iSurvey,sep="")]][,1],c("obs","model","sd"))
    res$class <- ac(res$class)
    res$data <- res$data/max(res$data,na.rm=T)
    resSort <- rbind(res,addToDF)
    resSort <- orderBy(~year+class,data=resSort)
    if(iSurvey == 1) tot <- cbind(resSort,rep(jjm.out$Index_names[iSurvey,1],nrow(resSort)))
    if(iSurvey != 1){
      res2  <- cbind(resSort,rep(jjm.out$Index_names[iSurvey,1],nrow(resSort)))
      tot   <- rbind(tot,res2)
    }
  }

  colnames(tot) <- c("year","data","classing","surveys")
  tot <- tot[duplicated(paste(tot$year,tot$classing,tot$surveys))==F,]
  tot$data[which(tot$data<0)] <- NA
  res <- tot

  ikey          <- simpleKey(text=c("Observed","Predicted"),
                             points=T,lines=T,columns = 2)

  ikey$lines$lty=c(0,1)
  ikey$lines$lwd=c(0,2)
  ikey$lines$col=c(0,1)
  ikey$points$pch=c(19,-1)
  ikey$points$col=c("grey","white")
  ikey$points$cex=0.9

  pic <- xyplot(data~year|as.factor(surveys),data=res,
        groups=classing,
        main="Predicted and observed indices",xlab="Years",ylab="Index value",
        key=ikey,ylim=c(0,2),
        panel=function(...){
          tmp   <- list(...)
          first <- which(tmp$groups[1:length(tmp$x)] == "model")
          second<- which(tmp$groups[1:length(tmp$x)] == "obs")
          third <- which(tmp$groups[1:length(tmp$x)] == "sd")
          panel.grid(h=-1, v= -1)
          panel.xyplot(tmp$x[first],tmp$y[first],type="l",col="black",lwd=3)
          panel.points(tmp$x[second],tmp$y[second],col="grey",pch=19,cex=0.6)
          panel.segments(tmp$x[third],c(tmp$y[second]+1.96*tmp$y[third]),tmp$x[third],c(tmp$y[second]-1.96*tmp$y[third]))
          panel.lines(tmp$x[first],tmp$y[first],col="black",lwd=3)
        })
  print(pic)

  # 16: Log residuals in survey
  for(iSurvey in 1:Nsurveys){
    if(any(jjm.out$Yr %in% jjm.out[[paste("Obs_Survey_",iSurvey,sep="")]][,1]==F)){
      addToDF <- jjm.out$Yr[which(!jjm.out$Yr %in% jjm.out[[paste("Obs_Survey_",iSurvey,sep="")]][,1])]
      addToDF <- data.frame(year=addToDF,obs=NA,model=NA,sd=NA)
      colnames(addToDF) <- c("year","obs","model","sd"); addToDF$year <- an(ac(addToDF$year))
    }

    res <- createDataFrame(jjm.out[[paste("Obs_Survey_",iSurvey,sep="")]][,-1],jjm.out[[paste("Obs_Survey_",iSurvey,sep="")]][,1],c("obs","model","sd"))
    res <- cbind(subset(res,class=="obs")[,1:2],subset(res,class=="model")$data,subset(res,class=="sd")$data)
    colnames(res) <- c("year","obs","model","sd")
    res <- rbind(res,addToDF)
    res <- orderBy(~year,data=res)
    if(iSurvey == 1) tot <- cbind(res,rep(jjm.out$Index_names[iSurvey,1],nrow(res)))
    if(iSurvey != 1){
      res2  <- cbind(res,rep(jjm.out$Index_names[iSurvey,1],nrow(res)))
      tot   <- rbind(tot,res2)
    }
  }
  colnames(tot) <- c("year","obs","model","sd","survey")
  tot           <- tot[duplicated(paste(tot$year,tot$survey)==F),]
  tot$obs[tot$obs<0] <- NA; tot$obs[tot$model<0] <- NA; tot$obs[tot$sd < 0] <- NA
  res <- tot

  resids        <- log(res$obs + 1) - log(res$model + 1)
  res$resids    <- resids
    scalar      <- 3/max(abs(resids),na.rm=T)
  resRange      <- range(resids,na.rm=T)
  ikey          <- simpleKey(text=as.character(round(seq(resRange[1],resRange[2],length.out=6),2)),
                             points=T,lines=F,columns = 2)
  ikey$points$cex <- abs(round(seq(resRange[1],resRange[2],length.out=6),2))*scalar
  ikey$points$col <- 1
  ikey$points$pch<- ifelse(round(seq(resRange[1],resRange[2],length.out=6),2)>0,19,1)


  pic <- xyplot(resids*scalar ~ year | as.factor(survey),data=res,
         xlab="Years",ylab="Residuals",main="Survey residuals",
         prepanel=function(...) {list(ylim=c(1,1))},
         layout=c(1,Nsurveys),
         type="p",col=cols,lwd=3,
         cex.axis=1.2,font=2,
         key=ikey,
         scales=list(y=list(draw=F)),
         panel=function(x,y){
          panel.grid(v= -1,h=1)
          panel.points(x,1,cex=abs(y),col=ifelse(y>0,"black","white"),pch=19)
          panel.points(x,1,cex=abs(y),col=1,pch=1)
         })

  print(pic)

  #-----------------------------------------------------------------------------
  #- Plots of selectivity in fleet and survey + F's
  #-----------------------------------------------------------------------------

  pic <- plot(1,1,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",main="Fleet & Survey sel & F",cex.main=3)
  print(pic)

  # 17: Selectivity at age in the fleet
  for(iFleet in 1:Nfleets){
    res <- createDataFrame(jjm.out[[paste("sel_fsh_",iFleet,sep="")]][,-c(1,2)],jjm.out$Yr,ages)
    res <- cbind(res,jjm.out$Fshry_names[iFleet])
    if(iFleet == 1) tot <- res
    if(iFleet != 1) tot <- rbind(tot,res)
  }
  colnames(tot) <- c("year","data","age","fleet"); res <- tot
  cols  <- brewer.pal(11, "PiYG")
  pic<- wireframe(data ~ age*year | fleet,data=res,
            scales = list(arrows = FALSE),
            col.regions=cols,cuts=10,
            main="Selectivity at age in fleets",
            zlab="",ylab="",colorkey=T,drape=T,screen=list(z=30,x=-60))
  print(pic)

  # 18: selecitivity at age in the survey
  for(iSurvey in 1:Nsurveys){
    res <- createDataFrame(jjm.out[[paste("sel_ind_",iSurvey,sep="")]][,-c(1,2)],jjm.out$Yr,ages)
    res <- cbind(res,jjm.out$Index_names[iSurvey])
    if(iSurvey == 1) tot <- res
    if(iSurvey != 1) tot <- rbind(tot,res)
  }
  colnames(tot) <- c("year","data","age","survey"); res <- tot

  cols  <- brewer.pal(11, "PiYG")
  pic <- wireframe(data ~ age*year | survey,data=res,
            scales = list(arrows = FALSE),
            col.regions=cols,cuts=10,
            main="Selectivity at age in surveys",
            zlab="",ylab="",colorkey=T,drape=T,screen=list(z=30,x=-60))
  print(pic)

  # 19: F at age

  res <- jjm.out$TotF
  dimnames(res) <- list(Years=jjm.out$Yr,Age=ages)
  pic <- levelplot(t(res),xlab="Age",ylab="Years",col.regions=rev(cols),cuts=10,main="F at age")
  print(pic)

  #20: Fisheries mean age
  if(an(ageFleets)[1] != 0){
    for(iFleet in an(ageFleets)){
      res <- data.frame(get("jjm.out")[[paste("EffN_Fsh_",iFleet,sep="")]][,c(1,4,5,7,8)])
      colnames(res) <- c("Year","Obs","Model","Obs5","Obs95")
      for(i in 2:5){
        tot <- data.frame(cbind(res[,1],res[,i]))
        tot$class <- names(res)[i]
        tot$Fleet <- jjm.out$Fshry_names[iFleet]
        if(iFleet == an(ageFleets)[1] & i == 2) total <- tot
        if(iFleet != an(ageFleets)[1] | i != 2) total <- rbind(total,tot)
      }
    }
    colnames(total) <- c("year","data","class","fleet")

    ikey           <- simpleKey(text=c("Observed","Modelled"),
                               points=T,lines=T,columns = 2)
    ikey$lines$col <- c("white","black")
    ikey$lines$lwd <- c(0,2)
    ikey$lines$lty <- c(0,1)
    ikey$lines$pch <- c(0,0)
    ikey$points$pch<- c(16,0)
    ikey$points$col<- c("grey","white")


    pic <- xyplot(data ~ year | fleet,data=total,
            type="l",lwd=3,lty=c(1,3),col=1,
            ylab="Age",xlab="Years",main="Fishery mean age",
            key=ikey,
            panel=function(x,y){
              panel.grid(h=-1, v= -1)
              idx <- mapply(seq(length(x)/4,length(x),length.out=4)-length(x)/4+1,
                            seq(length(x)/4,length(x),length.out=4),FUN=seq)
              obs   <- idx[,1]
              mod   <- idx[,2]
              obs5  <- idx[,3]
              obs95 <- idx[,4]

              panel.xyplot(x[obs],y[obs],type="p",col="grey",pch=19,cex=0.6)
              panel.segments(x[obs],y[obs5],x[obs],y[obs95])
              panel.xyplot(x[obs],y[mod],type="l",lwd=2,col="black")

            })
    print(pic)
  }

  #21: Survey mean age
  if(an(ageFleets)[1] != 0){
    for(iSurvey in 1:Nsurveys){
      res <- data.frame(get("jjm.out")[[paste("EffN_Survey_",iSurvey,sep="")]][,c(1,4,5,7,8)])
      if(nrow(res)>1){
        colnames(res) <- c("Year","Obs","Model","Obs5","Obs95")
        for(i in 2:5){
          tot <- data.frame(cbind(res[,1],res[,i]))
          tot$class <- names(res)[i]
          tot$Survey <- jjm.out$Index_names[iSurvey]
          if(iSurvey == 1 & i == 2) total <- tot
          if(iSurvey != 1 | i != 2) total <- rbind(total,tot)
        }
      }
    }
    colnames(total) <- c("year","data","class","survey")

    ikey           <- simpleKey(text=c("Observed","Modelled"),
                               points=T,lines=T,columns = 2)
    ikey$lines$col <- c("white","black")
    ikey$lines$lwd <- c(0,2)
    ikey$lines$lty <- c(0,1)
    ikey$lines$pch <- c(0,0)
    ikey$points$pch<- c(16,0)
    ikey$points$col<- c("grey","white")


    pic <- xyplot(data ~ year | survey,data=total,
            type="l",lwd=3,lty=c(1,3),col=1,
            ylab="Age",xlab="Years",main="Survey mean age",
            key=ikey,
            panel=function(x,y){
              panel.grid(h=-1, v= -1)
              idx <- mapply(seq(length(x)/4,length(x),length.out=4)-length(x)/4+1,
                            seq(length(x)/4,length(x),length.out=4),FUN=seq)
              obs   <- idx[,1]
              mod   <- idx[,2]
              obs5  <- idx[,3]
              obs95 <- idx[,4]

              panel.xyplot(x[obs],y[obs],type="p",col="grey",pch=19,cex=0.6)
              panel.segments(x[obs],y[obs5],x[obs],y[obs95])
              panel.xyplot(x[obs],y[mod],type="l",lwd=2,col="black")

            })
    print(pic)
  }

  #-----------------------------------------------------------------------------
  #- Plots of stock summary
  #-----------------------------------------------------------------------------

  pic <- plot(1,1,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",main="Stock summary",cex.main=3)
  print(pic)


  # 22: summary sheet with SSB, R, F and Biomass and Catch
  TotCatch    <- 0
  for(iFlt in grep("Obs_catch_",names(jjm.out)))
    TotCatch    <- jjm.out[[iFlt]]+TotCatch
  summaryData <- rbind(cbind(jjm.out$Yr,jjm.out$TotBiom[,-1],"Total biomass"),
                       cbind(jjm.out$SSB[which(jjm.out$SSB[,1] %in% jjm.out$Yr),1],jjm.out$SSB[which(jjm.out$SSB[,1] %in% jjm.out$Yr),-1],"Spawning Stock Biomass"),
                       cbind(jjm.out$Yr,jjm.out$R[,-1],"Recruitment"),
                       cbind(jjm.out$Yr,cbind(rowMeans(jjm.out$TotF[,-1]),rowMeans(jjm.out$TotF[,-1]),rowMeans(jjm.out$TotF[,-1]),rowMeans(jjm.out$TotF[,-1])),"Fishing mortality"),
                       cbind(jjm.out$Yr,TotCatch,TotCatch,TotCatch,TotCatch,"Catches"))

  summaryData <- rbind(cbind(summaryData[,c(1:2,6)],"point"),cbind(summaryData[,c(1,4,6)],"lower"),cbind(summaryData[,c(1,5,6)],"upper"))

  colnames(summaryData) <- c("year","data","class","estim"); summaryData <- data.frame(summaryData,stringsAsFactors=F)
  summaryData$year <- as.integer(summaryData$year); summaryData$data <- as.numeric(summaryData$data)

  pic <- xyplot(data ~ year | class,data=summaryData,
         groups=class,
         main="Summary sheet",
         prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
         layout=c(1,5),
         panel=function(x,y){
          panel.grid(h=-1, v= -1)
          point = 1                       :   length(jjm.out$Yr)
          lower = (  length(jjm.out$Yr)+1):(2*length(jjm.out$Yr))
          upper = (2*length(jjm.out$Yr)+1):(3*length(jjm.out$Yr))
          #catches
          if(panel.number()==1){
            panel.barchart(x[point],y[point],horizontal=FALSE,origin=0,box.width=1,col="grey")
          }
          #SSB
          if(panel.number()==4){
            panel.polygon(c(x[lower],rev(x[upper])),c(y[lower],rev(y[upper])),col="grey")
            panel.xyplot(x[point],y[point],type="l",lwd=6,lty=3,col=1)
          }
          #Recruitment
          if(panel.number()==3){
            panel.barchart(x[point],y[point],horizontal=FALSE,origin=0,box.width=1,col="grey")
            panel.segments(x[lower],y[lower],x[lower],y[upper])
          }
          #F
          if(panel.number()==2) panel.xyplot(x[point],y[point],lwd=4,lty=2,type="l",col=1)

          #Total biomass
          if(panel.number()==5){
            panel.polygon(c(x[lower],rev(x[upper])),c(y[lower],rev(y[upper])),col="grey")
            panel.xyplot(x[point],y[point],lwd=4,lty=1,type="l",col=1)
          }
        },
        scales=list(alternating=1,y=list(relation="free",rot=0)))

  print(pic)

  # 23: Mature - immature ratio
  N   <- jjm.out$N[,-1]
  Mat <- jjm.out$mature_a
  Wt  <- jjm.out$wt_a_pop

  MatureBiom    <- rowSums(sweep(N ,2,Mat*Wt,"*"))
  ImmatureBiom  <- rowSums(sweep(N ,2,(1-Mat)*Wt,"*"))

  res <- data.frame(rbind(cbind(jjm.out$Yr,MatureBiom,"Mature"),cbind(jjm.out$Yr,ImmatureBiom,"Immature")),stringsAsFactors=F)
  colnames(res) <- c("year","data","classing"); res$data <- as.numeric(res$data); res$year <- as.integer(res$year); res$classing <- as.factor(res$classing)

  ikey           <- simpleKey(text=c("Mature","Immature"),
                             points=F,lines=T,columns = 2)
  ikey$lines$col <- 1
  ikey$lines$lwd <- 3
  ikey$lines$lty <- c(3,1)
  ikey$points$col <- "white"

  pic <- xyplot(data ~ year,data=res,groups=classing,
          type="l",lwd=3,lty=c(1,3),col=1,
          ylab="Biomass in kt",xlab="Years",main="Mature - Immature fish",
          key=ikey,
          panel=function(...){
            panel.grid(h=-1, v= -1)
            panel.xyplot(...)
          })
  print(pic)
  
  # 24: Stock-recruitment
  res1 <- data.frame(get("jjm.out")[["Stock_Rec"]][,c(2,4)]); res1 <- res1[1:(nrow(res1)-1),]
  res1$class <- "observed"
  res2 <- data.frame(get("jjm.out")[["stock_Rec_Curve"]]); res2 <- res2[1:(nrow(res2)-1),]
  res2$class <- "modelled"
  res  <- rbind(res1,res2)
  colnames(res) <- c("SSB","Rec","class")
  
  ikey           <- simpleKey(text=c("Observed","Modelled"),
                             points=T,lines=T,columns = 2)
  ikey$lines$col <- c(1,rev(cols)[1])
  ikey$lines$lwd <- c(2,3)
  ikey$lines$lty <- c(3,2)

  ikey$points$pch <- c(19,0)
  ikey$points$col <- c("darkgrey","white")

  pic <- xyplot(Rec~SSB,data=res,groups=class,
                ylab="Recruitment",xlab="Spawning Stock Biomass",main="Stock Recruitment",
                key=ikey,
                panel=function(x,y){
                  panel.grid(h=-1, v= -1)
                  idxobs <- which(res$SSB %in% x & res$class == "observed")
                  idxmod <- which(res$SSB %in% x & res$class == "modelled")
                  panel.xyplot(x[idxobs],y[idxobs],type="l",lwd=3,col=1,lty=3)
                  panel.points(x[idxobs],y[idxobs],type="p",cex=0.6,pch=19,col="darkgrey")
                  panel.xyplot(x[idxmod],y[idxmod],type="l",lwd=5,col=rev(cols)[1],lty=2)
                })
  print(pic)
  
  # 25: SSB not fished over SSB fished
  BnoFish <- jjm.out$TotBiom_NoFish[,2]
  BFish   <- jjm.out$TotBiom[,2]
  res     <- as.data.frame(rbind(cbind(jjm.out$TotBiom[,1],BnoFish,"notfished"),
                                 cbind(jjm.out$TotBiom[,1],BFish,"fished")),stringsAsFactors=F)
  colnames(res) <- c("year","data","class")
  res$data      <- an(res$data)
  res$year      <- an(res$year)
  
  ikey           <- simpleKey(text=c("Fished","Unfished"),
                             points=F,lines=T,columns = 2)
  ikey$lines$col <- c(1,2)
  ikey$lines$lwd <- c(2,2)
  ikey$lines$lty <- c(1,1)

  pic <- xyplot(data ~ year,data=res,groups=class,
              ylab="Total biomass",xlab="Years",main="Fished vs. unfished biomass",
              key=ikey,col=c(1,2),lwd=2,lty=1,type="l",
              panel=function(...){
                panel.grid(h=-1, v= -1)
                panel.xyplot(...)
              })
  print(pic)
}
  #-----------------------------------------------------------------------------
  #- Plots of catch and ssb projections
  #-----------------------------------------------------------------------------
if("projections" %in% what){
  pic <- plot(1,1,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",main="Projections",cex.main=3)
  print(pic)

  # 25: SSB projections
  Nfutscen  <- length(grep("SSB_fut_",names(jjm.out)))
  scenarios <- c("F2012 SQ","F2012 0.75x","F2012 0.5x", "F2012 0.25x","F2012 0x")[1:Nfutscen]

  for(iScen in 1:length(scenarios)){
    idx <- nrow(get("jjm.out")[["SSB"]][,c(1,2,4,5)])
    tot <- rbind(get("jjm.out")[["SSB"]][-idx,c(1,2,4,5)],get("jjm.out")[[paste("SSB_fut_",iScen,sep="")]][,c(1,2,4,5)])
    colnames(tot) <- c("year","SSB","SSB5","SSB95")
    for(i in 2:4){
      if(iScen == 1 & i == 2){
        totres <- data.frame(cbind(tot[,1],tot[,2]))
        totres$class <- colnames(tot)[i]
        totres$scenario <- scenarios[iScen]
      }
      if(iScen != 1 | i != 2){
        res <- data.frame(cbind(tot[,1],tot[,i]))
        res$class <- colnames(tot)[i]
        res$scenario <- scenarios[iScen]
        totres <- rbind(totres,res)
      }
    }
  }
  colnames(totres) <- c("year","data","class","scenario")
  
  ikey           <- simpleKey(text=scenarios,
                             points=F,lines=T,columns = 2)
  ikey$lines$col <- 1:length(scenarios)
  ikey$lines$lwd <- 4
  ikey$lines$lty <- 1

  
#  pic <- xyplot(data ~ year,data=totres,type="l",groups=scenario,
#                ylab="Spawning Stock Biomass",xlab="Years",main="SSB prediction",
#                key=ikey,
#                prepanel=function(...) {list(ylim=c(0,max(totres$data,na.rm=T)))},
#                panel=function(x,y){
#                  panel.grid(h=-1, v= -1)
#                  idx <- mapply(seq(length(x)/Nfutscen,length(x),length.out=Nfutscen)-length(x)/Nfutscen+1,
#                          seq(length(x)/Nfutscen,length(x),length.out=Nfutscen),FUN=seq)
#                  idx2<- mapply(seq(length(idx[,1])/3,length(idx[,1]),length.out=3)-length(idx[,1])/3+1,
#                          seq(length(idx[,1])/3,length(idx[,1]),length.out=3),FUN=seq)
#                  #scen1 <- idx[,1]; scen2 <- idx[,2]; scen3 <- idx[,3]; scen4 <- idx[,4]; scen5 <- idx[,5]
#                  for(iScen in 2:Nfutscen){
#                    panel.xyplot(x[idx[,iScen][idx2[,1]]],y[idx[,iScen][idx2[,1]]],type="l",col=iScen,lwd=3)
#                    panel.xyplot(x[idx[,iScen][idx2[,2]]],y[idx[,iScen][idx2[,2]]],type="l",col=iScen,lwd=2,lty=3)
#                    panel.xyplot(x[idx[,iScen][idx2[,3]]],y[idx[,iScen][idx2[,3]]],type="l",col=iScen,lwd=2,lty=3)
#                  }
#                  panel.xyplot(x[idx[,1][idx2[,1]]],y[idx[,1][idx2[,1]]],type="l",col=1,lwd=4)
#                  panel.xyplot(x[idx[,1][idx2[,2]]],y[idx[,1][idx2[,2]]],type="l",col=1,lwd=2,lty=3)
#                  panel.xyplot(x[idx[,1][idx2[,3]]],y[idx[,1][idx2[,3]]],type="l",col=1,lwd=2,lty=3)
#                })
                
  pic <- xyplot(data ~ year,data=totres,type="l",groups=scenario,
                ylab="Spawning Stock Biomass",xlab="Years",main="SSB prediction",xlim=c(2000,max(totres$year)),
                key=ikey,
                prepanel=function(...) {list(ylim=c(0,max(totres$data,na.rm=T)))},
                panel=function(x,y){
                  panel.grid(h=-1, v= -1)
                  idx <- mapply(seq(length(x)/Nfutscen,length(x),length.out=Nfutscen)-length(x)/Nfutscen+1,
                          seq(length(x)/Nfutscen,length(x),length.out=Nfutscen),FUN=seq)
                  idx2<- mapply(seq(length(idx[,1])/3,length(idx[,1]),length.out=3)-length(idx[,1])/3+1,
                          seq(length(idx[,1])/3,length(idx[,1]),length.out=3),FUN=seq)
                  #scen1 <- idx[,1]; scen2 <- idx[,2]; scen3 <- idx[,3]; scen4 <- idx[,4]; scen5 <- idx[,5]
                  for(iScen in 2:Nfutscen){
                    panel.xyplot(x[idx[,iScen][idx2[,1]]],y[idx[,iScen][idx2[,1]]],type="l",col=iScen,lwd=3)
                    iCol  <- col2rgb(iScen)
                    iCol  <- rgb(iCol[1]/255,iCol[2]/255,iCol[3]/255,0.25)
                    panel.polygon(c(x[idx[,iScen][idx2[,2]]],rev(x[idx[,iScen][idx2[,3]]])),
                                  c(y[idx[,iScen][idx2[,2]]],rev(y[idx[,iScen][idx2[,3]]])),col=iCol,border=iCol)
                    panel.lines(x[idx[,iScen][idx2[,1]]],y[idx[,iScen][idx2[,1]]],col=iScen,lwd=3)
                  }
                  panel.xyplot(x[idx[,1][idx2[,1]]],y[idx[,1][idx2[,1]]],type="l",col=1,lwd=4)
                  iCol  <- col2rgb(1)
                  iCol  <- rgb(iCol[1]/255,iCol[2]/255,iCol[3]/255,0.15)
                  panel.polygon(c(x[idx[,1][idx2[,2]]],rev(x[idx[,1][idx2[,3]]])),
                                c(y[idx[,1][idx2[,2]]],rev(y[idx[,1][idx2[,3]]])),col=iCol,border=iCol)
                  panel.lines(x[idx[,1][idx2[,1]]],y[idx[,1][idx2[,1]]],col=1,lwd=4)
                })
  print(pic)
  
  # 26: Catch projections
  totCatch  <- 0
  for(iFlt in grep("Obs_catch_",names(jjm.out)))
    totCatch <- jjm.out[[iFlt]] + totCatch
  totCatch  <- cbind(jjm.out$Yr,totCatch); colnames(totCatch) <- c("year","catch")
  for(iScen in 1:length(scenarios)){
    tot <- rbind(totCatch,get("jjm.out")[[paste("Catch_fut_",iScen,sep="")]])
    colnames(tot) <- c("year","catch")
    if(iScen == 1){
        totres <- data.frame(tot)
        totres$scenario <- scenarios[iScen]
    } else {
        res <- data.frame(tot)
        res$scenario <- scenarios[iScen]
        totres  <- rbind(totres,res)
      }
  }
  
  colnames(totres) <- c("year","data","scenario")

  ikey           <- simpleKey(text=scenarios,
                             points=F,lines=T,columns = 2)
  ikey$lines$col <- 1:length(scenarios)
  ikey$lines$lwd <- 4
  ikey$lines$lty <- 1


  pic <- xyplot(data ~ year,data=totres,type="l",groups=scenario,
                ylab="Catch",xlab="Years",main="Catch prediction",
                key=ikey,
                prepanel=function(...) {list(ylim=c(0,max(totres$data,na.rm=T)))},
                panel=function(x,y){
                  panel.grid(h=-1, v= -1)
                  idx <- mapply(seq(length(x)/Nfutscen,length(x),length.out=Nfutscen)-length(x)/Nfutscen+1,
                          seq(length(x)/Nfutscen,length(x),length.out=Nfutscen),FUN=seq)
                  #scen1 <- idx[,1]; scen2 <- idx[,2]; scen3 <- idx[,3]; scen4 <- idx[,4]; scen5 <- idx[,5]
                  for(iScen in 2:Nfutscen) panel.xyplot(x[idx[,iScen]],y[idx[,iScen]],type="l",col=iScen,lwd=3)
                  panel.xyplot(x[idx[,1]],y[idx[,1]],type="l",col=1,lwd=4)

                })
  print(pic)

}
  #-----------------------------------------------------------------------------
  #- Plots of yield per recruit and yield biomass
  #-----------------------------------------------------------------------------
if("ypr" %in% what){
  pic <- plot(1,1,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",main="Yield & spawner per recruit",cex.main=3)
  print(pic)
  
  res <- rbind(data.frame(cbind(jjm.ypr$F,jjm.ypr$SSB),class="SSB"),
               data.frame(cbind(jjm.ypr$F,jjm.ypr$Yld),class="Yield"),
               data.frame(cbind(jjm.ypr$F,jjm.ypr$Recruit),class="Recruit"),
               data.frame(cbind(jjm.ypr$F,jjm.ypr$SPR),class="SPR"),
               data.frame(cbind(jjm.ypr$F,jjm.ypr$B),class="Biomass"),
               data.frame(cbind(jjm.ypr$F,jjm.ypr$Yld/jjm.ypr$Recruit),class="YPR"),
               data.frame(cbind(jjm.ypr$F,jjm.ypr$SSB/jjm.ypr$Recruit),class="SpawPR"))
  colnames(res) <- c("F","data","class")
  res           <- subset(res,class %in% c("YPR","SpawPR"))

  pic <- xyplot(data ~ F | class,data=res,type="l",
         main="Yield and spawing stock biomass per recruit",xlab="Fishing mortality",ylab="Spawing biomass / Yield per recruit",
         prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
         layout=c(1,2),
         panel=function(...){
            lst    <- list(...)
            panel.grid(h=-1, v= -1)
            if(panel.number()==1) panel.xyplot(lst$x,lst$y,type="l",lwd=3,lty=1,col=1)
            if(panel.number()==2) panel.xyplot(lst$x,lst$y,type="l",lwd=3,lty=1,col=1)
          },
          scales=list(alternating=1,y=list(relation="free",rot=0))
         )
  print(pic)
}

}


