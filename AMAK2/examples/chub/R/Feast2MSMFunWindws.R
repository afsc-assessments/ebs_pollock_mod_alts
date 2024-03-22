##Feast2MSMFun.R
## Functions for Feast 2 MSM R file

## read each of the feast txt files
FeastSlurpee<-function(wd,txtfilen,spp,file.yrs){
	nspp<-length(spp)
	fnall<-paste(wd,txtfilen,sep="")
	dat<-list()
	for(k in 1:nspp){
		tmp<-NA
		sp.file<-fnall[k]
		tmp<-FeastImport(file.yrs=file.yrs,sp.file=fnall[k])
		eval(parse(text=paste(spp[k],"<-tmp",sep="")))
		dat[[spp[k]]]<-tmp
	}
	return(dat)
} # end function

## import Feast txt files

FeastImport<-function(file.yrs=file.yrs,sp.file=sp.file){
	options("warn"=-1)
	for(y in 1:length(file.yrs)){
		if(y==1)
			A <- list()
		fn<-paste(sp.file,file.yrs[y],".txt",sep="")	
		ifile <- scan(fn, what="character", skip=0,flush=T,blank.lines.skip=FALSE, quiet=T)
		nlines<-length(ifile)
		ll<-1:nlines
		iflex<-which(is.na(ifile))
		idx <- sapply(as.double(ifile), is.na) # text lines
		vnam2 <- ifile[idx] #list names
		vnamt<-strsplit(vnam2,split="#")
		for(i in 1:length(vnam2))
			vnam2[i]<-strsplit(vnamt[[i]][2],split=":")[[1]][1]
		datnum<-which(idx==FALSE)
		nv <- length(vnam2) #number of objects
		ir <- 0
		yrdum<-paste("yrs_",vnam2,sep="")
		ndum<-paste("nyr_",vnam2,sep="")
		vnum<-ll[idx] # lines where text
		if(y==1)
			for (i in 1:(nv)){
				if(i>1)
					A[[ndum[i]]]<-"place.holder"
				if(i>2)
					A[[yrdum[i]]]<-"place.holder"
				A[[vnam2[i]]]<-"place.holder"	
			}		
		for (i in 1:(nv)){
			ir<-vnum[i]# first line of text
			SVM<-0  # if 0 then is a scaler,if 1 then is a vector, if M then is a matrix 
			dum<-NA
			ndumVal<-NA
			yrdumVal<-NA
			if(i<=2){
				dum <- as.double(scan(fn, skip=ir, nlines=1,quiet=TRUE, what=""))
				if(i==1)
					SVM<-0
				if(i==2){
					SVM<-1
					ndumVal<-1
				}	
			}else{
				dum1 <- as.numeric(na.omit(as.double(scan(fn, skip=ir, nlines=1,quiet=TRUE, what="")))) # scan the first data line for 0 or 1
				if(dum1==1) # if there is data
				{
					
					ndumVal<-1
					yrdumVal<-file.yrs[y]
					irr <- irr<-vnum[i+1]
					if(i==nv)
						irr <- length(ifile)+1 #next row
					dum <- NA
					if (is.na(irr))
						dum <- 000
					if (irr-ir==3)
						dum <- as.double(scan(fn, skip=ir+1, nlines=1,quiet=TRUE, what=""))
						SVM<-2
						if(length(dum)==1)
							SVM<-1
					if (irr-ir>3){
						SVM<-3  # make an array (nyr,rows,cols)
						nall<-(irr-ir-2)
						tmpl<-rep(0,nall)
						nrow<-nall
						for(j in 1:length(tmpl))
							tmpl[j] <- length(as.double(scan(fn, skip=ir+j, nlines=1,quiet=TRUE, what="")))
						if(all(tmpl==tmpl[1])){
							dum<-as.matrix(read.table(fn, skip=ir+1,nrow=nall, fill=TRUE))
						}else{	
							#dum<-list()
							#for(s in 1:nspp)
							#	dum[[s]]<-as.matrix(read.table(fn, skip=ir+1+(nrow*s-nrow),nrow=nrow, fill=TRUE))
						}
					}
				}	
			} #end if data exists	
			if(y==1){
				if(is.na(ndumVal)==FALSE)
					A[[ndum[i]]]<-ndumVal
				if(is.na(yrdumVal)==FALSE)
					A[[yrdum[i]]]<-yrdumVal		
				if(is.na(dum[1])==FALSE)	
					A[[vnam2[i]]] <- dum
				if(SVM==3){	#assign nages and nlengths to the matrix
					A[[paste("nrows_",vnam2[i],sep="")]]<- dim(dum)[1]	
					A[[paste("ncols_",vnam2[i],sep="")]]<- dim(dum)[2]		
				}	
			}else{
				if(SVM==1){
					if(is.na(ndumVal)==FALSE)
						A[[ndum[i]]]<-if(A[[ndum[i]]]=="place.holder"){ndumVal}else{A[[ndum[i]]]+ndumVal}
						#A[[ndum[i]]]<-A[[ndum[i]]]+ndumVal
					if(is.na(yrdumVal)==FALSE)
						A[[yrdum[i]]]<-if(A[[yrdum[i]]]=="place.holder"){yrdumVal}else{c(A[[yrdum[i]]],yrdumVal)}
						#A[[yrdum[i]]]<-c(A[[yrdum[i]]],yrdumVal)	
					if(is.na(dum[1])==FALSE)
						A[[vnam2[i]]]<-if(A[[vnam2[i]]]=="place.holder"){dum}else{c(A[[vnam2[i]]],dum)}
						#A[[vnam2[i]]]<-c(A[[vnam2[i]]],dum)
				}
				if(SVM==2|SVM==3){
					if(is.na(ndumVal)==FALSE)
						A[[ndum[i]]]<-if(A[[ndum[i]]]=="place.holder"){ndumVal}else{A[[ndum[i]]]+ndumVal}
						#A[[ndum[i]]]<-A[[ndum[i]]]+ndumVal
						#A[[ndum[i]]]<-A[[ndum[i]]]+ndumVal
					if(is.na(yrdumVal)==FALSE)
						A[[yrdum[i]]]<-if(A[[yrdum[i]]]=="place.holder"){yrdumVal}else{c(A[[yrdum[i]]],yrdumVal)}
						#A[[yrdum[i]]]<-c(A[[yrdum[i]]],yrdumVal)
						#A[[yrdum[i]]]<-c(A[[yrdum[i]]],yrdumVal)	
					if(is.na(dum[1])==FALSE){	
						if(SVM==2)
							A[[vnam2[i]]]<-if(A[[vnam2[i]]]=="place.holder"){dum}else{rbind(A[[vnam2[i]]],as.numeric(dum))}
							#A[[vnam2[i]]]<-rbind(A[[vnam2[i]]],as.numeric(dum))	
						if(SVM==3)
							A[[vnam2[i]]]<-if(A[[vnam2[i]]]=="place.holder"){dum}else{rbind(A[[vnam2[i]]],(dum))}
							#A[[vnam2[i]]]<-c(A[[vnam2[i]]],dum)
							#A[[vnam2[i]]]<-rbind(A[[vnam2[i]]],(dum))		
					}		
							
				}
			}# end make A
		}# end nv
	}# end y
	return(A)
	options("warn"=0)
	rm(list=ls())
}	








### reconfigure the FEAST A outputs:
reconfig<-function(dat=dataA){
	Surveytypes<-unique(dat$Survey_type)
	nd<-length(Surveytypes)
	dd<-dim(dat)
	tmptxt1<-paste("newdat<-list(",Surveytypes[1],"=",Surveytypes[1],sep="")
	for(i in 1:nd){
		sub.dat<-dat[dat$Survey_type==Surveytypes[i],]
		datatypes<-unique(sub.dat$Data_type)
		nds<-(tapply(sub.dat$Data_type,sub.dat$Data_type,length))
		if(any(is.na(as.numeric(nds))))
			nds<-nds[-which(is.na(as.numeric(nds)))]
		u.nds<-unique(nds)
		tmp.txt<-paste(Surveytypes[i],"<-list(new.dat1=new.dat1",sep="")
		for(j in 1:length(u.nds)){
			sub.types<-(names(nds)[nds==u.nds[j]])
			nst<-length(sub.types)
			for(jj in 1:nst){
				sub.sub<-sub.dat[sub.dat$Data_type==(sub.types[jj]),]
				if(jj==1)
					eval(parse(text=paste("new.dat",j,"<-data.frame(sub.sub[,1:(dd[2]-2)])",sep="")))
					
				eval(parse(text=paste("new.dat",j,"$",sub.types[jj],"<-sub.sub[,(dd[2])]",sep="")))					
			}
			if(j>1)
				tmp.txt<-paste(tmp.txt,",new.dat",j,"=","new.dat",j,sep="")
		}	
		eval(parse(text=paste(tmp.txt,")",sep="")))
		tmptxt1<-paste(tmptxt1,",",Surveytypes[i],"=",Surveytypes[i],sep="")
	}
	eval(parse(text=paste(tmptxt1,")",sep="")))
	# output list with the data types as named objects
	return(newdat)
}

## write .dat file
write.dat<-function(dat2wrt,datfilename="test.dat",spname=""){
	dat<-dat2wrt
	outfile<-datfilename
	if(file.exists(outfile)){0}else{file.create(outfile)}
	nd<-length(dat)
	i<-1
	labs<-paste("#",names(dat)," :",sep="")
	cat(paste("# .dat file generated from Feast outputs for",spname,format(Sys.time(), "%b %d %Y")),file=outfile,append=FALSE,sep="\n")
	cat("#______________________________________________________________",file=outfile,append=TRUE,sep="\n")
	
	for(i in 1:nd){
		cat(labs[i],file=outfile,append=TRUE,sep="\n")
		tmp<-dat[[i]]
		if(is.list(tmp)){
			ld<-length(tmp)
			for(s in 1:ld){
				dd<-dim(tmp[[s]])
				for(r in 1:dd[1]){
					cat(dat[[i]][[s]][r,],file=outfile,append=TRUE,sep=" ");cat("",file=outfile,append=TRUE,sep="\n")
				}
			}
		}else{
			if(is.null(dim(tmp))){
				cat(dat[[i]],file=outfile,append=TRUE,sep=" ");cat("",file=outfile,append=TRUE,sep="\n")
			}else{
				dd<-dim(tmp)
				for(r in 1:dd[1]){
					cat(dat[[i]][r,],file=outfile,append=TRUE,sep=" ");cat("",file=outfile,append=TRUE,sep="\n")
				}
			}	
		}
	}
}



