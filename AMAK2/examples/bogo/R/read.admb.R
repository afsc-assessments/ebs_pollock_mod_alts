read.admb <-
function(ifile)
{	
	ret=read.fit(ifile)
	
	fn=paste(ifile,'.rep', sep='')
	A=read.rep(fn)
	A$fit=ret
	
	pfn=paste(ifile,'.psv',sep='')
	if(file.exists(pfn))
		A$post.samp=read.psv(pfn)
	
	return(A)
}

read.fit <-
function(ifile)
{
	# __Example:             
	#	file <-("~/admb/simple")
	#	A <- reptoRlist(file)
	#	Note there is no extension on the file name.
	
	## The following is a contribution from:
	## Anders Nielsen that reads the par & cor files.
	ret<-list() 
	parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
	 what='', n=16, quiet=TRUE)[c(6,11,16)]) 
	ret$nopar<-as.integer(parfile[1]) 
	ret$nlogl<-parfile[2] 
	ret$maxgrad<-parfile[3] 
	file<-paste(ifile,'.cor', sep='') 
	lin<-readLines(file) 
	ret$npar<-length(lin)-2 
	ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
	sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
	ret$names<-unlist(lapply(sublin,function(x)x[2])) 
	ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
	ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
	ret$cor<-matrix(NA, ret$npar, ret$npar) 
	corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
	ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
	ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
	ret$cov<-ret$cor*(ret$std%o%ret$std)
	return(ret)
}

read.rep <- 
function(fn)
{
	# The following reads a report file
	# Then the 'A' object contains a list structure
	# with all the elemements in the report file.
	# In the REPORT_SECTION of the AMDB template use 
	# the following format to output objects:
	#  	report<<"object \n"<<object<<endl;
	#
	# The part in quotations becomes the list name.
	# Created By Steven Martell
	options(warn=-1)  #Suppress the NA message in the coercion to double
	
	
	ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
	idx=sapply(as.double(ifile),is.na)
	vnam=ifile[idx] #list names
	nv=length(vnam) #number of objects
	A=list()
	ir=0
	for(i in 1:nv)
	{
		ir=match(vnam[i],ifile)
		if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
		dum=NA
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
		if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))

		if(is.numeric(dum))#Logical test to ensure dealing with numbers
		{
			A[[vnam[i]]]=dum
		}
	}
	options(warn=0)
	
	return(A)
}

read.dat <- 
  function(fn)
  {
    # The following reads a report file
    # Then the 'A' object contains a list structure
    # with all the elemements in the report file.
    # In the REPORT_SECTION of the AMDB template use 
    # the following format to output objects:
    #  	report<<"object \n"<<object<<endl;
    #
    # The part in quotations becomes the list name.
    # Created By Steven Martell
    options(warn=-1)  #Suppress the NA message in the coercion to double
    ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
    idx=sapply(as.double(ifile),is.na)
    
    ifile[idx]=substr(ifile[idx],2,20) #list names but strip first character (# sign)
    vnam=ifile[idx] #list names but strip first character (# sign)
    nv=length(vnam) #number of objects
    print(vnam)
    A=list()
    ir=0
    for(i in 1:nv)
    {
      ir=match(vnam[i],ifile)
      if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
      dum=NA
      if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
      if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))
      
      if(is.numeric(dum))#Logical test to ensure dealing with numbers
      {
        A[[vnam[i]]]=dum
      }
    }
    options(warn=0)
    
    return(A)
  }

read.psv <-
function(fn, nsamples=10000)
{
	#This function reads the binary output from ADMB
	#-mcsave command line option.
	#fn = paste(ifile,'.psv',sep='')
	filen <- file(fn, "rb")
	nopar <- readBin(filen, what = integer(), n = 1)
	mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
	mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
	close(filen)
	return(mcmc)
}

gletter <-
function(i=1)
{
	usr=par("usr"); inset.x=0.05*(usr[2]-usr[1]); inset.y=0.05*(usr[4]-usr[3])
	text(usr[1]+inset.x,usr[4]-inset.y,paste("(",letters[i],")",sep=""),cex=1.,font=1)
}

WriteDat <- function (d, i1,i2, filen) {
  for (i in i1:i2) {
    cat(paste("\n #",names(d[i]),"\n",sep=""),file=filen,append=T)
    #cat(d[[i]],file=filen,append=T)
    write.table(d[[i]],file=filen,append=T,row.names=F,col.names=F)
  }
}

R2_admb_input<-function(d,filen,name="Run ",retro=0,removeall=T){
  objnames=names(d)
  len=length(objnames)
  termyr=d$endyr
  nyrs=d$endyr-d$styr+1-retro
  #Option to remove data only but keep total catch...
  if (removeall){
    d$endyr    = d$endyr-retro
    d$catch    = d$catch[1:nyrs]
    d$catch_cv = d$catch_cv[1:nyrs]
  }
  for (i in 1:retro){ #Pare this down till you get right set
    if (max(d$yrs_ages_fsh)>(termyr-i)) { 
      d$n_ages_fsh     = d$n_ages_fsh-1
      d$yrs_ages_fsh   = d$yrs_ages_fsh[1:d$n_ages_fsh]
      d$sample_ages_fsh= d$sample_ages_fsh[1:d$n_ages_fsh]
      d$page_fsh       = d$page_fsh[1:d$n_ages_fsh,]
    }
    if (max(d$yrs_ind)>(termyr-i))  {
      d$nobs_ind = d$nobs_ind-1
      d$yrs_ind  = d$yrs_ind[1:d$nobs_ind]
      d$biom_ind = d$biom_ind[1:d$nobs_ind]
      d$biom_cv  = d$biom_cv[1:d$nobs_ind]
    }
    if (max(d$yrs_ages_ind)>(termyr-i))  { 
      d$n_ages_ind      = d$n_ages_ind-1
      d$yrs_ages_ind    = d$yrs_ages_ind[1:d$n_ages_ind]
      d$sample_ages_ind = d$sample_ages_ind[1:d$n_ages_ind]
      d$page_ind        = d$page_ind[1:d$n_ages_ind,]
    }    
    if (removeall){
      d$wt_age_ind=d$wt_age_ind[1:nyrs,]
      d$wt_age_fsh=d$wt_age_fsh[1:nyrs,]
    }
  }  
  cat(paste("#","Retro=",retro,"\n"),file=filen)
  WriteDat(d,i1=1,i2=5,filen)
  cat("#FisheryName\n Atka",file=filen,append=T)
  WriteDat(d,i1=6,i2=13,filen)
  cat("#SurveyName\n NMFS_Bottom_Trawl",file=filen,append=T)
  WriteDat(d,i1=14,i2=len,filen)
  
}