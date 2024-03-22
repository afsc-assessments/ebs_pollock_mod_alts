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

read.psv <- function(fn, nsamples=10000) {
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
# amak.dat  
write_amak_dat <- function(retro=0,surv_dwnwt=0,lsryr=2017,nsel=41){
 cat (file="amak.dat","am2019b.dat \n Model_16.0b 
#selectivity_shar_matrix 
1 2  
1 1 
#Sr_type 
2  
#AgeError 
1  
#Retro \n ",
retro, 
"\n #surv_dwnwt \n",
surv_dwnwt, " 
#cv_inc
100
#n_comp
0.01
#Steepness 
0.8  300  -6  
#SigmaR 
0.6  15  4  
#yrs_sr 
1977  ", lsryr, "\n
#Linf
74.4  0.1 -4                                                                                  
#K
0.16  0.1 -4                                                                                  
#Lo_Len
18  0.1 -4                                                                                  
#Sigma_len
0.09 0.1  -4                                                                                  
#Natural_Mortality 
0.3  0.05  -4  
# NEW npars_mage
0
# NEW Mage_in
# phase_Mage
-5
#Nyrs_Random_walk_M 
0
#Random_walk_M_yrs blank if nyrs==0
#Random_walk_M_sigmas blank if nyrs==0
#Random_walk_q_phases 
-4
#catchability 
1  0.2  4  
#q_power 
1  0.2  -4  
#Random_walk_q_phases 
-1
#Nyrs_Random_walk_q
0
#Random_walk_q_yrs blank if nyrs==0

#Random_walk_q_sigmas blank if nyrs==0

#q_agemin 
4  
#q_agemax 
10  
#junk 
0.05  
#n_proj_yrs 
10  
#Fsh_selopt_1 
1  
#Fsh_nages_1 
10  
#Fsh_ph_1 
4  
#Fsh_curvpen_1 
0.946 
#Fsh_domepen_1 
10.3  
#Fsh_sel_change_1 
# nyrs fish selectivity changes
",
nsel, "
# Years of selectivity change fishery 
",
seq(1978,1978+nsel-1),"\n
# sigm fish selectivity changes
", rep(0.35,nsel),"\n",
"
#Fsh_sel_init_1 
0.02 0.7 1 1 1 1 1 1 1 1 1 
#Ind_selopt 
1  
# Selages
10  
# Phase Survey
5  
# sigma age-age
0.5  
# sigma Dome... 
0.4  
# Years of selectivity change survey 
# nyrs fish selectivity changes
0
# Yrs fish selectivity changes
#1987
# sigm fish selectivity changes
#0.7 
# Initial values for coefficitients at each change (one for every change plus 1) 
# 2 3 4 5 6 7 8 9 10 11 12  
0.2 0.7 1 1 1 1 1 1 1 1 1 
#Test  
123456789  
"
  )
}

write_proj_setup<- function(begyr=2019){
 cat (file="setup.dat","
# a new file
#Run_name
Std
#Tier
3
#nalts
1
#alts
3
#tac_abc
1
#srr
1
#rec_proj
1
#srr_cond
0
#srr_prior
0
#write_big
1
#nyrs_proj
3
#nsims
1000
#beg_yr_label
",
begyr)
}

write_spp_catch<- function(catch=9999,curyear){
	cat (file="spp_catch.dat","
#_Number_of_years with specified catch 
1
# Number of species                                         
1 
# OY Minimum                      
1343.248  # Note that this is for age-structured species    1330.148                
# OY Maximum                      
1943.248  # Note that this is for age-structured species    1930.148                
# data files for each species                 
amak.prj
#ABC Multipliers                                           
1
# scalars                   
1
0.6
  # Number of TAC model categories                
1 
  # TAC model indices (for aggregating)         (should be 1)
1
# Catch in each future year           
",
curyear,
catch
		)
}
#write_spp_catch(curyear=1999)