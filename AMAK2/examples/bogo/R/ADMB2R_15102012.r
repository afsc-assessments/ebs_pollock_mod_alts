setOutputNames <- function(out){

          names(out)     <- c("Years","Total Fishing mortality","Total biomass no Fishing","SSB no fishing","Total biomass",
                              "SSB in the future under scenario 1","SSB in the future under scenario 2","SSB in the future under scenario 3",
                              "SSB in the future under scenario 4","SSB in the future under scenario 5","Catch in the future under scenario 1",
                              "Catch in the future under scenario 2","Catch in the future under scenario 3","Catch in the future under scenario 4",
                              "Catch in the future under scenario 5","SSB","Recruitment","Numbers at age" ,"F at age fishery 1","F at age fishery 2",
                              "F at age fishery 3","F at age fishery 4","Fishery names","Indices names","Observations at age survey 1",
                              "Observations at age survey 2","Observations at age survey 3","Observations at age survey 4","Observations at age survey 5",
                              "Observations at age survey 6","Observations at age survey 7","Observations at age survey 8","Observations at age survey 9",
                              "Survey catchabilities","Proportions at age fishery 1 observed","Proportions at age fishery 2 observed","Proportions at age fishery 3 observed",
                              "Proportions at age fishery 4 observed","Proportions at age fishery 1 predicted","Proportions at age fishery 2 predicted","Proportions at age fishery 3 predicted",
                              "Proportions at age fishery 4 predicted","Proportions at age survey 1 observed","Proportions at age survey 4 observed",
                              "Proportions at age survey 1 predicted","Proportions at age survey 4 predicted","Total catch fishery 1 observed","Total catch fishery 1 predicted",
                              "Total catch fishery 2 observed","Total catch fishery 2 predicted","Total catch fishery 3 observed","Total catch fishery 3 predicted",
                              "Total catch fishery 4 observed","Total catch fishery 4 predicted","F_fsh_1","F_fsh_2","F_fsh_3","F_fsh_4","Selectivity fishery 1",
                              "Selectivity fishery 2","Selectivity fishery 3","Selectivity fishery 4","Selectivity survey 1","Selectivity survey 2","Selectivity survey 3",
                              "Selectivity survey 4","Selectivity survey 5","Selectivity survey 6","Selectivity survey 7","Selectivity survey 8","Selectivity survey 9",
                              "Stock recruitment","Stock recruitment curve","Likelihood composition","Likelihood composition names","Sel_Fshry_1","Sel_Fshry_2",
                              "Sel_Fshry_3","Sel_Fshry_4","Survey_Index_1","Survey_Index_2","Survey_Index_3","Survey_Index_4","Survey_Index_5","Survey_Index_6","Survey_Index_7",
                              "Survey_Index_8","Survey_Index_9","Age_Survey_1","Age_Survey_2","Age_Survey_3","Age_Survey_4","Age_Survey_5","Age_Survey_6","Age_Survey_7",
                              "Age_Survey_8","Age_Survey_9","Sel_Survey_1","Sel_Survey_2","Sel_Survey_3","Sel_Survey_4","Sel_Survey_5","Sel_Survey_6","Sel_Survey_7","Sel_Survey_8","Sel_Survey_9",
                              "Recruitment penalty","F penalty","Survey 1 catchability penalty","Survey 1 catchability power function","Survey 2 catchability penalty","Survey 2 catchability power function",
                              "Survey 3 catchability penalty","Survey 3 catchability power function","Survey 4 catchability penalty","Survey 4 catchability power function",
                              "Survey 5 catchability penalty","Survey 5 catchability power function","Survey 6 catchability penalty","Survey 6 catchability power function",
                              "Survey 7 catchability penalty","Survey 7 catchability power function","Survey 8 catchability penalty","Survey 8 catchability power function",
                              "Survey 9 catchability penalty","Survey 9 catchability power function","Natural mortality","Steepness of recruitment","Sigma recruitment","Number of parameters estimated",
                              "Steepness prior","Sigma recruitment prior","Rec_estimated_in_styr_endyr","SR_Curve_fit__in_styr_endyr","Model_styr_endyr","Natural mortality prior","q prior",
                              "q power prior","cv catch biomass","Projection year range","Fsh_sel_opt_fish","Survey_Sel_Opt_Survey","Phase_survey_Sel_Coffs","Fishery selectivity ages",
                              "Survey selectivity ages","Phase_for_age_spec_fishery","Phase_for_logistic_fishery","Phase_for_dble_logistic_fishery","Phase_for_age_spec_survey","Phase_for_logistic_survey",
                              "Phase_for_dble_logistic_srvy","EffN_Fsh_1","EffN_Fsh_2","EffN_Fsh_3","EffN_Fsh_4","C_fsh_1","C_fsh_2","C_fsh_3","C_fsh_4","Weight at age in the population",
                              "Maturity at age","Weight at age in fishery 1","Weight at age in fishery 2","Weight at age in fishery 3","Weight at age in fishery 4",
                              "Weight at age in survey 1","Weight at age in survey 2","Weight at age in survey 3","Weight at age in survey 4","Weight at age in survey 5",
                              "Weight at age in survey 6","Weight at age in survey 7","Weight at age in survey 8","Weight at age in survey 9","EffN_Survey_1","EffN_Survey_4")
                    return(out)}
                    
readYPR     <- function(fileName){
                 jjm.ypr            <- read.table(fileName,sep=" ",skip=4,header=T,fill=T);
                 jjm.ypr[1,]        <- jjm.ypr[1,c(1,8,2,3,4,5,6,7)]
                 colnames(jjm.ypr)  <- c("F","X","SSB","Yld","Recruit","SPR","B","X2");
                 jjm.ypr            <- jjm.ypr[,-grep("X",colnames(jjm.ypr))]
              return(jjm.ypr)}
              
              
an          <- function(x){return(as.numeric(x))}
ac          <- function(x){return(as.character(x))}
anf         <- function(x){return(as.numeric(as.character(x)))}
              
#-Code to read in final data
read.dat <- function(iFilename,iPath){
                ###-Read in the raw datafile-###
                res1      <- scan(file=paste(iPath,iFilename,sep=""),what='numeric', quiet=TRUE,sep="\n",comment.char="#",allowEscapes=T)
                res1      <- strsplit(res1,"\t")

                #- Get some initial dimensions
                nY        <- length(an(unlist(res1[[1]][1])):an(unlist(res1[[2]][1]))) #number of years
                Ys        <- an(unlist(res1[1:2])) #Years
                nA        <- length(an(unlist(res1[[3]][1])):an(unlist(res1[[4]][1]))) #number of ages
                As        <- an(unlist(res1[3:4])) #Ages
                Ls        <- an(unlist(res1[5:6])) #Lengths
                nL        <- na.omit(an(unlist(res1[[8]]))) #number of lengths
                nF        <- na.omit(an(unlist(res1[[9]]))) #number of fisheries

                #- Define storage object
                cols      <- list()

                ###-Fill cols with data from res1-###

                #-Common data
                cols$years                =matrix(NA,ncol=2,nrow=1 ,dimnames=list("years",c("first year","last year")))
                cols$years[]                          <- an(unlist(res1[1:2]))
                cols$ages                 =matrix(NA,ncol=2,nrow=1 ,dimnames=list("age",c("age recruit","oldest age")))
                cols$ages[]                           <- an(unlist(res1[3:4]))
                cols$lengths              =matrix(NA,ncol=2,nrow=1 ,dimnames=list("lengths",c("first length","last length")))
                cols$lengths[]                        <- an(unlist(res1[5:6]))
                cols$lengthbin            =numeric()
                cols$lengthbin                        <- an(unlist(res1[7]))

                #-Fisheries data
                cols$Fnum                 =numeric()
                cols$Fnum                             <- na.omit(an(unlist(res1[9])))
                
                #-Start of dynamic rows
                counter                               <- 10 #first dynamic row
                
                cols$Fnames               =list()
                cols$Fnames                           <- strsplit(unlist(res1[counter]),"%")[[1]];                                  counter <- counter + 1
                cols$Fcaton               =matrix(NA,ncol=nF,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("fishery",1:nF,sep="")))
                cols$Fcaton[]                         <- matrix(na.omit(an(unlist(res1[counter:(counter+nF-1)]))),ncol=nF,nrow=nY); counter <- counter + nF
                cols$Fcatonerr            =matrix(NA,ncol=nF,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("fishery",1:nF,sep="")))
                cols$Fcatonerr[]                      <- matrix(na.omit(an(unlist(res1[counter:(counter+nF-1)]))),ncol=nF,nrow=nY); counter <- counter + nF
                cols$FnumyearsA           =matrix(NA,ncol=nF,nrow=1,dimnames=list("years",paste("Fyears",1:nF,sep="")))
                cols$FnumyearsA[]                      <- na.omit(an(unlist(res1[counter:(counter+nF-1)])));                        counter <- counter + nF
                cols$FnumyearsL           =matrix(NA,ncol=nF,nrow=1,dimnames=list("years",paste("Fyears",1:nF,sep="")))
                cols$FnumyearsL[]                   <- na.omit(an(unlist(res1[counter:(counter+nF-1)])));                           counter <- counter + nF

                cols$Fageyears            =matrix(NA,ncol=nF,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("fishery",1:nF,sep="")))
                for(iFs in 1:nF){
                  if(cols$FnumyearsA[iFs]>0){
                    Fageyears <- c(na.omit(an(res1[[counter]])))
                    wFyears   <- pmatch(Fageyears,cols$years[1]:cols$years[2])
                    cols$Fageyears[wFyears,paste("fishery",iFs,sep="")]         <- Fageyears
                    counter   <- counter +1
                  }
                }
                cols$Flengthyears         =matrix(NA,ncol=nF,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("fishery",1:nF,sep="")))
                for(iFs in 1:nF){
                  if(cols$FnumyearsL[iFs]>0){
                    Flengthyears  <- c(na.omit(an(res1[[counter]])))
                    lFyears       <- pmatch(Flengthyears,cols$years[1]:cols$years[2])
                    cols$Flengthyears[lFyears,paste("fishery",iFs,sep="")]      <- Flengthyears
                    counter       <- counter +1
                  }
                }

                cols$Fagesample           =matrix(NA,ncol=nF,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("fishery",1:nF,sep="")))
                for(iFs in 1:nF){
                  if(cols$FnumyearsA[iFs]>0){
                    wFyears       <- rownames(cols$Fageyears)[which(is.na(cols$Fageyears[,paste("fishery",iFs,sep="")])==F)]
                    cols$Fagesample[wFyears,paste("fishery",iFs,sep="")]        <- na.omit(an(unlist(res1[counter])));
                    counter       <- counter + 1
                  }
                }

                cols$Flengthsample        =matrix(NA,ncol=nF,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("fishery",1:nF,sep="")))
                for(iFs in 1:nF){
                  if(cols$FnumyearsL[iFs]>0){
                    lFyears       <- rownames(cols$Flengthyears)[which(is.na(cols$Flengthyears[,paste("fishery",iFs,sep="")])==F)]
                    cols$Flengthsample[lFyears,paste("fishery",iFs,sep="")]        <- na.omit(an(unlist(res1[counter])));
                    counter       <- counter + 1
                  }
                }

                cols$Fagecomp             =array (NA,dim=c(nY,nA,nF),dimnames=list(years=Ys[1]:Ys[2],age=As[1]:As[2],paste("fishery",1:nF,sep="")))
                for(iFs in 1:nF){
                  if(cols$FnumyearsA[iFs]>0){
                    wFyears       <- rownames(cols$Fageyears)[which(is.na(cols$Fageyears[,paste("fishery",iFs,sep="")])==F)]
                    cols$Fagecomp[wFyears,,paste("fishery",iFs,sep="")]      <- matrix(na.omit(an(unlist(res1[counter:(counter+length(wFyears)-1)]))),ncol=nA,
                                                                                       nrow=length(wFyears),byrow=T)
                    counter       <- counter+length(wFyears)
                  }
                }

                cols$Flengthcomp          =array (NA,dim=c(nY,nL,nF),dimnames=list(years=Ys[1]:Ys[2],lengths=Ls[1]:Ls[2],paste("fishery",1:nF,sep="")))
                for(iFs in 1:nF){
                  if(cols$FnumyearsL[iFs]>0){
                    lFyears       <- rownames(cols$Flengthyears)[which(is.na(cols$Flengthyears[,paste("fishery",iFs,sep="")])==F)]
                    cols$Flengthcomp[lFyears,,paste("fishery",iFs,sep="")]      <- matrix(na.omit(an(unlist(res1[counter:(counter+length(lFyears)-1)]))),ncol=nL,
                                                                                       nrow=length(lFyears),byrow=T)
                    counter       <- counter+length(lFyears)
                  }
                }

                cols$Fwtatage              =array (NA,dim=c(nY,nA,nF),dimnames=list(years=Ys[1]:Ys[2],age=As[1]:As[2],paste("fishery",1:nF,sep="")))
                for(iFs in 1:nF){
                  cols$Fwtatage[,,iFs]                <-matrix(na.omit(an(unlist(res1[counter:(counter+nY-1)]))),ncol=nA,
                                                               nrow=nY,byrow=T)
                  counter                             <- counter + nY
                }



                #-Indices data
                nI                                    <- na.omit(an(res1[[counter]]))
                cols$Inum                 =numeric()
                cols$Inum                             <- na.omit(an(res1[[counter]]));                                  counter <- counter + 1
                cols$Inames               =list()
                cols$Inames                           <- strsplit(res1[[counter]],"%")[[1]];                            counter <- counter + 1
                cols$Inumyears            =matrix(NA,ncol=nI,nrow=1,dimnames=list("years",paste("index",1:nI,sep="")))
                cols$Inumyears[]                      <- na.omit(an(unlist(res1[counter:(counter+cols$Inum-1)])));      counter <- counter + cols$Inum

                cols$Iyears               =matrix(NA,ncol=nI,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("index",1:nI,sep="")))
                for(iSu in 1:nI){
                  if(cols$Inumyears[iSu]>0){
                    Iyears            <- na.omit(an(res1[[counter]]));  wIyears <- pmatch(Iyears,cols$years[1]:cols$years[2])
                    cols$Iyears[wIyears,paste("index",iSu,sep="")]              <- Iyears
                    counter           <- counter + 1
                  }
                }
                    
                cols$Imonths              =matrix(NA,ncol=nI,nrow=1,dimnames=list("month",paste("index",1:nI,sep="")))
                cols$Imonths[]                        <- na.omit(an(unlist(res1[counter:(counter+cols$Inum-1)]))); counter <- counter + cols$Inum
                cols$Index                =matrix(NA,ncol=nI,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("index",1:nI,sep="")))
                for(iSu in 1:nI){
                  if(cols$Inumyears[iSu]>0){
                    wIyears       <- rownames(cols$Iyears)[which(is.na(cols$Iyears[,paste("index",iSu,sep="")])==F)]
                    cols$Index[wIyears,paste("index",iSu,sep="")] <- na.omit(an(res1[[counter]]))
                    counter       <- counter + 1
                  }
                }

                cols$Indexerr             =matrix(NA,ncol=nI,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("index",1:nI,sep="")))
                for(iSu in 1:nI){
                  if(cols$Inumyears[iSu]>0){
                    wIyears       <- rownames(cols$Iyears)[which(is.na(cols$Iyears[,paste("index",iSu,sep="")])==F)]
                    cols$Indexerr[wIyears,paste("index",iSu,sep="")] <- na.omit(an(res1[[counter]]))
                    counter       <- counter + 1
                  }
                }
                
                cols$Inumageyears         =matrix(NA,ncol=nI,nrow=1,dimnames=list("years",paste("index",1:nI,sep="")))
                cols$Inumageyears[] <- na.omit(an(unlist(res1[counter:(counter+cols$Inum-1)]))); counter <- counter + cols$Inum

                cols$Iyearsage            =matrix(NA,ncol=nI,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("index",1:nI,sep="")))
                for(iSu in 1:nI){
                  if(cols$Inumageyears[iSu]>0){
                    Iyearsage         <- na.omit(an(res1[[counter]])); wIyearsage <- pmatch(Iyearsage,cols$years[1]:cols$years[2])
                    cols$Iyearsage[wIyearsage,iSu] <- Iyearsage
                    counter           <- counter + 1
                  }
                }

                cols$Iagesample           =matrix(NA,ncol=nI,nrow=nY,dimnames=list(years=Ys[1]:Ys[2],paste("index",1:nI,sep="")))
                for(iSu in 1:nI){
                  if(cols$Inumageyears[iSu]>0){
                    wIyears       <- rownames(cols$Iyearsage)[which(is.na(cols$Iyearsage[,paste("index",iSu,sep="")])==F)]
                    cols$Iagesample[wIyears,iSu] <- na.omit(an(res1[[counter]]))
                    counter       <- counter + 1
                  }
                }
                cols$Ipropage             =array (NA,dim=c(nY,nA,nI),dimnames=list(years=Ys[1]:Ys[2],age=As[1]:As[2],paste("index",1:nI,sep="")))
                for(iSu in 1:nI){
                  if(cols$Inumageyears[iSu]>0){
                    wIyears       <- rownames(cols$Iyearsage)[which(is.na(cols$Iyearsage[,paste("index",iSu,sep="")])==F)]
                    cols$Ipropage[wIyears,,iSu]   <- matrix(na.omit(an(unlist(res1[counter:(counter+cols$Inumageyears[iSu]-1)]))),ncol=nA,
                                                                nrow=cols$Inumageyears[iSu],byrow=T)
                    counter       <- counter + cols$Inumageyears[iSu]
                  }
                }

                cols$Iwtatage             =array (NA,dim=c(nY,nA,nI),dimnames=list(years=Ys[1]:Ys[2],age=As[1]:As[2],paste("index",1:nI,sep="")))

                for(iSu in 1:nI){
                  cols$Iwtatage[,,iSu]                <-matrix(na.omit(an(unlist(res1[counter:(counter+nY-1)]))),ncol=nA,
                                                               nrow=nY,byrow=T)
                  counter                             <- counter + nY
                }

                #-Population data
                cols$Pwtatage             =matrix(NA,ncol=1,nrow=nA,dimnames=list(age=As[1]:As[2],"weight"))
                cols$Pwtatage[]                       <- na.omit(an(res1[[counter]]));      counter <- counter + 1
                cols$Pmatatage            =matrix(NA,ncol=1,nrow=nA,dimnames=list(age=As[1]:As[2],"maturity"))
                cols$Pmatatage[]                      <- na.omit(an(res1[[counter]]));      counter <- counter + 1
                cols$Pspwn                =numeric()
                cols$Pspwn                            <- na.omit(an(res1[[counter]]));      counter <- counter + 1
                cols$Pageerr              =matrix(NA,ncol=nA,nrow=nA,dimnames=list(age=As[1]:As[2],age=As[1]:As[2]))
                cols$Pageerr[]                        <- matrix(na.omit(an(unlist(res1[counter:(counter+nA-1)]))),ncol=nA,nrow=nA,byrow=T); counter <- counter + nA
             return(cols)}