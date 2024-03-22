get_DSS<-function(species=201,area="AI"){
   if(area=="AI"){
    reg<-"539 and 544"
  }
  if(area=="GOA"){
    reg<-"600 and 699"
  }
  if(area=="BS"){
    reg<-"500 and 539"
  }

  test<-paste("SELECT OBSINT.DEBRIEFED_LENGTH.YEAR,\n ",
    "OBSINT.DEBRIEFED_LENGTH.HAUL_JOIN\n ",
    "FROM OBSINT.DEBRIEFED_LENGTH\n ",
    "WHERE OBSINT.DEBRIEFED_LENGTH.NMFS_AREA BETWEEN ",reg,"\n",
    "AND OBSINT.DEBRIEFED_LENGTH.SPECIES = ",species,"\n ",
    "ORDER BY OBSINT.DEBRIEFED_LENGTH.YEAR",sep="")

   DSS=sqlQuery(AFSC,test)


  DSS1<-aggregate(list(LENGTH=DSS$HAUL_JOIN),by=list(HAUL_JOIN=DSS$HAUL_JOIN,YEAR=DSS$YEAR),FUN=min)
  DSS2<-aggregate(list(LENGTH=DSS1$HAUL_JOIN),by=list(YEAR=DSS1$YEAR),FUN=length)
  DSS2<-data.frame(DSS2)
  DSS2$SS=100+((DSS2$LENGTH/mean(DSS2$LENGTH)))
  DSS2$SS<-apply(DSS2,1,min)
  DSS2
  }