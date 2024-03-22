get_Dlength<-function(species=201,area="AI"){
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
    "OBSINT.DEBRIEFED_LENGTH.NMFS_AREA,\n ",

    "OBSINT.DEBRIEFED_LENGTH.SPECIES,\n ",
    "OBSINT.DEBRIEFED_LENGTH.SEX,\n ",
    "OBSINT.DEBRIEFED_LENGTH.LENGTH,\n ",
    "OBSINT.DEBRIEFED_LENGTH.FREQUENCY,\n ",
    "OBSINT.DEBRIEFED_LENGTH.LATDD_END,\n ",
    "OBSINT.DEBRIEFED_LENGTH.LONDD_END,\n ",
    "OBSINT.DEBRIEFED_LENGTH.GEAR\n ",
  "FROM OBSINT.DEBRIEFED_LENGTH\n ",
  "WHERE OBSINT.DEBRIEFED_LENGTH.NMFS_AREA BETWEEN ",reg,"\n",
  "AND OBSINT.DEBRIEFED_LENGTH.SPECIES = ",species,"\n ",
  "ORDER BY OBSINT.DEBRIEFED_LENGTH.YEAR",sep="")
  
  Dlength=sqlQuery(AFSC,test)
  Dlength
  }