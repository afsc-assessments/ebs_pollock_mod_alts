get_Flength<-function(species=201,sy=1977,area="AI"){
   if(area=="AI"){
    reg<-"539 and 544"
  }
  if(area=="GOA"){
    reg<-"600 and 699"
  }
  if(area=="BS"){
    reg<-"500 and 539"
  }
  
  test<-paste("SELECT NORPAC.FOREIGN_HAUL.YEAR,\n ",
    "NORPAC.FOREIGN_HAUL.GENERIC_AREA,\n ",
    "NORPAC.FOREIGN_HAUL.HOOKS_PER_SKATE,\n ",
    "NORPAC.FOREIGN_HAUL.NUMBER_OF_POTS,\n ",
    "NORPAC.FOREIGN_LENGTH.SPECIES,\n ",
    "NORPAC.FOREIGN_LENGTH.SEX,\n ",
    "NORPAC.FOREIGN_LENGTH.SIZE_GROUP,\n ",
    "NORPAC.FOREIGN_LENGTH.FREQUENCY,\n ",
    "NORPAC.FOREIGN_HAUL.LATITUDE,\n ",
    "NORPAC.FOREIGN_HAUL.E_W,\n ",
    "NORPAC.FOREIGN_HAUL.LONGITUDE \n ",
  "FROM NORPAC.FOREIGN_HAUL\n ",
  "INNER JOIN NORPAC.FOREIGN_LENGTH\n ",
  "ON NORPAC.FOREIGN_HAUL.HAUL_JOIN = NORPAC.FOREIGN_LENGTH.HAUL_JOIN\n ",
  "WHERE NORPAC.FOREIGN_HAUL.YEAR   > ",sy,"\n ",
  "AND NORPAC.FOREIGN_HAUL.GENERIC_AREA BETWEEN ",reg,"\n ",
  "AND NORPAC.FOREIGN_LENGTH.SPECIES = ",species,"\n ",
  "ORDER BY NORPAC.FOREIGN_HAUL.YEAR",sep="")
  
  Flength=sqlQuery(AFSC,test)
    Flength$GEAR<-1
  Flength$GEAR[Flength$HOOKS_PER_SKATE > 0]  <- 6
  Flength$GEAR[Flength$NUMBER_OF_POTS > 0]<-8
  Flength<-subset(Flength,select=-c(HOOKS_PER_SKATE,NUMBER_OF_POTS))
  
  
  Flength
  }