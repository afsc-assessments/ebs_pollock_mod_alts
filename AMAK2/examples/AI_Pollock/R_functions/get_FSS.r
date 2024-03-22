get_FSS<-function(species=201,sy=1977,area="AI"){
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
    "NORPAC.FOREIGN_LENGTH.HAUL_JOIN\n ",
  "FROM NORPAC.FOREIGN_HAUL\n ",
  "INNER JOIN NORPAC.FOREIGN_LENGTH\n ",
  "ON NORPAC.FOREIGN_HAUL.HAUL_JOIN = NORPAC.FOREIGN_LENGTH.HAUL_JOIN\n ",
  "WHERE NORPAC.FOREIGN_HAUL.YEAR   > ",sy,"\n ",
  "AND NORPAC.FOREIGN_HAUL.GENERIC_AREA BETWEEN ",reg,"\n ",
  "AND NORPAC.FOREIGN_LENGTH.SPECIES = ",species,"\n ",
  "ORDER BY NORPAC.FOREIGN_HAUL.YEAR",sep="")

  FSS=sqlQuery(AFSC,test)
  

  FSS1<-aggregate(list(LENGTH=FSS$HAUL_JOIN),by=list(HAUL_JOIN=FSS$HAUL_JOIN,YEAR=FSS$YEAR),FUN=min)
  FSS2<-aggregate(list(LENGTH=FSS1$HAUL_JOIN),by=list(YEAR=FSS1$YEAR),FUN=length)
  FSS2<-data.frame(FSS2)
  FSS2$SS=100+((FSS2$LENGTH/mean(FSS2$LENGTH)))
  FSS2$SS<-apply(FSS2,1,min)
  FSS2
  }