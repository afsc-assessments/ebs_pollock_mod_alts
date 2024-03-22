GET_CATCH<-function(area="AI",species="'PLCK'",FYR=fyear,ADD_OLD=TRUE,OLD_FILE="OLD_CATCH.csv"){

   ##Define area
 if(area=="AI"){
    reg<-"'539' and '544'"
  }
  if(area=="GOA"){
    reg<-"'600' and '699'"
  }
  if(area=="BS"){
    reg<-"'500' and '539'"
  }
  
   YEARS<-c(1991:FYR)
   CATCH<-vector("list",length=length(YEARS))
   blnd<-c(as.character(91:99),"00","01","02")
   
  for(i in 1:12){
     test=paste("SELECT SUM(BLEND.BLEND",blnd[i],".TONS)AS TONS,\n ",
     "BLEND.BLEND",blnd[i],".ZONE,\n ",
     "BLEND.BLEND",blnd[i],".SPECN,\n ",
     "BLEND.BLEND",blnd[i],".SPEC,\n ",
     "BLEND.BLEND",blnd[i],".TYPE\n ",
     "FROM BLEND.BLEND",blnd[i],"\n ",
     "GROUP BY BLEND.BLEND",blnd[i],".ZONE,\n ",
     "BLEND.BLEND",blnd[i],".SPECN,\n ",
     "BLEND.BLEND",blnd[i],".SPEC,\n ",
      "BLEND.BLEND",blnd[i],".TYPE\n ",
     "HAVING BLEND.BLEND",blnd[i],".ZONE BETWEEN ",reg,"\n ",
     "AND BLEND.BLEND",blnd[i],".SPEC = ",species, sep="")
 
    CATCH[[i]]<-sqlQuery(AFSC,test)
    CATCH[[i]]$YEAR=YEARS[i]
      }


  for(i in 13:(FYR-1990)){

    test=paste("SELECT SUM(BLEND.CAS",YEARS[i],".WEIGHT_POSTED) AS TONS,\n ",
      "BLEND.CAS",YEARS[i],".REPORTING_AREA_CODE AS ZONE,\n ",
      "BLEND.CAS",YEARS[i],".AGENCY_SPECIES_CODE AS SPECN,\n ",
      "BLEND.CAS",YEARS[i],".AGENCY_SPECIES_ID AS SPEC,\n ",
      "BLEND.CAS",YEARS[i],".SOURCE_TABLE AS TYPE,\n ",
      "BLEND.CAS",YEARS[i],".YEAR\n ",
      "FROM BLEND.CAS",YEARS[i],"\n ",
      "GROUP BY BLEND.CAS",YEARS[i],".REPORTING_AREA_CODE,\n ",
      "BLEND.CAS",YEARS[i],".AGENCY_SPECIES_CODE,\n ",
      "BLEND.CAS",YEARS[i],".AGENCY_SPECIES_ID,\n ",
      "BLEND.CAS",YEARS[i],".SPECIES_GROUP_CODE,\n ",
      "BLEND.CAS",YEARS[i],".SOURCE_TABLE,\n ",
      "BLEND.CAS",YEARS[i],".YEAR\n ",
      "HAVING BLEND.CAS",YEARS[i],".REPORTING_AREA_CODE BETWEEN ",reg,"\n ",
      "AND BLEND.CAS",YEARS[i],".SPECIES_GROUP_CODE = ",species,sep="")

    CATCH[[i]]<-sqlQuery(AFSC,test)
    CATCH[[i]]$YEAR=YEARS[i]
    }



  CATCH1<-CATCH[[1]]
  for(i in 2:length(YEARS)){
    CATCH1<-rbind(CATCH1,CATCH[[i]])
  }
  CATCH<-CATCH1

  CATCH_TOTAL<-aggregate(list(TONS=CATCH$TONS),by=list(YEAR=CATCH$YEAR),FUN=sum)
  
  if(ADD_OLD){
    OLD_CATCH<-read.csv(OLD_FILE,header=T)
    CATCH_TOTAL<-rbind(OLD_CATCH,CATCH_TOTAL)
  }
  
  CATCH_TOTAL
}
