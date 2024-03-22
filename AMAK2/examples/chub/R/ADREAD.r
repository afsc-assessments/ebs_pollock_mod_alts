adread <- function(filename)
{
#
# For use with AD Model Builder output.
#
# expecting report format (e.g.):
#
# report << "Number: seq(8,17) " << endl;
# report << N << endl;
# report << "Catch:'a','b','c','d','e','f','g','h','i','j' " << endl;
# report << C << endl;
#
# Note that ":" will be treated as a special delimiter.
# 

xxx = system(paste("grep -n : ",filename,sep=""),intern=T)
xs=strsplit(xxx,":")
ll = unlist(lapply(xs,length))  #note that some have dimnames
id=c(1,1+cumsum(ll[-length(ll)]))
xsu=unlist(xs)
print(xsu)
the.breaks<-as.numeric(xsu[id])
the.names<- xsu[id+1]
the.cnames<-as.character(rep("",length(ll)))
the.cnames[ll>2]= xsu[id[ll>2]+2]

the.count<-length(the.names)
if(the.count!=length(unique(the.names))) stop("The names are not unique.")
if(the.count==1) 
   { 
   the.object<-read.table(filename,skip=1,
      col.names=eval(parse(text=paste("c(",the.cnames,")"))))
   }
if(the.count >1) 
   {
#
   the.nrecords<-diff(the.breaks)-1
#
   the.object<-list(100)
   for(i in 1:(the.count-1))
   {
      shell(paste("sed 1,",the.breaks[i],"d ",filename," > file1.tmp",sep=""))
      shell(paste("head -",the.nrecords[i]," file1.tmp > file2.tmp",sep=""))
      if (the.nrecords[i] == 1)
          the.object[[i]]<-scan("file2.tmp")
      else
      {
      xx<-read.table("file2.tmp")
       if (ll[i]>2) colnames(xx) = eval(parse(text=paste("c(",the.cnames[i],")")))
       the.object[[i]] <- as.matrix(xx)
      }
   }
   shell(paste("sed 1,",the.breaks[the.count],"d ",filename," > file1.tmp",sep=""))
#
   xx <- read.table("file1.tmp")
   if (ll[the.count] > 2) colnames(xx)=eval(parse(text=paste("c(",the.cnames[the.count],")")))
   xx  <-as.matrix(xx)
   the.object[[the.count]]= xx
   shell("del file1.tmp")
   shell("del file2.tmp")
#
   } # end if the.count > 1
#
names(the.object)<-the.names
#
return(the.object)
}



