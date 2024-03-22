read.fit<-function(file){
  # Function to read a basic fit
  ret<-list()
  parfile<-as.numeric(scan(paste(file,'.par', sep=''), 
                      what='', n=16, quiet=TRUE)[c(6,11,16)])
  ret$nopar<-as.integer(parfile[1])
  ret$nlogl<-parfile[2]
  ret$maxgrad<-parfile[3]
  file<-paste(file,'.cor', sep='')
  lin<-readLines(file)
  ret$npar<-length(lin)-2
  ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2])
  sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!=''])
  ret$names<-unlist(lapply(sublin,function(x)x[2]))
  ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3])))
  ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4])))

  ret$cor<-matrix(NA, ret$npar, ret$npar)
  for(i in 1:ret$npar){
    ret$cor[1:i,i]<-as.numeric(unlist(lapply(sublin[i],
      function(x)x[5:(4+i)])))
    ret$cor[i,1:i]<-ret$cor[1:i,i]
  }
  ret$cov<-ret$cor*(ret$std%o%ret$std)
  return(ret)
}

read.outfile<-function(path.string)
{
	tt.chars <- scan(path.string, sep = "", what = "")
	tt.nums <- as.numeric(tt.chars)
	list.length <- length(tt.nums[is.na(tt.nums)])
	out <- vector("list", list.length)
	object.counter <- 1
	element.counter <- 1
	while(object.counter <= list.length) {
		names(out[object.counter]) <- tt.chars[element.counter]
		print(tt.chars[element.counter])
		element.counter <- element.counter + 1
		tt.rows <- tt.nums[element.counter]
		element.counter <- element.counter + 1
		tt.cols <- tt.nums[element.counter]
		element.counter <- element.counter + 1
		if(tt.rows > 1) {
			out[[object.counter]] <- matrix(data = tt.nums[element.counter:(element.counter + tt.rows * tt.cols - 1)], nrow = 
				tt.rows, ncol = tt.cols, byrow = T)
			element.counter <- element.counter + (tt.rows * tt.cols)
		}
		if(tt.rows == 1) {
			out[[object.counter]] <- c(tt.nums[element.counter:(element.counter + tt.cols - 1)])
			element.counter <- element.counter + tt.cols
		}
		object.counter <- object.counter + 1
	}
	out
}

