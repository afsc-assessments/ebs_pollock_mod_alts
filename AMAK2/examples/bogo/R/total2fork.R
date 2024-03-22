
total2fork = function(numbers, lengths, factor, bounds=NULL, bin=NULL) {  
  # Conversion from total to fork length 
  # Make the conversion (for 'factor' bigger than 1)
  # or fork to length ('factor' lower than 1)
  # numbers : number at total length
  # lengths : mark class for numbers
  # factor  : conversion factor from fork to total length
  # bounds  : bounds for fork length output
  # bin     : size of the bin for fork length output
  # complaints to: Ricardo Oliveros-Ramos (roliveros@imarpe.gob.pe)
  
  if(!all(diff(diff(lengths)) == 0)) stop("length bins must be equal")
  xbin = diff(lengths)[1]
  if(is.null(bounds)) bounds = floor(range(lengths)/factor)
  if(length(bounds)!=2) stop("'range' must be length 2")
  if(bounds[2]<bounds[1]) stop("upper range must be greater than lower")
  if(is.null(bin)) bin = xbin
  if(bin<=0) stop("'bin' must be strictly positive")
  
  .total2fork = Vectorize(
    function(x, y, factor, bin) {
      if(length(bin)==1L) bin = rep(bin,2)
      x = c(x-bin[1]/2, x+bin[1]/2)/factor
      y = c(y-bin[2]/2, y+bin[2]/2)
      ylow = min(max(x[1], y[1]),x[2])
      yupp = max(min(x[2], y[2]), x[1])
      ybar = c(ylow, yupp)
      out = abs(diff(ybar))/(bin[1]/factor)
      return(out)
    } , vectorize.args = c("x","y"))
  
  olengths = seq(from=bounds[1], to=bounds[2]+bin, by=bin)  
  out = outer(X=lengths, Y=olengths, FUN=.total2fork, 
              factor=factor, bin=c(xbin, bin))
  out = apply(out*numbers, 2, sum, na.rm=TRUE)
  names(out) = olengths
  return(out)
}

