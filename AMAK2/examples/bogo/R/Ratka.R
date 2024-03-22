setwd("c:/users/jim/documents/dropbox/models/atka/prg")
source("../Ratka/read.admb.R")
system("dir")
d=read.dat("am2013.dat")
names(d)
setwd("c:/users/jim/documents/dropbox/models/atka/Ratka")
package.skeleton(name="test",code_files="read.admb.R")
require(test,lib.loc=".")
# Still need to fix the ctl file SR window...
R2_admb_input(d,"shit.shit",retro=8,removeall=F)
detach(d)
names(d[[12]])  
d[[12]]
names(d)
d$ind_names