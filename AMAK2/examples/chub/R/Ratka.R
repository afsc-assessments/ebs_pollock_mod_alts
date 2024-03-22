setwd("/users/jim/_mymods/AMAK/examples/chub")
source("R/read.admb.R")
getwd()
d=read.dat("chub.dat")
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
require(PBSadmb)
source("R/ADfunctions2.r")
source("~/dropbox/R_common/ADfunctions.r")
 mod1<- readList("arc/mod1_R.rep")
 mod2<- readList("arc/mod2_R.rep")
 mod3<- readList("arc/mod3_R.rep")

# show projections
x11()
plt_proj(mod2)

# show fit to catch biomass
CatchFit(mod2) 

# show spawning biomass relative to population with no fishing
spwn_ratio(mod1) 

# example of writing multiple plots to pdf file:
pdf("agefits.pdf",width=9, height=7)
  AgeFits(mod1,f=1,main="All catches")
dev.off()

# another example of writing multiple plots to pdf file:
x11()
pdf("indices.pdf",width=9, height=7)
    Indices(mod1,"Model 1",fy=1970,ly=2014)
    Indices(mod2,"Model 2",fy=1970,ly=2014)
    Indices(mod3,"Model 3",fy=1970,ly=2014)
detach()
dev.off()

pdf("selectivity.pdf",width=9, height=7)
Mntns(mod1,"Model 1")
dev.off()

Mntns_srv(mod1,"Model 1")

# Stock recruitment curve
 p.stock.rec(mod1)
# recruitment hist w/ errors
 rec_age=1
 p.rec.hist(mod1,ylab="Age 1 recruitment",fy=1950,ly=2014,main="Model 1")
 p.rec.hist(mod2,ylab="Age 1 recruitment",fy=1950,ly=2014,main="Model 2")
 p.rec.hist(mod3,ylab="Age 1 recruitment",fy=1950,ly=2014,main="Model 3")

# Fishing mortality 
 p.full.f(mod1,f=1)

# sample sizes
 p.eff.n(mod1,typ="f")

# spawning biomass and last year's estimates 
 p.biom.pol(mod1)
 p.biom.pol(mod2)
 p.biom.pol(mod3)
 detach()

 lines(mod1a$TotBiom[,1],mod1a$TotBiom[,2],lty=3,lw=3)

# Survey fit
 p.sur.stk(mod1c,f=1)

# Rec
 p.biom.stk(mod1,typ="R")

# Numbers at age
 p.bub.age(mod1c,siz=3000)

