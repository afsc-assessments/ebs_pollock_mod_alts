require(PBSadmb)
source("adfunctions.r")

mod1spl<- readList("for_r.dat")
mod1smu<- readList("for_r.dat")

mod1b<- readList("arc/mod1b_r.dat")
mod1c<- readList("arc/mod1c_r.dat")
mod2 <- readList("arc/mod2_r.dat")
mod3 <- readList("arc/mod3_r.dat")
mod4 <- readList("arc/mod4_r.dat")
mod5 <- readList("arc/mod5_r.dat")

# show projections
plt_proj(mod1c)

# show fit to catch biomass
CatchFit(mod1c) 

# show spawning biomass relative to population with no fishing
spwn_ratio(mod1c) 

# example of writing multiple plots to pdf file:
pdf("agefits.pdf",width=9, height=7)
  AgeFits(mod1c,f=1)
  AgeFits(mod1c,f=2)
  AgeFits(mod1c,f=3)
  AgeFits(mod1c,f=4)
dev.off()

# another example of writing multiple plots to pdf file:
pdf("indices.pdf",width=9, height=7)
    Indices(mod1a,"Model 1a")
    Indices(mod1b,"Model 1b")
    Indices(mod1c,"Model 1c")
    Indices(mod2,"Model 2")
    Indices(mod3,"Model 3")
    Indices(mod4,"Model 4")
    Indices(mod5,"Model 5")
dev.off()

pdf("selectivity.pdf",width=9, height=7)
Mntns(mod1a,"Model 1a")
Mntns(mod1b,"Model 1b")
Mntns(mod1c,"Model 1c")
Mntns(mod2 ,"Model 2")
Mntns(mod3 ,"Model 3")
Mntns(mod4 ,"Model 4")
Mntns(mod5 ,"Model 5")
dev.off()

Mntns_srv(mod1c,"Model 1c")

# Stock recruitment curve
 p.stock.rec(mod1c)
# recruitment hist w/ errors
 p.rec.hist(mod1c,ylab="Age 2 recruitment",fy=1970,ly=2010)

# Fishing mortality 
 p.full.f(mod1c,f=2)

# sample sizes
 p.eff.n(mod1c,typ="f")
 p.eff.n(mod1c,typ="f",f=2)
 p.eff.n(mod1c,typ="f",f=3)
 p.eff.n(mod1c,typ="f",f=4)

# spawning biomass and last year's estimates 
 p.biom.pol(mod1c,typ="TB")
 lines(mod1a$TotBiom[,1],mod1a$TotBiom[,2],lty=3,lw=3)

# Survey fit
 p.sur.stk(mod1c,f=1)

# Rec
 p.biom.stk(mod1,typ="R")

# Numbers at age
 p.bub.age(mod1c,siz=3000)

