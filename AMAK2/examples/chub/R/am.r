require(PBSadmb)
source("/Users/Jim/Dropbox/r_common/adfunctions.r")
setwd("/Users/Jim/_mymods/AMAK/examples/chub/")

mod1 <- readList("arc/mod1_R.rep")
mod2 <- readList("arc/mod2_R.rep")
mod3 <- readList("arc/mod3_R.rep")
length(mod3$Like_Comp)
rbind(mod2$Like_Comp,mod3$Like_Comp)
# show projections
plt_proj(mod1,fy=1970,ly=2020)
plt_proj(mod2)

# show fit to catch biomass
CatchFit(mod1) 
CatchFit
p.catch.fit(mod1,f=1,ylab="Catch biomass (t)" ,ylim=c(0,1700000))
p.catch.fit(mod1,f=2,ylab="Catch biomass (t)" ,ylim=c(0,1700000))
p.catch.fit(mod1,f=3,ylab="Catch biomass (t)" ,ylim=c(0,1700000))

sel.age.mountain(mod2, f=1, xvec=NULL, yvec=NULL, zscale=3, fy=1977,ly=2014)
,nshades=50, xaxs="i", yaxs="i", xlab="", ylab="", las=1, new=TRUE, addbox=FALSE, cex.xax=1, cex.yax=1, main="", ...)
p.catch.fit(mod2,f=1,ylab="Catch biomass (t)" ,ylim=c(0,500000))
p.catch.fit(mod2,f=2,ylab="Catch biomass (t)" ,ylim=c(0,500000))


# example of writing multiple plots to pdf file:
pdf("agefits.pdf",width=9, height=7)
  AgeFits(mod1,f=1)
  AgeFits(mod2,f=1)
dev.off()

# another example of writing multiple plots to pdf file:
pdf("indices.pdf",width=9, height=7)
    Indices(mod1,"Model 1",fy=1970,ly=2015)
    Indices(mod2,"Model 2",fy=1970,ly=2015)
    Indices(mod3,"Model 3",fy=1970,ly=2015)
dev.off()

pdf("selectivity.pdf",width=9, height=7)
Mntns(mod1,"Model 1")
Mntns(mod2 ,"Model 2")
dev.off()

Mntns_srv(mod1,"Model 1")
detach()
# Stock recruitment curve
styr=1950
 p.stock.rec(mod1)
 p.stock.rec(mod2,main="Model 2")
# recruitment hist w/ errors
 vpa_r <- read.table("clipboard")
 rec_age=0
 p.rec.hist(mod1,ylab="Age 0 recruitment",fy=1950,ly=2014,main="model 1")
 p.rec.hist(mod2,ylab="Age 0 recruitment",fy=1950,ly=2014,main="model 2")
  base=mod1$Stock_Rec[26:64,4]
length(base)  
  wchina=mod2$Stock_Rec[20:65,4]
  vpa_r=(t(vpa_r))
  plot(base,vpa_r,xlab="Base",ylab="VPA",pch=19)
  plot(mod1$Stock_Rec[20:65,4],mod2$Stock_Rec[20:65,4],ylab="Base",xlab="Including Chinese catches",pch=19)
  lines(supsmu(wchina,base))
  ?supsmu


# Fishing mortality 
 p.full.f(mod1)
 p.full.f(mod2)

# sample sizes
 p.eff.n(mod1,typ="F")
 p.eff.n(mod2,typ="F",main="Model 2")

# spawning biomass and last year's estimates 
 p.biom.pol(mod1,main="Model 1")
 p.biom.pol(mod2,main="Model 2",fy=1950,ly=2014)
 p.biom.pol(mod3,main="Model 3")
 lines(mod1a$TotBiom[,1],mod1a$TotBiom[,2],lty=3,lw=3)

# Survey fit
 p.sur.stk(mod1,S=1)
detach()
# Rec
 p.biom.stk(mod1,typ="R")

# Numbers at age
 p.bub.age(mod1c,siz=3000)


# show spawning biomass relative to population with no fishing
spwn_ratio(mod1) 