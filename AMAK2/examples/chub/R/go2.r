setwd("c:/Users/jianelli/Desktop/Chile SPRFMO/prg/Aug")

require(PBSadmb)
source('ADfunctions2.r')
source('ADfunctions.r')
source('IndexFit.r')

mod1<-readList("For_R.dat")

mod2<-readList("For_R.dat")
mod3<-readList("For_R.dat")
mod4<-readList("For_R.dat")
# comparing effective sample sizes for 2 fisheries in 2 models
windows(,9,6)
par(mfrow=c(2,2))
for(i in 1:length(mod1$Fshry_names)) p.eff.n(mod1, typ='F', f=i, new=FALSE)
for(i in 1:length(mod2$Fshry_names)) p.eff.n(mod2, typ='F', f=i, new=FALSE)

# 1 model
p.biom.pol(mod1, n.mod=1, typ='SSB')
p.biom.pol(mod3, n.mod=1, typ='Tot')

# 2 models - Spawning Stock Biomass
p.biom.pol(mod1, mod2, n.mod=2, typ='SSB')
legend('topright', legend=paste(c(2,3), 'Fisheries'), fill=c('light blue', 'light green'), bty='n')

# 2 models - Total Biomass
p.biom.pol(mod1, mod2, n.mod=2, typ='TB')
legend('topright', legend=paste(c(2,3), 'Fisheries'), fill=c('light blue', 'light green'), bty='n')

# construct plot of selectivies over time and age
# Model 1
# Fisheries selectivies
windows(,9,6.5)
par(mfrow=c(1, length(mod1$Fshry_names)))
for(i in 1:length(mod1$Fshry_names))
{
	sel.age.mountain(mod1, f=i, typ='F', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2.5, new=F, cex.yax=0.6)
}
Mntns(mod4)
windows(,9,6.5)
par(mfrow=c(1, length(mod3$Fshry_names)))
for(i in 1:length(mod3$Fshry_names))
{
	sel.age.mountain(mod3, f=i, typ='F', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2.5, new=F, cex.yax=2.)
}
sel.age.mountain(mod4, f=4, typ='S', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2.5, new=F, cex.yax=1)

windows(,9,6.5)
par(mfrow=c(1, length(mod3$Fshry_names)))
for(i in 1:length(mod3$Fshry_names))
{
	sel.age.mountain(mod3, f=i, typ='F', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2.5, new=F, cex.yax=0.6)
}
# Index selectivies
windows(,11, 6.5)
par(mfrow=c(1, length(mod3$Index_names)))
for(i in 1:length(mod3$Index_names))
{
	sel.age.mountain(mod3, f=i, typ='S', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2, new=F, cex.yax=0.6)
}



# Model 2
# Fisheries selectivies
windows(,9,6.5)
par(mfrow=c(1, length(mod2$Fshry_names)))
for(i in 1:length(mod2$Fshry_names))
{
	sel.age.mountain(mod2, f=i, typ='F', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2.5, new=F, cex.yax=0.6)
}

# Index selectivies
windows(,11, 6.5)
par(mfrow=c(1, length(mod2$Index_names)))
for(i in 1:length(mod2$Index_names))
{
	sel.age.mountain(mod2, f=i, typ='S', xvec=c(2:12), xlab='Age', ylab='Year', zscale=2, new=F, cex.yax=0.6)
}
