### This file demonstrates how to run Bayesian inference on the pollock
### assessment using the adnuts R package.

radian
## as well. This folder gets copied during parallel runs.


## Demonstrate adnuts on Hake model

library(adnuts)
library(matrixcalc)
packageVersion('adnuts')                # needs to be 1.0.1
library(snowfall)
library(rstan)
library(shinystan)
chains <- parallel::detectCores()-2 # chains to run in parallel
## Reproducible seeds are passed to ADMB
set.seed(352)
seeds <- sample(1:1e4, size=chains)
chains
d<-'.'
m <- 'test' # the model name, folder is also assumed to be called this runs/base/mcmc
sample_admb(model='t1', path='.', iter=500, algorithm='RWM', warmup=250, seeds=seeds, chains=1) 
system('test -nox')
sample_admb(model='t1', path='.', iter=500, algorithm='RWM', warmup=250, seeds=seeds, chains=1,extra.args="-binp t1.bar")
d<-'mod'
m <- 'jjms' # the model name, folder is also assumed to be called this runs/base/mcmc
setwd(d)
getwd()
system(paste(m, '-nox -iprint 200 -mcmc 10'))
setwd('..')
sample_admb(model='t', path='.', iter=500, algorithm='RWM', warmup=250, seeds=seeds, chains=1)

## Run single RWM chain to make sure it works in serial
test <- sample_admb(model=m, path=d, iter=500, algorithm='RWM', warmup=250,
                    seeds=seeds, chains=1)
test$cmd[1]
## Now with parallel
test <- sample_admb(model=m, path=d, iter=1000, algorithm='RWM', warmup=500, seeds=seeds, chains=chains, parallel=TRUE, cores=chains,extra.args="-binp jjms.bar")
launch_shinyadmb(test)
## Diagnose pilot analsyis
mon <- monitor(test$samples, warmup=test$warmup, print=FALSE)
ess <- mon[,'n_eff']
(slow <- names(sort(ess))[1:6]) # slowest mixing parameters
pairs_admb(fit=test, pars=slow)
(fast <- names(sort(ess, decr=TRUE))[1:6]) # fastest
pairs_admb(fit=test, pars=fast)


## Now with more thinning
thin <- 200
iter <- 1000*thin
mrw <- 'jjms'
#test <- sample_admb(model=m, path=d, iter=1000, algorithm='RWM', warmup=500, seeds=seeds, chains=chains, parallel=TRUE, cores=chains,extra.args="-binp jjms.bar")
  fit <- sample_admb(model=mrw, path=d, iter=iter, algorithm='RWM', warmup=iter*.2, seeds=seeds, chains=chains, thin=thin, parallel=TRUE, cores=chains, extra.args="-binp jjms.bar")
launch_shinyadmb(fit)

## Diagnose pilot analsyis
mon <- monitor(fit$samples, warmup=fit$warmup, print=FALSE)
ess <- mon[,'n_eff']
(slow <- names(sort(ess))[1:6]) # slowest mixing parameters
pairs_admb(fit=fit, pars=slow)
(fast <- names(sort(ess, decr=TRUE))[1:6]) # fastest
pairs_admb(fit=fit, pars=fast)

## Another trick is you can run for a certain duration if you
## have a time limit. Specify the warmup period and then a really
## large value for iter and set a duration


### End of pilot exploration.


### Explore using NUTS
setwd(d)
getwd()
system(paste0(m, ' -binp ',m,'.bar -phase 22 -nox -iprint 200 -mcmc 10 -hbf 1'))
setwd('..')

## Never thin NUTS and use 500-1000 iterations per chain w/
## approximately 1/4 warmup
iter <- 1000

## Two different ways to gets NUTS working. First is to use the
## Hessian (metric) just like with the RMW above. Note that
## control argument.
fit.mle <- sample_admb(model=m, path=d, iter=iter, algorithm='NUTS', warmup=iter/4, seeds=seeds, parallel=TRUE, chains=chains, cores=chains, control=list(metric='mle'),extra.args='-binp jjms.bar')
launch_shinyadmb(fit.mle)
mon <- monitor(fit.mle$samples, warmup=fit.mle$warmup, print=FALSE)
ess <- mon[,'n_eff']
(slow <- names(sort(ess))[1:6]) # slowest mixing parameters
pairs_admb(fit.mle, pars=slow)
(fast <- names(sort(ess, decr=TRUE))[1:6]) # fastest
pairs_admb(fit=fit.mle, pars=fast)


## Look at high correlations
library(corrplot)
x <- fit.mle$mle$cor
dimnames(x) <- list('par'=fit.mle$mle$par.names, 'par2'=fit.mle$mle$par.names)
ind <- sort(unique(which(abs(x)>.6 & x!= 1, arr.ind=TRUE)[,1]))
y <- x[ind, ind]
corrplot(y, method='color', type='upper')


## Alternatively if no Hessian is available (e.g., b/c of
## hierarchical model), then adapt a diagonal one during
## warmup. This is much slower b/c it doesn't know the shape of
## the posterior

fit.diag <- sample_admb(model=m, path=d, iter=iter, algorithm='NUTS', warmup=iter/4,
                   seeds=seeds, parallel=TRUE, chains=chains, cores=chains)

## Now the samples from fit.diag can be used to estimate the

## covariance and that can be used directly.
iter <- 3000
fit.updated <- sample_admb(model=m, path=d, iter=iter, algorithm='NUTS', thin=3, warmup=iter/4,
                   seeds=seeds, parallel=TRUE, chains=chains,
                   cores=chains, control=list(metric=fit.mle$covar.est),extra.args='-binp jjms.bar')
launch_shinyadmb(fit.updated)
mon <- monitor(fit.updated$samples, warmup=fit.updated$warmup, print=FALSE)
ess <- mon[,'n_eff']
(slow <- names(sort(ess))[1:6]) # slowest mixing parameters
pairs_admb(fit.updated, pars=slow)
(fast <- names(sort(ess, decr=TRUE))[1:6]) # fastest
pairs_admb(fit=fit.updated, pars=fast)


library(ggthemes)
library(tidyverse)
library(scales)
library(PBSmodelling)
library(ggplot2)
M <- readList("For_R.rep")
names(M)
M
names(ssb.mle) <- c("Year","SSB", "SD","LB","UB")
ssb.mle <- ssb.mle %>% filter(Year > 1990, Year<2020)
c("Year","SSB", "SD","LB","UB")
mc.df <- data.frame(read.table(paste0("mceval.dat"),header=TRUE,as.is=TRUE))
mc1.df  <- mc.df %>%filter(type=="sel",Year >= 1991, Year<2000) %>% transmute(Sam, Year=as.factor(Year), Age=as.numeric(Age), Selectivity=value)
tail(mc1.df)
ggplot(mc1.df,aes(x=Age,y=Selectivity,fill=Year)) + expand_limits(y = 0) +
        ggplot2::stat_summary(fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975), geom = "ribbon", alpha = 0.25, colour = NA) +
        ggplot2::stat_summary(fun.y = function(x) quantile(x, 0.5), geom = "line", lwd = .5) +
        theme_few()+ labs(x="Year",y='Partial F')  + facet_wrap(.~Year) #+ geom_ribbon(data=ssb.mle,aes(x=Year, ymax=UB, ymin=LB), alpha=.2,fill="blue") + geom_line(data=ssb.mle,aes(x=Year, y=SSB), color="blue")
mc1.df  <- mc.df %>%filter(type=="sel",Year >= 1990, Year<=2020) %>% transmute(Sam, Year=as.factor(Year), Age=as.numeric(Age), Selectivity=value)
mc1.df %>% filter(between(Sam,20,40)) %>% mutate(Sam=as.factor(Sam)) %>% arrange((Age)) %>% ggplot(aes(x=Age,y=Selectivity,col=Year,group=Sam)) + expand_limits(y = 0) + geom_line(alpha=.5) + 
        theme_few()+ labs(x="Year",y='Partial F')  + facet_wrap(.~Year) #+ geom_ribbon(data=ssb.mle,aes(x=Year, ymax=UB, ymin=LB), alpha=.2,fill="blue") + geom_line(data=ssb.mle,aes(x=Year, y=SSB), color="blue")

mc1.df %>% filter(Sam==2|Sam==3,Year==1997)
head(mc1.df)
head(mc1.df)
mc2.df  <- mc.df %>%filter(type=="SSB",Year > 1990, Year<2020) %>% select(Sam, Year, SSB=value)
names(mc2.df)
mc2.df

ggplot(mc2.df,aes(x=Year,y=SSB)) + 
        expand_limits(y = 0) +
        ggplot2::stat_summary(fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975), geom = "ribbon", alpha = 0.25, colour = NA) +
        ggplot2::stat_summary(fun.y = function(x) quantile(x, 0.5), geom = "line", lwd = .5) +
        theme_few()+ labs(x="Year",y='SSB (t)')    + geom_ribbon(data=ssb.mle,aes(x=Year, ymax=UB, ymin=LB), alpha=.2,fill="blue") + geom_line(data=ssb.mle,aes(x=Year, y=SSB), color="blue")

ggplot(mc.df, aes(x=Year,y=SSB)) + geom_point() + theme_few() 

save.image()
#-------------------------------------------------------------------------------
#Selectivity
library(ggridges)
m0<- readList('For_R.rep')
df <- data.frame(m0$sel_fsh_1[,2:23] ); names(df) <- c("yr",1:21)
sdf <- gather(df,age,sel,2:22) %>% filter(yr>1976) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .9,color="blue",fill="yellow",size=.5)+ xlim(c(1,21))+ ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr)))) + theme_few()

ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .4,color="blue",fill="yellow",size=.5)+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr))))
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .4,color="black",fill="orange",size=.5)+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr))))
ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .2,fill="orange")+ xlim(c(1,11))+ mytheme + ylab("Year") + xlab("Age (years)")+ xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr)))) 

 ssb.mle <- data.frame(m0$SSB)
ssb.mle <- ssb.mle %>% filter(Year > 1990, Year<2020)

 ssb.parsel <- data.frame(m2$SSB,model="Parametric")
 ssb.mle    <- data.frame(m0$SSB,model="non-Parametric")
ssb <- rbind(ssb.mle,ssb.parsel)
names(ssb) <- c("Year","SSB", "SD","LB","UB","Model")
ggplot(ssb,aes(x=Year,y=SSB, ymin=LB,ymax=UB,fill=Model)) + geom_ribbon( alpha = .3) + xlab("Year") + theme_few()


lstOut1  <- list( "base"= m0, "param1"=m1,"parametric"= m2)
# Likelihood table
lstOuts <- lstOut1
tab       <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))
do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)}))
do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)}))
tab  <- rbind(
cbind("q",do.call(cbind,lapply(lstOuts,function(x){round(x[["q_1"]][1,2],2)}))),
cbind("Npars",do.call(cbind,lapply(lstOuts,function(x){round(x[["Num_parameters_Est"]],2)}))),
cbind("M",do.call(cbind,lapply(lstOuts,function(x){round(x[["Mest"]][2],2)}))),
cbind("SigmaR",do.call(cbind,lapply(lstOuts,function(x){round(x[["Sigmar"]][2],2)}))),
cbind("EffN_Fish",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["EffN_Fsh_1"]][,2]),2)}))),
cbind("EffN_Surv",do.call(cbind,lapply(lstOuts,function(x){round(mean(x[["EffN_Survey_1"]][,2]),2)}))),
tab,
cbind("F2016",do.call(cbind,lapply(lstOuts,function(x){round((x[["F_fsh_1"]][40,2]),2)}))),
cbind("F2016/F40%",do.call(cbind,lapply(lstOuts,function(x){round((x[["F_fsh_1"]][40,2])/(x[["F40_est"]]),2)}))),
cbind("B 1977",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][1,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][1,3])/(x[["TotBiom"]][1,2])*100,0)}))),
cbind("B 2016",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][40,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["TotBiom"]][40,3])/(x[["TotBiom"]][38,2])*100,0)}))),
cbind("2006 YC",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][31,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][31,3])/(x[["R"]][31,2])*100,0)}))),
cbind("2012 YC",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][37,2]),0)}))),
cbind("CV",do.call(cbind,lapply(lstOuts,function(x){round((x[["R"]][37,3])/(x[["R"]][37,2])*100,0)})))
)
tab
i=2
names(lstOuts)
for (i in 1:2){
  df <- data.table( lstOuts[[i]]$Obs_Survey_1[-1,] )
  names(df) <- c("Year","obs","pre","sd","pr","lnpr")
  df$pre <- as.numeric(df$pre)
  df$obs <- as.numeric(df$obs)
  df <- df[df$obs!="NA",]
  print( df %>% summarize(sqrt(mean(log(obs/pre)^2))) )
}
plt_srr<-function(M,xlab="Spawning biomass",ylab="Recruits",main="Model 16.0b") {
  df <- data.frame( yrs  = M$Stock_Rec[,1] , ssb  = M$Stock_Rec[,2], rec  = M$Stock_Rec[,4])
  df2 <- data.frame(stock = M$stock_Rec_Curve[,1], rec = M$stock_Rec_Curve[,2])
  ggplot() + geom_path(data=df,aes(x=ssb,y=rec) ) +geom_text(data=df,aes(x=ssb,y=rec, label=yrs)) + 
        scale_y_continuous(labels = comma) +
        scale_x_continuous(labels = comma) +
  xlab(xlab) + ylab(ylab) + ggtitle(main)  + theme_few(base_size=16) + geom_line(data=df2,aes(x=stock, y=rec),col="salmon",size=2) +
  theme(legend.position="none" ) + theme(axis.text.x = element_text(angle = 0))
}
plt_srr(m0,main="Base 2019")
plt_srr(m2,main="Parametric selectivity")
AgeFitsSrv(mod16.0b,rec_age=1,case_label="2018 assessment")

#Rec ============================================================
rdf <- cbind(data.frame(m0$R,case="Base"))
rdf <- rbind(rdf,cbind(data.frame(m2$R,case="Parametric")))
names(rdf) <- c("yr","R","se","lb","ub","case")
rdf  <- rdf %>% filter(yr>1970,yr<2020)
mnR <- mean(m0$R)
dodge <- position_dodge(width=0.8)
ggplot(rdf,aes(x=yr-1,y=R,fill=case)) + xlab("Year class") + ylab("Age 1 recruits (thousands)") + 
       geom_bar(width=0.75,position="dodge",stat="identity",color="black") + 
       scale_x_continuous(breaks=seq(1970,2019,5)) +
       geom_errorbar(aes(ymin=lb,ymax=ub),width=.3,colour="blue",position=dodge) + theme_few() + geom_hline(aes(yintercept=mnR))
source("../plot_ind.R")
plot_ind(M)