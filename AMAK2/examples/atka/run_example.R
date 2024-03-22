## Demonstrate adnuts on Hake model

library(adnuts)
packageVersion('adnuts')                # needs to be 1.0.1
library(snowfall)
library(rstan)
library(shinystan)
chains <- parallel::detectCores()-1 # chains to run in parallel
## Reproducible seeds are passed to ADMB
set.seed(352)
seeds <- sample(1:1e4, size=chains)
chains

m <- 'amak'
d <- 'ws1'
setwd(d)
file.copy('../../../src/amak.tpl', 'amak.tpl')
system(paste('admb ',m))
system(paste(m, '-nox -iprint 200 -mcmc 10'))
setwd('..')

## Run single RWM chain to make sure it works in serial
test <- sample_admb(model=m, path=d, iter=500, algorithm='RWM', warmup=250,
                    seeds=seeds, chains=1)
test$cmd[1]
## Now with parallel
test <- sample_admb(model=m, path=d, iter=1000, algorithm='RWM', warmup=500,
                    seeds=seeds, chains=chains,
                    parallel=TRUE, cores=chains)
launch_shinyadmb(test)


## Now with more thinning
thin <- 10
iter <- 1000*thin
fit <- sample_admb(model=m, path=d, iter=iter, algorithm='RWM', warmup=iter/4,
                    seeds=seeds, chains=chains, thin=thin,
                    parallel=TRUE, cores=chains)
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
system(paste0(m, ' -binp ',m,'.bar -nox -iprint 200 -mcmc 10 -hbf 1'))
setwd('..')

## Never thin NUTS and use 500-1000 iterations per chain w/
## approximately 1/4 warmup
iter <- 500

## Two different ways to gets NUTS working. First is to use the
## Hessian (metric) just like with the RMW above. Note that
## control argument.
fit.mle <- sample_admb(model=m, path=d, iter=iter, algorithm='NUTS', warmup=iter/4,
                   seeds=seeds, parallel=TRUE, chains=chains,
                   cores=chains, control=list(metric='mle'))
launch_shinyadmb(fit.mle)
mon <- monitor(fit.mle$samples, warmup=fit.mle$warmup, print=FALSE)
ess <- mon[,'n_eff']
(slow <- names(sort(ess))[1:6]) # slowest mixing parameters
pairs_admb(fit.mle, pars=slow)

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
iter <- 100
fit.diag <- sample_admb(model=m, path=d, iter=iter, algorithm='NUTS', warmup=iter/4,
                   seeds=seeds, parallel=TRUE, chains=chains, cores=chains)

## Now the samples from fit.diag can be used to estimate the
## covariance and that can be used directly.
fit.updated <- sample_admb(model=m, path=d, iter=iter, algorithm='NUTS', warmup=iter/4,
                   seeds=seeds, parallel=TRUE, chains=chains,
                   cores=chains, control=list(metric=fit$covar.est))


save.image()
