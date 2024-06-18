library(RTMB)
library(dplyr) 
library(tidyr)
library(ggplot2)

load("data/am2022.RData")
load("data/sizeage_matrix.RData")
input$sizeage <- sizeage; rm(sizeage)
source("R/helper.R")

head(input$obsdf, 5) # long format with all observations
# obs_type # 0=catch, 1=index, 2=agecom, 3=lencomp
# nll_type # 0=dnorm, 1=dmultinom
# fit_data # 1/0=TRUE/FALSE
# fleet    # 1=fishery, 2=survey
# obs      # transformed appropriately for nll_type (becomes keep vec)
# obserror # if nll_type obs error is an input (note this is Neff for dmultinom)

#Remove length comps for now till we implement them
input$obsdf <- input$obsdf[input$obsdf$obs_type!=3,]

# data list ----
dat <- list()
dat$obs <- input$obsdf$obs
dat$aux <- input$obsdf
dat$aux <- get_id(dat$aux)
dat$aux <- get_likelihood_index(dat$aux)
dat$year <- input$years
dat$minYear <- min(dat$year)
dat$age <-  input$ages
dat$len <-  input$lens
dat$minAge <- min(dat$age)
dat$sampleTimes <- input$srv_frac
dat$spawnTimes <- input$sp_frac
dat$waa <- input$waa
dat$mature <- input$maturity
dat$sizeage <- input$sizeage
dat$fleetTypes <- unique(input$obsdf$fleet) 
dat$srmode <- 0 #
dat$logN_mode <- 0 # 0 = deterministic SCAA, 1 = sigR estimated, logR estimated as ranef, 2 = full state-space with/ shared sigN for a > 1 (same as n_NAA_sigma in WHAM)

# prediction data frame
dat$aux <- get_pred(dat$aux, input)

# parameter ----
par <- list()
par$logsigR <- log(input$sigr)
par$logsigN <- if(dat$logN_mode==2){log(0.5)}else{numeric(0)}
par$logQ <- 0
# is M a constant in FIMS or by year/age?
par$logM <- matrix(log(input$natmort), nrow=length(dat$year), ncol=length(dat$age))
par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)}
par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)}
par$logN <- matrix(10, nrow=length(dat$year), ncol=length(dat$age)) # Tim suggested initializing at 10 rather than zero
par$logFmort <- matrix(0, nrow=length(dat$year), ncol=1)
par$logfshslx <- log(input$fsh_slx) # need parametric selectivity
par$logsrvslx <- log(input$srv_slx)

# assumes vectors at age are supplied to function for N, F, M, W, and mature
calc_ssb <- function(Naa, Faa, M, waa, mature, spawnTimes){
  sum(Naa*exp((-Faa-M)*spawnTimes)*mature*waa)/1e3
}

# model ----
f<-function(par){ # note dat isn't an argument in the fxn
  getAll(par, dat) # RTMB's attach
  obs <- OBS(obs) # access simulation, OSA residuals
  predObs <- rep(0, nrow(aux))
  nobs <- length(obs) 
  nyear <- length(year)
  nage <- length(age)
  
  sigR <- exp(logsigR)
  sigN <- exp(logsigN)
  M <- exp(logM)
  Q <- exp(logQ)
  Fmort <- exp(logFmort)
  fshslx <- exp(logfshslx)
  srvslx <- exp(logsrvslx)
  Faa <- matrix(data = 1, nrow = nyear, ncol = nage) 
  for(y in 1:nyear) Faa[y,] = Fmort[y,] * fshslx 
  Z <- Faa+M
  
  jnll <- 0

  # Recruitment ----
  ssb <- rep(0, nyear)
  for(y in 1:nyear) ssb[y] <- calc_ssb(exp(logN[y,]),Faa[y,],M[y,],waa,mature,spawnTimes)
  predlogR <- rep(0, nyear)
  for(y in 1:nyear){
    thisSSB <- ifelse((y-minAge-1)>(-.5),ssb[y-minAge],ssb[1]) 
    if(srmode==0){ # RW
      if (y == 1){
        predlogR[y] <- logN[y,1] # need to fix this later
      }else{
        predlogR[y] <- logN[y-1,1]
      }
    }
    if(srmode==1){ # Ricker
      predlogR[y] <- rickerpar[1]+log(thisSSB)-exp(rickerpar[2])*thisSSB
    }
    if(srmode==2){ # BH
      predlogR[y] <- bhpar[1]+log(thisSSB)-log(1.0+exp(bhpar[2])*thisSSB)
    }
    if(!(srmode %in% c(0, 1, 2))){
      stop(paste("srmode", srmode, "not implemented yet"))
    }      
    jnll <- jnll - dnorm(logN[y,1],predlogR[y],sigR,log=TRUE)
  }  
  
  # N matrix
  predlogN <- matrix(0, nrow=nyear, ncol=nage)
  for(y in 2:nyear){
    for(a in 2:nage){
      # full state-space model 
      predlogN[y,a] <- logN[y-1,a-1]-Faa[y-1,a-1]-M[y-1,a-1]
      if(a==nage){
        predlogN[y,a] <- log(exp(predlogN[y,a])+exp(logN[y-1,a]-Faa[y-1,a]-M[y-1,a]))
      }
      if(logN_mode == 0){ # fixed or random effects recruitment with deterministic N matrix
        logN[y,a] <- predlogN[y,a]
      }else{
        jnll <- jnll - dnorm(logN[y,a],predlogN[y,a],sigN,log=TRUE)
      }
    }
  }

  # predicted catch ----
  
  # get_pred_logCaa(idx, # lkup vector linking to aux, aux with appropriate flt opts
  #              Z, logN, Faa) # Z, logN, Faa are also vectors
  # need to modify this to allow for multiple fleets
  logpredcatchatage <- logN-log(Z)+log(1-exp(-Z))+log(Faa)

  # obs_type 0 is aggregate catch in weight (need to figure out how we want to input units)
  for (i in which(aux$obs_type == 0)){
    y <- which(year == aux$year[i])
    predObs[i] <- log(sum(exp(logpredcatchatage[y,]) * waa)/1e6) # waa in g and aggregate catch in t
  }
  
  # predicted survey biomass ----
  logpredindexatage <- logQ + logsrvslx + logN - Z * sampleTimes

  # obs_type 1 is survey biomass in weight (need to figure out units)
  for (i in which(aux$obs_type == 1)){
    y <- which(year == aux$year[i])
    predObs[i] <- log(sum(exp(logpredindexatage[y,]) * waa)/1e6) # waa in g and aggregate srv biom in t
  }
  
  # age comps (age error not included) ----
  
  # fishery
  tmp <- exp(logpredcatchatage)
  tmptot <- rowSums(tmp)
  tmp <- tmp/tmptot

  # survey
  tmp2 <- exp(logpredindexatage)
  tmptot2 <- rowSums(tmp2)
  tmp2 <- tmp2/tmptot2

  # combine and vectorize
  tmp3 <- rbind(tmp,tmp2)
  out <- tmp3[1,]
  
  # wow!
  for(i in seq_along(tmp3[,1])[-1]) out <- c(out, tmp3[i,])
  predObs[which(aux$obs_type == 2)] <- out 
  
  # length comps ----
  
  tmp5 <- exp(logpredcatchatage)
  tmptot5 <- rowSums(tmp5)
  tmp5 <- tmp5/tmptot5

  predcatchatlength <- tmp5[,rep(1,length(sizeage[1,]))]
  
  # This is a matrix version of the elementwise loop below
  #
  # predcatchatlength2 <- tmp5[,rep(1,length(sizeage[1,]))]
  # for(i in seq_along(year)){
  #   predcatchatlength2[i,]<- tmp5[i,] %*% sizeage 
  # }
  
  for(i in seq_along(year)){
    for(j in seq_along(len)){
      predcatchatlength[i,j] <- sum(sizeage[,j]*tmp5[i,])
    }
  }
  
  tmp6 <- exp(logpredcatchatage)
  tmptot6 <- rowSums(tmp6)
  tmp6 <- tmp6/tmptot6

  predcatchatlength2 <- tmp6[,rep(1,length(sizeage[1,]))]
  
  for(i in seq_along(year)){
    for(j in seq_along(len)){
      predcatchatlength2[i,j] <- sum(sizeage[,j]*tmp6[i,])
    }
  }
  
  # combine and vectorize
  tmp7 <- rbind(predcatchatlength,predcatchatlength2)
  out2 <- tmp7[1,]
  
  # wow!
  for(i in seq_along(tmp7[,1])[-1]) out2 <- c(out2, tmp7[i,])
  predObs[which(aux$obs_type == 3)] <- out2 
  
  # observational likelihoods ----
  
  for (i in unique(aux$likelihood_index[!is.na(aux$likelihood_index)])){
    
    tmp <- aux[which(aux$likelihood_index==i), ] 
    tmppred <- predObs[which(aux$likelihood_index==i)]
    
    unique_nll_type <- unique(tmp$nll_type)
    if(length(unique_nll_type)>1) stop("multiple nll types within tmp")
    
    # dnorm for catches, indices
    if(unique_nll_type==0) {
      #browser()
      jnll <- jnll - RTMB::dnorm(tmp$obs, tmppred, tmp$obserror, log=TRUE)
    }
    # multinomial for comps
    if(unique_nll_type==1) {
      jnll <- jnll - RTMB::dmultinom(x=tmp$obserror * tmp$obs, size=sum(tmp$obserror * tmp$obs), prob=tmppred, log=TRUE)
    }
  }

  REPORT(predObs)
  
  ADREPORT(predlogR)
  logssb<-log(ssb)
  ADREPORT(logssb)
  jnll
}    

fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}
map <- list()
map$logsigR <- if(dat$logN_mode==0){fill_vals(par$logsigR, NA)}else{factor(1)}
# map$logQ <- fill_vals(par$logQ, NA)
map$logM <- fill_vals(par$logM, NA)
map$logfshslx <- fill_vals(par$logfshslx, NA)
map$logsrvslx <- fill_vals(par$logsrvslx, NA)

nyr <- length(dat$year)
nage <- length(dat$age)
tmp <- matrix(data = NA, ncol = nage, nrow = nyr)
tmp[,1] <- 1:nyr
tmp[1,2:nage] <- (nyr+1):(nyr+nage-1)
map$logN <- as.factor(as.vector(tmp))

obj <- MakeADFun(f, par, 
                 map=map,
                 random=NULL,
                 silent=FALSE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
opt$objective

names(opt$par)
opt$par

sdr <- sdreport(obj)
sdr
plr <- as.list(sdr,report=TRUE, "Est")
plrsd <- as.list(sdr,report=TRUE, "Std")
load('data/orig_am2022.Rdata')

Rec <- as.data.frame(arep$R)
names(Rec) <- c('year', 'R', 'Rec_sd', 'Rec_lci', 'Rec_uci')
Rec <- Rec %>%  mutate(version = 'amak') %>% 
  bind_rows(data.frame(year = dat$year,
               R = exp(plr$predlogR)/1e6,
               Rec_uci = exp(plr$predlogR+2*plrsd$predlogR),
               Rec_lci = exp(plr$predlogR-2*plrsd$predlogR),
               version = 'babyFIMS'))
#View(Rec)
ggplot(Rec %>% filter(year >= 1978), aes(year, (R), col = version)) +
  geom_point() +
  geom_line() + ylim(0,NA) + ggthemes::theme_few()


ssb <- as.data.frame(arep$SSB)
names(ssb) <- c('year', 'ssb', 'ssb_sd', 'ssb_lci', 'ssb_uci')
ssb <- ssb %>% 
  mutate(version = 'amak') %>% 
  bind_rows(data.frame(year = dat$year,
               ssb = exp(plr$logssb),
               ssb_uci = exp(plr$logssb+2*plrsd$logssb),
               ssb_lci = exp(plr$logssb-2*plrsd$logssb),
               version = 'babyFIMS'))

ggplot(ssb %>% filter(year >= 1978), aes(year, (ssb), col = version)) +
  geom_point() +
  geom_line() + ylim(0,NA) + ggthemes::theme_few()


