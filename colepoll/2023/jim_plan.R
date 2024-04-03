## A quick script to reproduce the final 2023 assessment for GOA pollock

## devtools::install_github('afsc-assessments/GOApollock', ref='v0.1.2')
library(GOApollock)
library(TMB)

## run in ADMB
system("admb gp.tpl")
system("gp")
arep <- read_pk_rep(model='gp', version='EBS test admb', styr=1964, endyr=2023)
asdrep <- read_pk_std(model='gp', version='EBS test admb', styr=1964, endyr=2023)
## all sdreported variables
unique(asdrep$name)
ssb <- filter(asdrep, name=='Espawnbio')
ggplot(ssb, aes(year, est, ymin=lwr, ymax=upr)) + geom_line() +
  geom_ribbon(alpha=.5) + labs(y='SSB (Mt)')
rec <- filter(asdrep, name=='recruit')
ggplot(rec, aes(year, est, ymin=lwr, ymax=upr)) + geom_pointrange() +
  labs(y='Recruits (billions)')

#--Need to bring all model objects into a single list and name them tidy way-----------------
# e.g., dataframe with "Model, TYPE (N, SSB,likelihood), Age, Value"
#

arep$Total_catch
plot(1:10, arep$Fishery_selectivity[50,], type='b')
matplot(arep$Fishery_selectivity)
plot(1:10, arep$Survey_1_selectivity, type='b')
plot(1:10, arep$Survey_2_selectivity, type='b')

filter(asdrep, grepl('q1|q2|q3', name))

## can parse into R and check stuff as needed
dat <- read_pk_dat('gp.dat')
matplot(arep$year,dat$wt_fsh)
rowSums(dat$srvp1)
rowSums(dat$srvp2)
rowSums(dat$srvp3)

## Try in TMB
input <- prepare_pk_input(path=getwd(), datfile='gp.dat', version='EBS test')
input$dat$catp <- input$dat$catp/rowSums(input$dat$catp)
str(input$dat)
str(input$par)

fit <- fit_pk(input, getsd=1, newtonsteps=0)
fit$rep$years <- fit$rep$styr: fit$rep$endyr

str(fit$opt)
str(fit$rep)
str(fit$sd)

round(fit$rep$loglik,2)

tsdrep <- fit$sd
fit$sdrep

x <- get_rep(fit, 'Espawnbio')
names(fit)
ggplot(x, aes(year, value))  + geom_line()
x <- get_rep(fit, 'recruit')
ggplot(x, aes(year, value))  + geom_line()
plot_pk_ssb(fit, add_uncertainty=0)
