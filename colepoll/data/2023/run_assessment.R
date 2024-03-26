## A quick script to reproduce the final 2023 assessment for GOA pollock

## devtools::install_github('afsc-assessments/GOApollock', ref='v0.1.2')
library(GOApollock)
library(TMB)

input <- prepare_pk_input(path=getwd(), datfile='pk23_10.txt', version='2023 final')
input <- prepare_pk_input(path=getwd(), datfile='gp.dat', version='2023 final')
str(input$dat)
str(input$par)

fit_pk
fit <- fit_pk(input)

str(fit$opt)
str(fit$rep)
str(fit$sd)

plot_data_overview(input)
read_pk_rep("goa.rep")
dat<-read_pk_dat("pk23_10.txt")
dat<-read_pk_dat("gp.dat")
data
