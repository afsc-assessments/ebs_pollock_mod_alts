#!/bin/bash
amak -nox -iprint 200 -ind mod$1.ctl -mcmc 1000000 -mcsave 200 
amak -nox -iprint 200 -ind mod$1.ctl -mceval
cp For_R.rep arc/mod_$1_R.rep
cp mceval.dat arc/mod_$1_mc.rep
cp amak.rep arc/mod_$1.rep
cp proj.dat arc/mod$1.prj
cp mod$1.ctl arc/mod$1.ctl
cp amak.std arc/mod$1.std
cp amak.par arc/mod$1.par
cleanad
