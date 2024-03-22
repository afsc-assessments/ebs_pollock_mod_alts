#!/bin/bash
jjm -nox -iprint 200 -ind mod$1.ctl
cp For_R.rep arc/mod$1_R.rep
cp amak.rep arc/mod$1.rep
cp proj.dat arc/mod$1.prj
cp mod$1.ctl arc/mod$1.ctl
cp amak.std arc/mod$1.std
cp Fprof.yld arc/mod$1.yld
cp amak.par arc/mod$1.par
cleanad
