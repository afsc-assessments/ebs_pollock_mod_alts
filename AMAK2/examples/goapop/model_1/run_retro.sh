#!/bin/bash
# set modn=%1
# amak -ind ..\examples\atka\%1.ctl -nox -iprint 100
for i in `seq 0 15`;
do
  awk -v rrr=$i 'NR==12{print rrr} NR!=12 {print $0}' mod1.ctl >amak.dat
  amak -nox -iprint 200
  cp For_R.rep retro/r_$i.rep
	cp amak.std retro/r_$i.std
	cp amak.bar retro/r_$i.bar
done    
