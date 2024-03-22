@echo off
set modn=%1
amak -ind %1.ctl -nox -iprint 100

grep q_srv amak.std >>amak.rep
grep totbiom amak.std >>amak.rep
grep recruits amak.std >>amak.rep
grep Sp_Biom amak.std >>amak.rep
copy for_r.rep retro\%modn%_R.rep
copy mod.ctl   retro\%modn%.ctl
copy proj.dat   retro\%modn%.prj
copy amak.bar  retro\%modn%.bar
copy amak.rep  retro\%modn%.rep
copy amak.std  retro\%modn%.std
:: copy extra_sd.rep retro\mod%1_extra_sd.rep
call cleanad
