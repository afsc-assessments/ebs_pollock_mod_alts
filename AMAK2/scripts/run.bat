@echo off
set modn=%1
amak -ind %1.ctl -nox -iprint 100

grep q_srv amak.std >>amak.rep
grep totbiom amak.std >>amak.rep
grep recruits amak.std >>amak.rep
grep Sp_Biom amak.std >>amak.rep
copy for_r.rep arc\%modn%_R.rep
copy mod.ctl   arc\%modn%.ctl
copy proj.dat   arc\%modn%.prj
copy amak.bar  arc\%modn%.bar
copy amak.par  arc\%modn%.par
copy amak.cor  arc\%modn%.cor
copy amak.rep  arc\%modn%.rep
copy amak.std  arc\%modn%.std
:: copy extra_sd.rep arc\mod%1_extra_sd.rep
call cleanad
