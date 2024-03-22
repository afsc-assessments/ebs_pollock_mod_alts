@echo off
:: if NOT EXIST arc mkdir arc

:: del amak.std
:: with model 0 in matrix
:: set /a modn=%1-1
:: without model 0 in matrix
set modn=_%1

:: awk "{print $%1}" models.dat >mod.ctl

:: amak -ind mod.ctl -binp arc\mod%modn%.bar -phase 22 -nox -iprint 100
amak -ind %1.ctl -nox -iprint 100

grep q_srv amak.std >>amak.rep
grep totbiom amak.std >>amak.rep
grep recruits amak.std >>amak.rep
grep Sp_Biom amak.std >>amak.rep
copy for_r.rep arc\mod%modn%_R.rep
copy mod.ctl arc\mod%modn%.ctl
copy amak.bar arc\mod%modn%.bar
copy amak.rep arc\mod%modn%.rep
copy amak.std arc\mod%modn%.std
:: copy extra_sd.rep arc\mod%1_extra_sd.rep

