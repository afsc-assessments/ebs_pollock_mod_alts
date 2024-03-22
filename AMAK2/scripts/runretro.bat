@echo off
:: Args are first retro (0) and second to run retrospectives
FOR /L %%y IN (%1,1,%2) DO call rr.bat %3_r%%y

