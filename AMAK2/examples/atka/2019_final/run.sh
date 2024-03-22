cp $1.ctl amak.dat
./amak -nox -iprint 200
cp amak.prj proj/data
cd proj
./main
cp bigfile.out ../arc/$1_proj.rep
cd ..
cp amak.rep arc/$1.rep
cp For_R.rep arc/$1_R.rep
cp cum_NLL.rep arc/$1_NLL.rep
cp amak.std arc/$1.std
cp amak.par arc/$1.par
cp amak.cor arc/$1.cor
