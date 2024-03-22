:: batch job to do projections
cls
copy arc\%1.prj proj\data\BSAI_atka.dat
pushd proj
call run BSAI_Atka
cd BSAI_Atka_out
call vi percentiles.out
popd

