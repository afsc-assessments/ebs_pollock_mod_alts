#!/bin/bash
cp arc/mod$1.prj proj/data/bsai_atka.dat
pushd proj
run bsai_atka
subl bsai_atka_out/*
popd
