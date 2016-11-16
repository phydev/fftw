#!/bin/bash
rm -fr mod 
mkdir mod
dir='./src'
ifort -r8 -mkl -O3 -module mod $dir/global.F90  $dir/init.F90 $dir/main.F90 -o poison_fftw
