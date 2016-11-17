#!/bin/bash
rm -fr mod 
mkdir mod
dir='./src'
ifort -r8 -i8 -O3 -module mod $dir/global.F90  $dir/init.F90 $dir/main.F90 -o poisson_fftw -L/usr/local/fftw/intel/3.3.4/lib  -lfftw3
