#!/bin/sh
# create standalone ifit distribution
#
# usage: make.sh <target>

# requires packages: realpath, Matlab Compiler

# location of the make script
export m=$PWD
# location of the iFit directory
export p=`realpath $PWD/../..`

# create the help pages
matlab -nojvm -nosplash -r "addpath(genpath('$p')); help_create; exit"

mv edit.m.org edit.m
mv web.m.org  web.m
mkdir -p $1
cd $1
echo Creating the stand-alone version 
echo from $p 
echo into $1
mcc -m ifit -a $p
mv ifit run_ifit
rm run_ifit.sh
cp $m/ifit .
cp $p/README.txt .
cp $p/COPYING .
cd $m
mv edit.m edit.m.org
mv web.m  web.m.org
