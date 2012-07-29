#!/bin/sh
# create standalone ifit distribution
#
# usage: make.sh <target>

echo usage: make.sh target

# location of the make script
export m=$PWD
# location of the iFit directory
export p=`realpath $PWD/../../..`
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
cd $m
mv edit.m edit.m.org
mv web.m  web.m.org
