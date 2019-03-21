#!/bin/bash
if [ "$1" != "" ]
then
scale=$1
else
scale=3
fi

if [ "$2" != "" ]
then
model=$2
else
model=2
fi

if [ "$3" != "" ]
then
k=$3
else
k=3
fi

#nohup ./analyze_fits_0.sh $scale $model $k 2>&1 > /dev/null &
#nohup ./analyze_fits_1.sh $scale $model $k 2>&1 > /dev/null &
#nohup ./analyze_fits_2.sh $scale $model $k 2>&1 > /dev/null &
#nohup ./analyze_fits_3.sh $scale $model $k 2>&1 > /dev/null &
#nohup ./analyze_fits_4.sh $scale $model $k 2>&1 > /dev/null &
nohup ./analyze_fits_5.sh $scale $model $k 2>&1 > /dev/null &
nohup ./analyze_fits_6.sh $scale $model $k 2>&1 > /dev/null &
nohup ./analyze_fits_7.sh $scale $model $k 2>&1 > /dev/null &
nohup ./analyze_fits_8.sh $scale $model $k 2>&1 > /dev/null &
nohup ./analyze_fits_9.sh $scale $model $k 2>&1 > /dev/null &
