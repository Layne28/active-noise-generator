#!/bin/bash

nx=200
ny=200
nseed=1

#taus=(1.000000 10.000000 100.000000 1000.000000 10000.000000)
#lambdas=(1.000000 3.000000 10.000000 30.000000)

taus=(1000.000000)
lambdas=(10.000000)

basedir=$PROJECT

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for (( i=1; i<=${nseed}; i++ ))
        do
            echo $i
            sbatch -J "noise_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}" -o "log/noise_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}.o%j" -e "log/noise_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}.e%j" scripts/submit_sims.sh "/jet/home/lfrechet/active-noise-generator/input_files/active_noise_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}.in" 100000 $i 
        done
    done
done
