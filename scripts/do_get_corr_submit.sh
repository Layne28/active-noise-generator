#!/bin/bash

nx=200
ny=200
nseed=1

for tau in 1.000000 10.000000 100.000000 1000.000000 10000.000000
do
    for lambda in 1.000000 3.000000 10.000000 30.000000
    do
        #seedfile="seeds/seeds_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}.txt"
        for (( i=1; i<=${nseed}; i++ ))
        do
            echo $i
            sbatch -J "noise_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}" -o "log/noise_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}.o%j" -e "log/noise_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}.e%j" scripts/submit_sims.sh "/home/laynefrechette/active-noise/active-noise-generator/input_files/active_noise_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}.in" 100000 $i #${seedfile}
        done
    done
done
