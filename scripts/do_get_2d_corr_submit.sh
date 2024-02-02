#!/bin/bash

nx=200
ny=200
nseed=1
dim=2

for tau in 1.000000 10.000000 100.000000 1000.000000 10000.000000
do
    for lambda in 1.000000 3.000000 10.000000 30.000000
    do
        #seedfile="seeds/seeds_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}.txt"
        for (( i=1; i<=${nseed}; i++ ))
        do
            echo $i
            mydir="/scratch0/laynefrechette/active-noise-results/2d/nx\=${nx}_ny\=${ny}/tau\=${tau}/lambda\=${lambda}/seed\=$i/"
            sbatch -J "get_corr_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}" -o "log/get_corr_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}.o%j" -e "log/get_corr_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}_seed=${i}.e%j" scripts/submit_get_corr.sh "${mydir}noise_traj.h5" ${dim}
        done
    done
done
