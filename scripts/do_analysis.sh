#!/bin/bash

#$1 = nx, ny, nz

if [ "$1" == "" ]; then
    echo "Error: need system linear size"
    exit 0
fi

for tau in 0.010000 0.100000 1.000000 10.000000
do
    for lambda in 1.000000 2.000000 4.000000 8.000000
    do
        sbatch -J "analyze_noise_nx=$1_ny=$1_nz=$1_tau=${tau}_lambda=${lambda}" -o "analyze_noise_nx=$1_ny=$1_nz=$1_tau=${tau}_lambda=${lambda}.o%j" -e "analyze_noise_nx=$1_ny=$1_nz=$1_tau=${tau}_lambda=${lambda}.e%j" scripts/submit_analysis.sh "/home/laynefrechette/active-noise/active-noise-generator/conf/noise_analysis_nx=$1_ny=$1_nz=$1_tau=${tau}_lambda=${lambda}.conf"
    done
done
